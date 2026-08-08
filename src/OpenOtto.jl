module OpenOtto

using DataFrames
using DifferentialEquations
using Plots

export EngineGeometry,
       GasProperties,
       Reservoir,
       Valve,
       WiebeLaw,
       EngineConfig,
       CylinderState,
       CycleResult,
       default_config,
       piston_position,
       cylinder_volume,
       dvolume_dangle,
       chamber_area,
       valve_area,
       wiebe_fraction,
       wiebe_derivative,
       orifice_mass_flow,
       simulate_cycle,
       simulate_to_periodic,
       plot_cycle,
       cycle_summary

const TWO_PI = 2π

"""Geometria de um cilindro. Todas as grandezas devem estar em SI."""
struct EngineGeometry
    bore::Float64
    stroke::Float64
    rod_length::Float64
    compression_ratio::Float64
    cylinders::Int

    function EngineGeometry(bore, stroke, rod_length, compression_ratio, cylinders=1)
        bore > 0 || throw(ArgumentError("bore deve ser positivo"))
        stroke > 0 || throw(ArgumentError("stroke deve ser positivo"))
        rod_length > stroke / 2 || throw(ArgumentError("rod_length deve ser maior que o raio da manivela"))
        compression_ratio > 1 || throw(ArgumentError("compression_ratio deve ser maior que 1"))
        cylinders >= 1 || throw(ArgumentError("cylinders deve ser pelo menos 1"))
        new(Float64(bore), Float64(stroke), Float64(rod_length), Float64(compression_ratio), Int(cylinders))
    end
end

"""Propriedades constantes do gás ideal usado na primeira abordagem."""
struct GasProperties
    R::Float64
    cp::Float64
    cv::Float64
    gamma::Float64

    function GasProperties(; R=287.0, gamma=1.35)
        R > 0 || throw(ArgumentError("R deve ser positivo"))
        gamma > 1 || throw(ArgumentError("gamma deve ser maior que 1"))
        cv = R / (gamma - 1)
        cp = gamma * cv
        new(Float64(R), Float64(cp), Float64(cv), Float64(gamma))
    end
end

"""Reservatório de pressão e temperatura constantes."""
struct Reservoir
    pressure::Float64
    temperature::Float64

    function Reservoir(pressure, temperature)
        pressure > 0 || throw(ArgumentError("pressure deve ser positiva"))
        temperature > 0 || throw(ArgumentError("temperature deve ser positiva"))
        new(Float64(pressure), Float64(temperature))
    end
end

"""
Válvula com janela angular e área efetiva senoidal.

Os ângulos são dados em graus no ciclo [-360, 360], com 0 grau no PMS de
combustão. A vazão usa `discharge_coefficient * area`.
"""
struct Valve
    opens_deg::Float64
    closes_deg::Float64
    max_area::Float64
    discharge_coefficient::Float64

    function Valve(opens_deg, closes_deg, max_area, discharge_coefficient=0.70)
        closes_deg > opens_deg || throw(ArgumentError("closes_deg deve ser maior que opens_deg"))
        closes_deg - opens_deg < 720 || throw(ArgumentError("a janela da válvula deve ser menor que 720 graus"))
        max_area >= 0 || throw(ArgumentError("max_area não pode ser negativa"))
        0 < discharge_coefficient <= 1 || throw(ArgumentError("discharge_coefficient deve pertencer a (0, 1]"))
        new(Float64(opens_deg), Float64(closes_deg), Float64(max_area), Float64(discharge_coefficient))
    end
end

"""Lei de Wiebe simples, normalizada para atingir exatamente a unidade."""
struct WiebeLaw
    start_deg::Float64
    duration_deg::Float64
    a::Float64
    shape::Float64

    function WiebeLaw(start_deg=-15.0, duration_deg=50.0, a=5.0, shape=2.0)
        duration_deg > 0 || throw(ArgumentError("duration_deg deve ser positiva"))
        a > 0 || throw(ArgumentError("a deve ser positivo"))
        shape > -1 || throw(ArgumentError("shape deve ser maior que -1"))
        new(Float64(start_deg), Float64(duration_deg), Float64(a), Float64(shape))
    end
end

"""Conjunto de parâmetros do modelo aberto zero-dimensional."""
struct EngineConfig
    geometry::EngineGeometry
    gas::GasProperties
    intake::Reservoir
    exhaust::Reservoir
    intake_valve::Valve
    exhaust_valve::Valve
    combustion::WiebeLaw
    rpm::Float64
    wall_temperature::Float64
    heat_transfer_coefficient::Float64
    equivalence_ratio::Float64
    stoichiometric_afr::Float64
    lower_heating_value::Float64
    combustion_efficiency::Float64

    function EngineConfig(
        geometry,
        gas,
        intake,
        exhaust,
        intake_valve,
        exhaust_valve,
        combustion;
        rpm=2000.0,
        wall_temperature=430.0,
        heat_transfer_coefficient=350.0,
        equivalence_ratio=1.0,
        stoichiometric_afr=14.7,
        lower_heating_value=44.0e6,
        combustion_efficiency=0.95,
    )
        rpm > 0 || throw(ArgumentError("rpm deve ser positiva"))
        wall_temperature > 0 || throw(ArgumentError("wall_temperature deve ser positiva"))
        heat_transfer_coefficient >= 0 || throw(ArgumentError("heat_transfer_coefficient não pode ser negativo"))
        equivalence_ratio >= 0 || throw(ArgumentError("equivalence_ratio não pode ser negativa"))
        stoichiometric_afr > 0 || throw(ArgumentError("stoichiometric_afr deve ser positiva"))
        lower_heating_value >= 0 || throw(ArgumentError("lower_heating_value não pode ser negativo"))
        0 <= combustion_efficiency <= 1 || throw(ArgumentError("combustion_efficiency deve pertencer a [0, 1]"))
        new(
            geometry,
            gas,
            intake,
            exhaust,
            intake_valve,
            exhaust_valve,
            combustion,
            Float64(rpm),
            Float64(wall_temperature),
            Float64(heat_transfer_coefficient),
            Float64(equivalence_ratio),
            Float64(stoichiometric_afr),
            Float64(lower_heating_value),
            Float64(combustion_efficiency),
        )
    end
end

"""Estado conservativo no início ou no fim de um ciclo."""
struct CylinderState
    mass::Float64
    internal_energy::Float64
end

"""Resultados tabulares, estado final e histórico de convergência."""
struct CycleResult
    data::DataFrame
    summary::DataFrame
    initial_state::CylinderState
    final_state::CylinderState
    cycles::Int
    convergence::DataFrame
end

crank_radius(g::EngineGeometry) = g.stroke / 2
piston_area(g::EngineGeometry) = π * g.bore^2 / 4
displacement(g::EngineGeometry) = piston_area(g) * g.stroke
clearance_volume(g::EngineGeometry) = displacement(g) / (g.compression_ratio - 1)
angular_speed(cfg::EngineConfig) = TWO_PI * cfg.rpm / 60

"""Deslocamento do pistão medido a partir do PMS (Equação 5 do TCC)."""
function piston_position(g::EngineGeometry, angle_rad)
    R = crank_radius(g)
    L = g.rod_length
    return R * (1 - cos(angle_rad)) + L - sqrt(L^2 - R^2 * sin(angle_rad)^2)
end

"""Volume instantâneo do cilindro (Equações 7 e 8 do TCC)."""
cylinder_volume(g::EngineGeometry, angle_rad) = clearance_volume(g) + piston_area(g) * piston_position(g, angle_rad)

"""Derivada analítica dV/dalpha, em m^3/rad."""
function dvolume_dangle(g::EngineGeometry, angle_rad)
    R = crank_radius(g)
    L = g.rod_length
    s = sin(angle_rad)
    c = cos(angle_rad)
    dx = R * s + R^2 * s * c / sqrt(L^2 - R^2 * s^2)
    return piston_area(g) * dx
end

"""Área molhada aproximada: cabeçote + pistão + camisa exposta."""
chamber_area(g::EngineGeometry, angle_rad) = 2 * piston_area(g) + π * g.bore * piston_position(g, angle_rad)

"""Área geométrica instantânea, com perfil sin^2 entre abertura e fechamento."""
function valve_area(v::Valve, angle_rad)
    angle_deg = rad2deg(angle_rad)
    # A janela é repetida a cada 720 graus. Isso representa corretamente uma
    # válvula que abre antes de -360 graus ou fecha depois de +360 graus.
    for cycle_shift in -1:1
        shifted_angle = angle_deg + 720 * cycle_shift
        if v.opens_deg < shifted_angle < v.closes_deg
            phase = (shifted_angle - v.opens_deg) / (v.closes_deg - v.opens_deg)
            return v.max_area * sin(π * phase)^2
        end
    end
    return 0.0
end

function wiebe_fraction(w::WiebeLaw, angle_rad)
    angle_deg = rad2deg(angle_rad)
    z = (angle_deg - w.start_deg) / w.duration_deg
    z <= 0 && return 0.0
    z >= 1 && return 1.0
    return (1 - exp(-w.a * z^(w.shape + 1))) / (1 - exp(-w.a))
end

"""Derivada da fração queimada normalizada em relação a alpha [1/rad]."""
function wiebe_derivative(w::WiebeLaw, angle_rad)
    angle_deg = rad2deg(angle_rad)
    z = (angle_deg - w.start_deg) / w.duration_deg
    (z <= 0 || z >= 1) && return 0.0
    duration_rad = deg2rad(w.duration_deg)
    numerator = w.a * (w.shape + 1) * z^w.shape * exp(-w.a * z^(w.shape + 1))
    return numerator / (duration_rad * (1 - exp(-w.a)))
end

"""
Vazão mássica compressível em um orifício, sempre positiva do montante para
jusante. A expressão alterna automaticamente entre os regimes bloqueado e não
bloqueado.
"""
function orifice_mass_flow(area, discharge_coefficient, p_up, T_up, p_down, gas::GasProperties)
    area <= 0 && return 0.0
    p_up <= p_down && return 0.0
    pressure_ratio = clamp(p_down / p_up, 0.0, 1.0)
    critical_ratio = (2 / (gas.gamma + 1))^(gas.gamma / (gas.gamma - 1))

    flow_factor = if pressure_ratio <= critical_ratio
        sqrt(gas.gamma * (2 / (gas.gamma + 1))^((gas.gamma + 1) / (gas.gamma - 1)))
    else
        term = 2 * gas.gamma / (gas.gamma - 1) *
               (pressure_ratio^(2 / gas.gamma) - pressure_ratio^((gas.gamma + 1) / gas.gamma))
        sqrt(max(term, 0.0))
    end

    return discharge_coefficient * area * p_up / sqrt(gas.R * T_up) * flow_factor
end

function signed_port_flow(reservoir::Reservoir, valve::Valve, angle_rad, p_cyl, T_cyl, gas)
    area = valve_area(valve, angle_rad)
    if reservoir.pressure >= p_cyl
        return orifice_mass_flow(area, valve.discharge_coefficient,
            reservoir.pressure, reservoir.temperature, p_cyl, gas)
    end
    return -orifice_mass_flow(area, valve.discharge_coefficient,
        p_cyl, T_cyl, reservoir.pressure, gas)
end

"""Calor total liberado por ciclo e por cilindro [J]."""
function nominal_heat_release(cfg::EngineConfig)
    g = cfg.geometry
    m_charge = cfg.intake.pressure * cylinder_volume(g, π) / (cfg.gas.R * cfg.intake.temperature)
    fuel_fraction = cfg.equivalence_ratio / (cfg.stoichiometric_afr + cfg.equivalence_ratio)
    return cfg.combustion_efficiency * cfg.lower_heating_value * fuel_fraction * m_charge
end

function cylinder_ode!(du, u, cfg::EngineConfig, angle_rad)
    mass = max(u[1], 1e-12)
    internal_energy = max(u[2], 1e-9)
    gas = cfg.gas
    g = cfg.geometry
    omega = angular_speed(cfg)

    volume = cylinder_volume(g, angle_rad)
    temperature = internal_energy / (mass * gas.cv)
    pressure = mass * gas.R * temperature / volume

    mdot_intake = signed_port_flow(cfg.intake, cfg.intake_valve, angle_rad, pressure, temperature, gas)
    mdot_exhaust = signed_port_flow(cfg.exhaust, cfg.exhaust_valve, angle_rad, pressure, temperature, gas)
    mdot_net = mdot_intake + mdot_exhaust

    hdot_intake = mdot_intake * gas.cp * (mdot_intake >= 0 ? cfg.intake.temperature : temperature)
    hdot_exhaust = mdot_exhaust * gas.cp * (mdot_exhaust >= 0 ? cfg.exhaust.temperature : temperature)
    hdot = hdot_intake + hdot_exhaust

    dV_dangle = dvolume_dangle(g, angle_rad)
    dQcomb_dangle = nominal_heat_release(cfg) * wiebe_derivative(cfg.combustion, angle_rad)
    qwall_dot = cfg.heat_transfer_coefficient * chamber_area(g, angle_rad) *
                (temperature - cfg.wall_temperature)

    du[1] = mdot_net / omega
    du[2] = dQcomb_dangle - qwall_dot / omega - pressure * dV_dangle + hdot / omega
    du[3] = pressure * dV_dangle
    du[4] = dQcomb_dangle
    du[5] = qwall_dot / omega
    du[6] = max(mdot_intake, 0.0) / omega
    du[7] = max(-mdot_exhaust, 0.0) / omega
    du[8] = hdot / omega
    return nothing
end

function default_config()
    # Mesmo ponto de partida geométrico do TCC: 1667 cm^3, 6 cilindros,
    # motor quadrado, L/R = 3.7 e taxa de compressão 13.
    total_displacement = 1667e-6
    cylinders = 6
    unit_displacement = total_displacement / cylinders
    bore = (4 * unit_displacement / π)^(1 / 3)
    stroke = bore
    crank = stroke / 2
    geometry = EngineGeometry(bore, stroke, 3.7 * crank, 13.0, cylinders)

    gas = GasProperties(R=287.0, gamma=1.35)
    intake = Reservoir(101.35e3, 303.15)
    exhaust = Reservoir(110.0e3, 750.0)
    intake_valve = Valve(-380.0, -130.0, 3.0e-4, 0.70)
    exhaust_valve = Valve(130.0, 380.0, 2.5e-4, 0.72)
    combustion = WiebeLaw(-15.0, 50.0, 5.0, 2.0)

    return EngineConfig(
        geometry,
        gas,
        intake,
        exhaust,
        intake_valve,
        exhaust_valve,
        combustion;
        rpm=2000.0,
        wall_temperature=430.0,
        heat_transfer_coefficient=350.0,
        equivalence_ratio=1.0,
        stoichiometric_afr=14.7,
        lower_heating_value=44.0e6,
        combustion_efficiency=0.95,
    )
end

function state_from_pressure_temperature(cfg::EngineConfig, pressure, temperature; angle_rad=-TWO_PI)
    volume = cylinder_volume(cfg.geometry, angle_rad)
    mass = pressure * volume / (cfg.gas.R * temperature)
    return CylinderState(mass, mass * cfg.gas.cv * temperature)
end

function state_temperature(cfg::EngineConfig, state::CylinderState)
    return state.internal_energy / (state.mass * cfg.gas.cv)
end

function integration_stops(cfg::EngineConfig)
    valve_events = [
        cfg.intake_valve.opens_deg,
        cfg.intake_valve.closes_deg,
        cfg.exhaust_valve.opens_deg,
        cfg.exhaust_valve.closes_deg,
    ]
    wrapped_valve_events = [event + 720 * shift for event in valve_events for shift in -1:1
                            if -360 < event + 720 * shift < 360]
    combustion_events = [
        cfg.combustion.start_deg,
        cfg.combustion.start_deg + cfg.combustion.duration_deg,
    ]
    return sort(unique(deg2rad.([wrapped_valve_events; combustion_events])))
end

function solve_one_cycle(cfg::EngineConfig, initial::CylinderState; step_deg=0.25, reltol=1e-8)
    u0 = [initial.mass, initial.internal_energy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    angle_span = (-TWO_PI, TWO_PI)
    save_angles = range(angle_span[1], angle_span[2]; step=deg2rad(step_deg))
    problem = ODEProblem(cylinder_ode!, u0, angle_span, cfg)
    return solve(
        problem,
        Tsit5();
        reltol=reltol,
        abstol=[1e-12, 1e-6, 1e-7, 1e-7, 1e-7, 1e-12, 1e-12, 1e-7],
        saveat=save_angles,
        tstops=integration_stops(cfg),
    )
end

function solution_dataframe(solution, cfg::EngineConfig)
    n = length(solution.t)
    angle_deg = Vector{Float64}(undef, n)
    volume_cm3 = Vector{Float64}(undef, n)
    pressure_kPa = Vector{Float64}(undef, n)
    temperature_K = Vector{Float64}(undef, n)
    mass_g = Vector{Float64}(undef, n)
    mdot_intake = Vector{Float64}(undef, n)
    mdot_exhaust = Vector{Float64}(undef, n)

    for i in eachindex(solution.t)
        angle = solution.t[i]
        u = solution.u[i]
        mass = u[1]
        temperature = u[2] / (mass * cfg.gas.cv)
        volume = cylinder_volume(cfg.geometry, angle)
        pressure = mass * cfg.gas.R * temperature / volume
        angle_deg[i] = rad2deg(angle)
        volume_cm3[i] = volume * 1e6
        pressure_kPa[i] = pressure / 1e3
        temperature_K[i] = temperature
        mass_g[i] = mass * 1e3
        mdot_intake[i] = signed_port_flow(cfg.intake, cfg.intake_valve, angle, pressure, temperature, cfg.gas)
        mdot_exhaust[i] = signed_port_flow(cfg.exhaust, cfg.exhaust_valve, angle, pressure, temperature, cfg.gas)
    end

    return DataFrame(
        angle_deg=angle_deg,
        volume_cm3=volume_cm3,
        pressure_kPa=pressure_kPa,
        temperature_K=temperature_K,
        mass_g=mass_g,
        work_J=[u[3] for u in solution.u],
        heat_combustion_J=[u[4] for u in solution.u],
        heat_wall_J=[u[5] for u in solution.u],
        intake_mass_g=[u[6] * 1e3 for u in solution.u],
        exhaust_mass_g=[u[7] * 1e3 for u in solution.u],
        enthalpy_flow_J=[u[8] for u in solution.u],
        mdot_intake_kg_s=mdot_intake,
        mdot_exhaust_kg_s=mdot_exhaust,
    )
end

function make_summary(data::DataFrame, cfg::EngineConfig, initial::CylinderState, final::CylinderState, cycles)
    work = data.work_J[end]
    qcomb = data.heat_combustion_J[end]
    qwall = data.heat_wall_J[end]
    hflow = data.enthalpy_flow_J[end]
    delta_u = final.internal_energy - initial.internal_energy
    energy_residual = delta_u - (qcomb - qwall - work + hflow)
    unit_displacement = displacement(cfg.geometry)
    cycle_frequency = cfg.rpm / 120
    indicated_power = work * cycle_frequency * cfg.geometry.cylinders
    reference_mass = cfg.intake.pressure * displacement(cfg.geometry) /
                     (cfg.gas.R * cfg.intake.temperature)

    return DataFrame(
        metric=[
            "cycles_to_convergence",
            "net_indicated_work_J_per_cylinder",
            "gross_heat_release_J_per_cylinder",
            "indicated_efficiency",
            "imep_kPa",
            "total_indicated_power_kW",
            "volumetric_efficiency",
            "peak_pressure_kPa",
            "peak_temperature_K",
            "wall_heat_loss_J_per_cylinder",
            "energy_residual_J",
        ],
        value=[
            Float64(cycles),
            work,
            qcomb,
            qcomb == 0 ? NaN : work / qcomb,
            work / unit_displacement / 1e3,
            indicated_power / 1e3,
            data.intake_mass_g[end] / (reference_mass * 1e3),
            maximum(data.pressure_kPa),
            maximum(data.temperature_K),
            qwall,
            energy_residual,
        ],
        unit=["-", "J", "J", "-", "kPa", "kW", "-", "kPa", "K", "J", "J"],
    )
end

"""Executa exatamente um ciclo de 720 graus."""
function simulate_cycle(
    cfg::EngineConfig=default_config();
    initial_state=nothing,
    initial_pressure=cfg.exhaust.pressure,
    initial_temperature=750.0,
    step_deg=0.25,
    reltol=1e-8,
)
    initial = isnothing(initial_state) ?
        state_from_pressure_temperature(cfg, initial_pressure, initial_temperature) : initial_state
    solution = solve_one_cycle(cfg, initial; step_deg=step_deg, reltol=reltol)
    final = CylinderState(solution.u[end][1], solution.u[end][2])
    data = solution_dataframe(solution, cfg)
    convergence = DataFrame(cycle=[1], mass_g=[final.mass * 1e3],
        temperature_K=[state_temperature(cfg, final)], relative_change=[NaN])
    summary = make_summary(data, cfg, initial, final, 1)
    return CycleResult(data, summary, initial, final, 1, convergence)
end

"""Repete ciclos até que massa e temperatura no início/fim sejam periódicas."""
function simulate_to_periodic(
    cfg::EngineConfig=default_config();
    initial_pressure=cfg.exhaust.pressure,
    initial_temperature=750.0,
    step_deg=0.25,
    reltol=1e-8,
    periodic_tolerance=1e-6,
    min_cycles=3,
    max_cycles=100,
)
    state = state_from_pressure_temperature(cfg, initial_pressure, initial_temperature)
    history = DataFrame(cycle=Int[], mass_g=Float64[], temperature_K=Float64[], relative_change=Float64[])
    final_solution = nothing
    cycle_initial = state

    for cycle in 1:max_cycles
        cycle_initial = state
        solution = solve_one_cycle(cfg, cycle_initial; step_deg=step_deg, reltol=reltol)
        state = CylinderState(solution.u[end][1], solution.u[end][2])
        old_temperature = state_temperature(cfg, cycle_initial)
        new_temperature = state_temperature(cfg, state)
        relative_change = max(
            abs(state.mass - cycle_initial.mass) / max(state.mass, 1e-12),
            abs(new_temperature - old_temperature) / max(new_temperature, 1e-12),
        )
        push!(history, (cycle, state.mass * 1e3, new_temperature, relative_change))
        final_solution = solution

        if cycle >= min_cycles && relative_change <= periodic_tolerance
            data = solution_dataframe(solution, cfg)
            summary = make_summary(data, cfg, cycle_initial, state, cycle)
            return CycleResult(data, summary, cycle_initial, state, cycle, history)
        end
    end

    @warn "O estado periódico não foi atingido" max_cycles periodic_tolerance last_change=history.relative_change[end]
    data = solution_dataframe(final_solution, cfg)
    summary = make_summary(data, cfg, cycle_initial, state, max_cycles)
    return CycleResult(data, summary, cycle_initial, state, max_cycles, history)
end

cycle_summary(result::CycleResult) = result.summary

"""Gera painel com p-alpha, T-alpha, p-V e vazões nas válvulas."""
function plot_cycle(result::CycleResult; output_path=nothing)
    d = result.data
    p1 = plot(d.angle_deg, d.pressure_kPa; xlabel="Ângulo [graus]", ylabel="Pressão [kPa]",
        label=false, linewidth=2, title="Pressão no cilindro")
    p2 = plot(d.angle_deg, d.temperature_K; xlabel="Ângulo [graus]", ylabel="Temperatura [K]",
        label=false, linewidth=2, title="Temperatura no cilindro")
    p3 = plot(d.volume_cm3, d.pressure_kPa; xlabel="Volume [cm³]", ylabel="Pressão [kPa]",
        label=false, linewidth=2, title="Diagrama p-V")
    p4 = plot(d.angle_deg, 1e3 .* d.mdot_intake_kg_s; xlabel="Ângulo [graus]",
        ylabel="Vazão [g/s]", label="Admissão", linewidth=2, title="Vazões nas válvulas")
    plot!(p4, d.angle_deg, 1e3 .* d.mdot_exhaust_kg_s; label="Escape", linewidth=2)
    panel = plot(p1, p2, p3, p4; layout=(2, 2), size=(1100, 760))
    !isnothing(output_path) && savefig(panel, output_path)
    return panel
end

# Solver avançado: composição reativa, propriedades cp(T) e atrito.
include("AdvancedModel.jl")

# Baterias paramétricas e gráficos comparativos.
include("ParametricStudies.jl")

end

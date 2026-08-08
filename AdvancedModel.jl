export PhysicsOptions,
       NASAPolynomial,
       ReactiveFuel,
       AdvancedEngineConfig,
       AdvancedState,
       AdvancedCycleResult,
       methane_fuel,
       default_advanced_config,
       with_physics,
       reconfigure_advanced,
       ideal_physics,
       reaction_stoichiometry,
       reaction_heat_per_mol,
       intake_mole_fractions,
       exhaust_mole_fractions,
       species_cp,
       species_enthalpy,
       species_molar_mass,
       mixture_mass,
       mixture_internal_energy,
       mixture_temperature,
       advanced_initial_state,
       simulate_advanced_cycle,
       simulate_advanced_to_periodic,
       plot_advanced_cycle,
       advanced_cycle_summary

const RU = 8.31446261815324
const T_REF = 298.15
const ADV_SPECIES = (:fuel, :O2, :N2, :CO2, :H2O, :CO, :H2)
const NSPECIES = length(ADV_SPECIES)

const I_FUEL = 1
const I_O2 = 2
const I_N2 = 3
const I_CO2 = 4
const I_H2O = 5
const I_CO = 6
const I_H2 = 7
const I_U = 8
const I_WIND = 9
const I_WFRIC = 10
const I_QRELEASE = 11
const I_QWALL = 12
const I_HFLOW = 13
const I_MIN = 14
const I_MOUT = 15
const I_BURNABLE_FUEL = 16
const ADV_NSTATES = 16

"""Polinômio NASA de sete coeficientes, com transição em `T_switch`."""
struct NASAPolynomial
    low::NTuple{7, Float64}
    high::NTuple{7, Float64}
    T_switch::Float64
end

function NASAPolynomial(low, high, T_switch=1000.0)
    length(low) == 7 || throw(ArgumentError("low deve possuir sete coeficientes"))
    length(high) == 7 || throw(ArgumentError("high deve possuir sete coeficientes"))
    return NASAPolynomial(Tuple(Float64.(low)), Tuple(Float64.(high)), Float64(T_switch))
end

"""Combustível genérico CcHhOoNn com propriedades NASA fornecidas."""
struct ReactiveFuel
    name::Symbol
    carbon::Float64
    hydrogen::Float64
    oxygen::Float64
    nitrogen::Float64
    molar_mass::Float64
    thermo::NASAPolynomial

    function ReactiveFuel(name, carbon, hydrogen, oxygen, nitrogen, molar_mass, thermo)
        carbon >= 0 || throw(ArgumentError("carbon não pode ser negativo"))
        hydrogen >= 0 || throw(ArgumentError("hydrogen não pode ser negativo"))
        molar_mass > 0 || throw(ArgumentError("molar_mass deve ser positiva"))
        new(Symbol(name), Float64(carbon), Float64(hydrogen), Float64(oxygen),
            Float64(nitrogen), Float64(molar_mass), thermo)
    end
end

"""
Seleção modular das hipóteses físicas.

- `combustion_model`: `:isochoric`, `:wiebe_source` ou `:tcc_air_fuel`;
- `friction_model`: `:none` ou `:linear_velocity`;
- `heat_transfer_model`: `:adiabatic` ou `:constant_h`;
- `gas_exchange_model`: `:closed` ou `:compressible_valves`;
- `property_model`: `:constant` ou `:nasa`.
"""
struct PhysicsOptions
    combustion_model::Symbol
    friction_model::Symbol
    heat_transfer_model::Symbol
    gas_exchange_model::Symbol
    property_model::Symbol

    function PhysicsOptions(;
        combustion_model=:tcc_air_fuel,
        friction_model=:linear_velocity,
        heat_transfer_model=:constant_h,
        gas_exchange_model=:compressible_valves,
        property_model=:nasa,
    )
        combustion_model in (:isochoric, :wiebe_source, :tcc_air_fuel) ||
            throw(ArgumentError("combustion_model inválido"))
        friction_model in (:none, :linear_velocity) ||
            throw(ArgumentError("friction_model inválido"))
        heat_transfer_model in (:adiabatic, :constant_h) ||
            throw(ArgumentError("heat_transfer_model inválido"))
        gas_exchange_model in (:closed, :compressible_valves) ||
            throw(ArgumentError("gas_exchange_model inválido"))
        property_model in (:constant, :nasa) ||
            throw(ArgumentError("property_model inválido"))
        new(combustion_model, friction_model, heat_transfer_model,
            gas_exchange_model, property_model)
    end
end

ideal_physics() = PhysicsOptions(
    combustion_model=:isochoric,
    friction_model=:none,
    heat_transfer_model=:adiabatic,
    gas_exchange_model=:closed,
    property_model=:constant,
)

"""Configuração do solver avançado."""
struct AdvancedEngineConfig
    geometry::EngineGeometry
    intake::Reservoir
    exhaust::Reservoir
    intake_valve::Valve
    exhaust_valve::Valve
    combustion::WiebeLaw
    fuel::ReactiveFuel
    physics::PhysicsOptions
    rpm::Float64
    equivalence_ratio::Float64
    wall_temperature::Float64
    heat_transfer_coefficient::Float64
    friction_coefficient::Float64
    atmospheric_n2_o2_ratio::Float64

    function AdvancedEngineConfig(
        geometry,
        intake,
        exhaust,
        intake_valve,
        exhaust_valve,
        combustion,
        fuel,
        physics;
        rpm=2000.0,
        equivalence_ratio=1.0,
        wall_temperature=430.0,
        heat_transfer_coefficient=350.0,
        friction_coefficient=16.0,
        atmospheric_n2_o2_ratio=3.76,
    )
        rpm > 0 || throw(ArgumentError("rpm deve ser positiva"))
        0 < equivalence_ratio <= 1 ||
            throw(ArgumentError("esta versão do modelo TCC requer 0 < equivalence_ratio <= 1"))
        wall_temperature > 0 || throw(ArgumentError("wall_temperature deve ser positiva"))
        heat_transfer_coefficient >= 0 ||
            throw(ArgumentError("heat_transfer_coefficient não pode ser negativo"))
        friction_coefficient >= 0 ||
            throw(ArgumentError("friction_coefficient não pode ser negativo"))
        atmospheric_n2_o2_ratio > 0 ||
            throw(ArgumentError("atmospheric_n2_o2_ratio deve ser positiva"))
        new(geometry, intake, exhaust, intake_valve, exhaust_valve, combustion,
            fuel, physics, Float64(rpm), Float64(equivalence_ratio),
            Float64(wall_temperature), Float64(heat_transfer_coefficient),
            Float64(friction_coefficient), Float64(atmospheric_n2_o2_ratio))
    end
end

struct AdvancedState
    moles::Vector{Float64}
    internal_energy::Float64

    function AdvancedState(moles, internal_energy)
        length(moles) == NSPECIES ||
            throw(ArgumentError("moles deve ter $(NSPECIES) componentes"))
        new(Float64.(moles), Float64(internal_energy))
    end
end

struct AdvancedCycleResult
    data::DataFrame
    summary::DataFrame
    initial_state::AdvancedState
    final_state::AdvancedState
    cycles::Int
    convergence::DataFrame
end

# Coeficientes do banco termoquímico GRI-Mech 3.0.
const NASA_O2 = NASAPolynomial(
    (3.78245636, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9,
     3.24372837e-12, -1.06394356e3, 3.65767573),
    (3.28253784, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10,
     -2.16717794e-14, -1.08845772e3, 5.45323129),
)
const NASA_N2 = NASAPolynomial(
    (3.53100528, -1.23660987e-4, -5.02999433e-7, 2.43530612e-9,
     -1.40881235e-12, -1.04697628e3, 2.96747468),
    (2.95257626, 1.39690040e-3, -4.92631603e-7, 7.86010195e-11,
     -4.60755204e-15, -9.23948645e2, 5.87188762),
)
const NASA_CO2 = NASAPolynomial(
    (2.35677352, 8.98459677e-3, -7.12356269e-6, 2.45919022e-9,
     -1.43699548e-13, -4.83719697e4, 9.90105222),
    (3.85796028, 4.41437026e-3, -2.21481404e-6, 5.23490188e-10,
     -4.72084164e-14, -4.87591660e4, 2.27163806),
)
const NASA_H2O = NASAPolynomial(
    (4.19864056, -2.03643410e-3, 6.52040211e-6, -5.48797062e-9,
     1.77197817e-12, -3.02937267e4, -0.849032208),
    (3.03399249, 2.17691804e-3, -1.64072518e-7, -9.70419870e-11,
     1.68200992e-14, -3.00042971e4, 4.96677010),
)
const NASA_CO = NASAPolynomial(
    (3.57953347, -6.10353680e-4, 1.01681433e-6, 9.07005884e-10,
     -9.04424499e-13, -1.43440860e4, 3.50840928),
    (2.71518561, 2.06252743e-3, -9.98825771e-7, 2.30053008e-10,
     -2.03647716e-14, -1.41518724e4, 7.81868772),
)
const NASA_H2 = NASAPolynomial(
    (2.34433112, 7.98052075e-3, -1.94781510e-5, 2.01572094e-8,
     -7.37611761e-12, -9.17935173e2, 0.683010238),
    (3.33727920, -4.94024731e-5, 4.99456778e-7, -1.79566394e-10,
     2.00255376e-14, -9.50158922e2, -3.20502331),
)
const NASA_CH4 = NASAPolynomial(
    (5.14987613, -1.36709788e-2, 4.91800599e-5, -4.84743026e-8,
     1.66693956e-11, -1.02466476e4, -4.64130376),
    (7.48514950e-2, 1.33909467e-2, -5.73285809e-6, 1.22292535e-9,
     -1.01815230e-13, -9.46834459e3, 18.4373180),
)

const ADV_FIXED_THERMO = (NASA_O2, NASA_N2, NASA_CO2, NASA_H2O, NASA_CO, NASA_H2)
const ADV_FIXED_MOLAR_MASS = (31.9988e-3, 28.0134e-3, 44.0095e-3,
    18.01528e-3, 28.0101e-3, 2.01588e-3)

methane_fuel() = ReactiveFuel(:CH4, 1, 4, 0, 0, 16.04246e-3, NASA_CH4)

function default_advanced_config(;
    physics=PhysicsOptions(),
    rpm=nothing,
    equivalence_ratio=nothing,
    wall_temperature=nothing,
    heat_transfer_coefficient=nothing,
    friction_coefficient=16.0,
)
    base = default_config()
    return AdvancedEngineConfig(
        base.geometry,
        base.intake,
        base.exhaust,
        base.intake_valve,
        base.exhaust_valve,
        base.combustion,
        methane_fuel(),
        physics;
        rpm=isnothing(rpm) ? base.rpm : rpm,
        equivalence_ratio=isnothing(equivalence_ratio) ? base.equivalence_ratio : equivalence_ratio,
        wall_temperature=isnothing(wall_temperature) ? base.wall_temperature : wall_temperature,
        heat_transfer_coefficient=isnothing(heat_transfer_coefficient) ?
            base.heat_transfer_coefficient : heat_transfer_coefficient,
        friction_coefficient=friction_coefficient,
        atmospheric_n2_o2_ratio=3.76,
    )
end

function with_physics(cfg::AdvancedEngineConfig, physics::PhysicsOptions)
    return AdvancedEngineConfig(
        cfg.geometry, cfg.intake, cfg.exhaust, cfg.intake_valve,
        cfg.exhaust_valve, cfg.combustion, cfg.fuel, physics;
        rpm=cfg.rpm,
        equivalence_ratio=cfg.equivalence_ratio,
        wall_temperature=cfg.wall_temperature,
        heat_transfer_coefficient=cfg.heat_transfer_coefficient,
        friction_coefficient=cfg.friction_coefficient,
        atmospheric_n2_o2_ratio=cfg.atmospheric_n2_o2_ratio,
    )
end

"""Cria uma configuração avançada alterando apenas os parâmetros indicados."""
function reconfigure_advanced(cfg::AdvancedEngineConfig;
    rpm=cfg.rpm,
    equivalence_ratio=cfg.equivalence_ratio,
    compression_ratio=cfg.geometry.compression_ratio,
    property_model=cfg.physics.property_model,
)
    geometry = EngineGeometry(
        cfg.geometry.bore,
        cfg.geometry.stroke,
        cfg.geometry.rod_length,
        compression_ratio,
        cfg.geometry.cylinders,
    )
    physics = PhysicsOptions(
        combustion_model=cfg.physics.combustion_model,
        friction_model=cfg.physics.friction_model,
        heat_transfer_model=cfg.physics.heat_transfer_model,
        gas_exchange_model=cfg.physics.gas_exchange_model,
        property_model=property_model,
    )
    return AdvancedEngineConfig(
        geometry, cfg.intake, cfg.exhaust, cfg.intake_valve,
        cfg.exhaust_valve, cfg.combustion, cfg.fuel, physics;
        rpm=rpm,
        equivalence_ratio=equivalence_ratio,
        wall_temperature=cfg.wall_temperature,
        heat_transfer_coefficient=cfg.heat_transfer_coefficient,
        friction_coefficient=cfg.friction_coefficient,
        atmospheric_n2_o2_ratio=cfg.atmospheric_n2_o2_ratio,
    )
end

advanced_angular_speed(cfg::AdvancedEngineConfig) = TWO_PI * cfg.rpm / 60

function species_thermo(cfg::AdvancedEngineConfig, index::Int)
    index == I_FUEL && return cfg.fuel.thermo
    return ADV_FIXED_THERMO[index - 1]
end

function species_molar_mass(cfg::AdvancedEngineConfig, index::Int)
    index == I_FUEL && return cfg.fuel.molar_mass
    return ADV_FIXED_MOLAR_MASS[index - 1]
end

function nasa_coefficients(model::NASAPolynomial, temperature)
    return temperature <= model.T_switch ? model.low : model.high
end

function nasa_cp(model::NASAPolynomial, temperature)
    T = Float64(temperature)
    a = nasa_coefficients(model, T)
    return RU * (a[1] + a[2] * T + a[3] * T^2 + a[4] * T^3 + a[5] * T^4)
end

function nasa_enthalpy(model::NASAPolynomial, temperature)
    T = Float64(temperature)
    a = nasa_coefficients(model, T)
    return RU * T * (a[1] + a[2] * T / 2 + a[3] * T^2 / 3 +
        a[4] * T^3 / 4 + a[5] * T^4 / 5 + a[6] / T)
end

function species_cp(cfg::AdvancedEngineConfig, index::Int, temperature)
    model = species_thermo(cfg, index)
    cfg.physics.property_model == :nasa && return nasa_cp(model, temperature)
    return nasa_cp(model, T_REF)
end

function species_enthalpy(cfg::AdvancedEngineConfig, index::Int, temperature)
    model = species_thermo(cfg, index)
    cfg.physics.property_model == :nasa && return nasa_enthalpy(model, temperature)
    href = nasa_enthalpy(model, T_REF)
    return href + nasa_cp(model, T_REF) * (temperature - T_REF)
end

species_internal_energy(cfg, index, temperature) =
    species_enthalpy(cfg, index, temperature) - RU * temperature

function mixture_internal_energy(cfg::AdvancedEngineConfig, moles, temperature)
    return sum(max(moles[i], 0.0) * species_internal_energy(cfg, i, temperature)
        for i in 1:NSPECIES)
end

function mixture_enthalpy(cfg::AdvancedEngineConfig, mole_fractions, temperature)
    return sum(mole_fractions[i] * species_enthalpy(cfg, i, temperature)
        for i in 1:NSPECIES)
end

function mixture_molar_mass(cfg::AdvancedEngineConfig, mole_fractions)
    return sum(mole_fractions[i] * species_molar_mass(cfg, i) for i in 1:NSPECIES)
end

function mixture_mass(cfg::AdvancedEngineConfig, moles)
    return sum(max(moles[i], 0.0) * species_molar_mass(cfg, i) for i in 1:NSPECIES)
end

function normalized_mole_fractions(moles)
    positive_moles = max.(moles, 0.0)
    total = sum(positive_moles)
    total > 0 || throw(DomainError(total, "a mistura não possui matéria"))
    return positive_moles ./ total
end

function mixture_heat_capacities(cfg::AdvancedEngineConfig, moles, temperature)
    x = normalized_mole_fractions(moles)
    cp_molar = sum(x[i] * species_cp(cfg, i, temperature) for i in 1:NSPECIES)
    cv_molar = cp_molar - RU
    return cp_molar, cv_molar, cp_molar / cv_molar
end

"""Inverte U(T,n) por bisseção monotônica."""
function mixture_temperature(cfg::AdvancedEngineConfig, moles, internal_energy)
    lower = 150.0
    upper = 6500.0
    f_lower = mixture_internal_energy(cfg, moles, lower) - internal_energy
    f_upper = mixture_internal_energy(cfg, moles, upper) - internal_energy
    f_lower <= 0 || throw(DomainError(internal_energy, "energia abaixo do limite térmico"))
    f_upper >= 0 || throw(DomainError(internal_energy, "energia acima do limite térmico"))

    for _ in 1:48
        middle = (lower + upper) / 2
        residual = mixture_internal_energy(cfg, moles, middle) - internal_energy
        if residual > 0
            upper = middle
        else
            lower = middle
        end
    end
    return (lower + upper) / 2
end

function reaction_stoichiometry(cfg::AdvancedEngineConfig)
    f = cfg.fuel
    oxygen_requirement = f.carbon + f.hydrogen / 4 - f.oxygen / 2
    oxygen_requirement > 0 || throw(ArgumentError("combustível sem demanda positiva de O2"))
    return [-1.0, -oxygen_requirement, f.nitrogen / 2,
        f.carbon, f.hydrogen / 2, 0.0, 0.0]
end

function reaction_heat_per_mol(cfg::AdvancedEngineConfig)
    nu = reaction_stoichiometry(cfg)
    delta_h = sum(nu[i] * nasa_enthalpy(species_thermo(cfg, i), T_REF)
        for i in 1:NSPECIES)
    return -delta_h
end

function intake_mole_fractions(cfg::AdvancedEngineConfig)
    f = cfg.fuel
    stoich_o2 = f.carbon + f.hydrogen / 4 - f.oxygen / 2
    oxygen = stoich_o2 / cfg.equivalence_ratio
    amounts = zeros(NSPECIES)
    amounts[I_FUEL] = 1.0
    amounts[I_O2] = oxygen
    amounts[I_N2] = cfg.atmospheric_n2_o2_ratio * oxygen
    return amounts ./ sum(amounts)
end

function exhaust_mole_fractions(cfg::AdvancedEngineConfig)
    f = cfg.fuel
    stoich_o2 = f.carbon + f.hydrogen / 4 - f.oxygen / 2
    supplied_o2 = stoich_o2 / cfg.equivalence_ratio
    amounts = zeros(NSPECIES)
    amounts[I_CO2] = f.carbon
    amounts[I_H2O] = f.hydrogen / 2
    amounts[I_O2] = supplied_o2 - stoich_o2
    amounts[I_N2] = cfg.atmospheric_n2_o2_ratio * supplied_o2 + f.nitrogen / 2
    return amounts ./ sum(amounts)
end

function mixture_flow_properties(cfg::AdvancedEngineConfig, composition, temperature)
    cp_molar = sum(composition[i] * species_cp(cfg, i, temperature) for i in 1:NSPECIES)
    cv_molar = cp_molar - RU
    molar_mass = mixture_molar_mass(cfg, composition)
    return RU / molar_mass, cp_molar / cv_molar, molar_mass
end

function advanced_orifice_mass_flow(area, discharge_coefficient, p_up, T_up, p_down,
    R_specific, gamma)
    area <= 0 && return 0.0
    p_up <= p_down && return 0.0
    pressure_ratio = clamp(p_down / p_up, 0.0, 1.0)
    critical_ratio = (2 / (gamma + 1))^(gamma / (gamma - 1))
    factor = if pressure_ratio <= critical_ratio
        sqrt(gamma * (2 / (gamma + 1))^((gamma + 1) / (gamma - 1)))
    else
        term = 2 * gamma / (gamma - 1) *
            (pressure_ratio^(2 / gamma) - pressure_ratio^((gamma + 1) / gamma))
        sqrt(max(term, 0.0))
    end
    return discharge_coefficient * area * p_up / sqrt(R_specific * T_up) * factor
end

function advanced_port_flow(cfg::AdvancedEngineConfig, reservoir::Reservoir, valve::Valve,
    reservoir_composition, angle, cylinder_pressure, cylinder_temperature, cylinder_moles)
    area = valve_area(valve, angle)
    cylinder_composition = normalized_mole_fractions(cylinder_moles)

    if reservoir.pressure >= cylinder_pressure
        upstream_composition = reservoir_composition
        upstream_temperature = reservoir.temperature
        p_up = reservoir.pressure
        p_down = cylinder_pressure
        direction = 1.0
    else
        upstream_composition = cylinder_composition
        upstream_temperature = cylinder_temperature
        p_up = cylinder_pressure
        p_down = reservoir.pressure
        direction = -1.0
    end

    R_specific, gamma, molar_mass = mixture_flow_properties(
        cfg, upstream_composition, upstream_temperature)
    mass_rate = advanced_orifice_mass_flow(area, valve.discharge_coefficient,
        p_up, upstream_temperature, p_down, R_specific, gamma)
    molar_rate = mass_rate / molar_mass
    species_rates = direction * molar_rate .* upstream_composition
    enthalpy_rate = direction * mass_rate / molar_mass *
        mixture_enthalpy(cfg, upstream_composition, upstream_temperature)
    return species_rates, direction * mass_rate, enthalpy_rate
end

function nominal_external_heat(cfg::AdvancedEngineConfig)
    volume_bdc = cylinder_volume(cfg.geometry, π)
    total_moles = cfg.intake.pressure * volume_bdc / (RU * cfg.intake.temperature)
    fuel_moles = total_moles * intake_mole_fractions(cfg)[I_FUEL]
    return reaction_heat_per_mol(cfg) * fuel_moles
end

function advanced_initial_state(cfg::AdvancedEngineConfig;
    pressure=nothing, temperature=nothing)
    closed = cfg.physics.gas_exchange_model == :closed
    initial_pressure = isnothing(pressure) ?
        (closed ? cfg.intake.pressure : cfg.exhaust.pressure) : pressure
    initial_temperature = isnothing(temperature) ?
        (closed ? cfg.intake.temperature : cfg.exhaust.temperature) : temperature
    initial_angle = closed ? -π : -TWO_PI
    composition = closed ? intake_mole_fractions(cfg) : exhaust_mole_fractions(cfg)
    volume = cylinder_volume(cfg.geometry, initial_angle)
    total_moles = initial_pressure * volume / (RU * initial_temperature)
    moles = total_moles .* composition
    internal_energy = mixture_internal_energy(cfg, moles, initial_temperature)
    return AdvancedState(moles, internal_energy)
end

function advanced_rhs!(du, u, cfg::AdvancedEngineConfig, angle)
    fill!(du, 0.0)
    moles = view(u, 1:NSPECIES)
    temperature = mixture_temperature(cfg, moles, u[I_U])
    volume = cylinder_volume(cfg.geometry, angle)
    total_moles = sum(max.(moles, 0.0))
    pressure = total_moles * RU * temperature / volume
    omega = advanced_angular_speed(cfg)

    intake_species_rate = zeros(NSPECIES)
    exhaust_species_rate = zeros(NSPECIES)
    intake_mass_rate = 0.0
    exhaust_mass_rate = 0.0
    enthalpy_rate = 0.0

    if cfg.physics.gas_exchange_model == :compressible_valves
        intake_species_rate, intake_mass_rate, intake_enthalpy_rate = advanced_port_flow(
            cfg, cfg.intake, cfg.intake_valve, intake_mole_fractions(cfg), angle,
            pressure, temperature, moles)
        exhaust_species_rate, exhaust_mass_rate, exhaust_enthalpy_rate = advanced_port_flow(
            cfg, cfg.exhaust, cfg.exhaust_valve, exhaust_mole_fractions(cfg), angle,
            pressure, temperature, moles)
        enthalpy_rate = intake_enthalpy_rate + exhaust_enthalpy_rate
    end

    for i in 1:NSPECIES
        du[i] = (intake_species_rate[i] + exhaust_species_rate[i]) / omega
    end

    released_energy_rate = 0.0
    external_heat_rate = 0.0
    if cfg.physics.combustion_model == :tcc_air_fuel
        burn_derivative = wiebe_derivative(cfg.combustion, angle)
        if burn_derivative > 0 && u[I_BURNABLE_FUEL] > 0
            extent_rate = u[I_BURNABLE_FUEL] * burn_derivative
            nu = reaction_stoichiometry(cfg)
            for i in 1:NSPECIES
                du[i] += nu[i] * extent_rate
            end
            released_energy_rate = reaction_heat_per_mol(cfg) * extent_rate
        end
    elseif cfg.physics.combustion_model == :wiebe_source
        external_heat_rate = nominal_external_heat(cfg) * wiebe_derivative(cfg.combustion, angle)
        released_energy_rate = external_heat_rate
    end

    dV_dangle = dvolume_dangle(cfg.geometry, angle)
    indicated_work_rate = pressure * dV_dangle
    heat_wall_rate = if cfg.physics.heat_transfer_model == :constant_h
        cfg.heat_transfer_coefficient * chamber_area(cfg.geometry, angle) *
            (temperature - cfg.wall_temperature) / omega
    else
        0.0
    end

    dx_dangle = dV_dangle / piston_area(cfg.geometry)
    friction_work_rate = if cfg.physics.friction_model == :linear_velocity
        cfg.friction_coefficient * omega * dx_dangle^2
    else
        0.0
    end

    du[I_U] = external_heat_rate - heat_wall_rate - indicated_work_rate + enthalpy_rate / omega
    du[I_WIND] = indicated_work_rate
    du[I_WFRIC] = friction_work_rate
    du[I_QRELEASE] = released_energy_rate
    du[I_QWALL] = heat_wall_rate
    du[I_HFLOW] = enthalpy_rate / omega
    du[I_MIN] = max(intake_mass_rate, 0.0) / omega
    du[I_MOUT] = max(-exhaust_mass_rate, 0.0) / omega
    return nothing
end

function advanced_integration_stops(cfg::AdvancedEngineConfig)
    valve_events = [cfg.intake_valve.opens_deg, cfg.intake_valve.closes_deg,
        cfg.exhaust_valve.opens_deg, cfg.exhaust_valve.closes_deg]
    wrapped = [event + 720 * shift for event in valve_events for shift in -1:1
        if -360 < event + 720 * shift < 360]
    combustion_events = [cfg.combustion.start_deg,
        cfg.combustion.start_deg + cfg.combustion.duration_deg, 0.0]
    return sort(unique(deg2rad.([wrapped; combustion_events])))
end

function solve_advanced_cycle(cfg::AdvancedEngineConfig, initial::AdvancedState;
    step_deg=0.25, reltol=1e-8)
    u0 = zeros(ADV_NSTATES)
    u0[1:NSPECIES] .= initial.moles
    u0[I_U] = initial.internal_energy
    angle_span = cfg.physics.gas_exchange_model == :closed ? (-π, π) : (-TWO_PI, TWO_PI)
    save_angles = range(angle_span[1], angle_span[2]; step=deg2rad(step_deg))
    problem = ODEProblem(advanced_rhs!, u0, angle_span, cfg)
    common = (
        reltol=reltol,
        abstol=[fill(1e-11, NSPECIES); 1e-7; fill(1e-8, 7); 1e-11],
        saveat=save_angles,
        tstops=filter(t -> angle_span[1] < t < angle_span[2],
            advanced_integration_stops(cfg)),
        isoutofdomain=(u, p, t) -> any(view(u, 1:NSPECIES) .< -1e-10),
    )

    if cfg.physics.combustion_model == :isochoric
        heat_input = nominal_external_heat(cfg)
        affect! = integrator -> begin
            integrator.u[I_U] += heat_input
            integrator.u[I_QRELEASE] += heat_input
        end
        callback = PresetTimeCallback([0.0], affect!; save_positions=(true, true))
        return solve(problem, Tsit5(); common..., callback=callback)
    elseif cfg.physics.combustion_model == :tcc_air_fuel
        stoichiometric_o2 = -reaction_stoichiometry(cfg)[I_O2]
        capture_fuel! = integrator -> begin
            available_fuel = max(integrator.u[I_FUEL], 0.0)
            oxygen_limited_fuel = max(integrator.u[I_O2], 0.0) / stoichiometric_o2
            integrator.u[I_BURNABLE_FUEL] = min(available_fuel, oxygen_limited_fuel)
        end
        ignition_angle = deg2rad(cfg.combustion.start_deg)
        callback = PresetTimeCallback([ignition_angle], capture_fuel!;
            save_positions=(true, true))
        return solve(problem, Tsit5(); common..., callback=callback)
    end
    return solve(problem, Tsit5(); common...)
end

function advanced_dataframe(solution, cfg::AdvancedEngineConfig)
    rows = length(solution.t)
    angle_deg = rad2deg.(solution.t)
    volume_cm3 = cylinder_volume.(Ref(cfg.geometry), solution.t) .* 1e6
    pressure_kPa = Vector{Float64}(undef, rows)
    temperature_K = Vector{Float64}(undef, rows)
    mass_g = Vector{Float64}(undef, rows)
    cp_molar = Vector{Float64}(undef, rows)
    cv_molar = Vector{Float64}(undef, rows)
    gamma = Vector{Float64}(undef, rows)

    species_umol = [Vector{Float64}(undef, rows) for _ in 1:NSPECIES]
    for i in eachindex(solution.t)
        u = solution.u[i]
        n = view(u, 1:NSPECIES)
        temperature = mixture_temperature(cfg, n, u[I_U])
        volume = cylinder_volume(cfg.geometry, solution.t[i])
        pressure = sum(max.(n, 0.0)) * RU * temperature / volume
        cp_i, cv_i, gamma_i = mixture_heat_capacities(cfg, n, temperature)
        pressure_kPa[i] = pressure / 1e3
        temperature_K[i] = temperature
        mass_g[i] = mixture_mass(cfg, n) * 1e3
        cp_molar[i] = cp_i
        cv_molar[i] = cv_i
        gamma[i] = gamma_i
        for j in 1:NSPECIES
            species_umol[j][i] = u[j] * 1e6
        end
    end

    return DataFrame(
        angle_deg=angle_deg,
        volume_cm3=volume_cm3,
        pressure_kPa=pressure_kPa,
        temperature_K=temperature_K,
        mass_g=mass_g,
        cp_molar_J_molK=cp_molar,
        cv_molar_J_molK=cv_molar,
        gamma=gamma,
        burned_fraction=wiebe_fraction.(Ref(cfg.combustion), solution.t),
        fuel_umol=species_umol[I_FUEL],
        O2_umol=species_umol[I_O2],
        N2_umol=species_umol[I_N2],
        CO2_umol=species_umol[I_CO2],
        H2O_umol=species_umol[I_H2O],
        indicated_work_J=[u[I_WIND] for u in solution.u],
        friction_work_J=[u[I_WFRIC] for u in solution.u],
        brake_work_J=[u[I_WIND] - u[I_WFRIC] for u in solution.u],
        released_energy_J=[u[I_QRELEASE] for u in solution.u],
        wall_heat_J=[u[I_QWALL] for u in solution.u],
        enthalpy_flow_J=[u[I_HFLOW] for u in solution.u],
        intake_mass_g=[u[I_MIN] * 1e3 for u in solution.u],
        exhaust_mass_g=[u[I_MOUT] * 1e3 for u in solution.u],
    )
end

function advanced_summary(data, cfg, initial, final, cycles)
    indicated_work = data.indicated_work_J[end]
    friction_work = data.friction_work_J[end]
    brake_work = data.brake_work_J[end]
    released_energy = data.released_energy_J[end]
    wall_heat = data.wall_heat_J[end]
    enthalpy_flow = data.enthalpy_flow_J[end]
    delta_u = final.internal_energy - initial.internal_energy
    external_heat = cfg.physics.combustion_model == :tcc_air_fuel ? 0.0 : released_energy
    energy_residual = delta_u - (external_heat - wall_heat - indicated_work + enthalpy_flow)
    unit_displacement = displacement(cfg.geometry)
    cycle_frequency = cfg.rpm / 120

    xint = intake_mole_fractions(cfg)
    _, _, intake_molar_mass = mixture_flow_properties(cfg, xint, cfg.intake.temperature)
    reference_mass = cfg.intake.pressure * intake_molar_mass * unit_displacement /
        (RU * cfg.intake.temperature)

    return DataFrame(
        metric=[
            "cycles_to_convergence",
            "indicated_work_J_per_cylinder",
            "friction_work_J_per_cylinder",
            "brake_work_J_per_cylinder",
            "released_energy_J_per_cylinder",
            "indicated_efficiency",
            "brake_efficiency",
            "imep_kPa",
            "bmep_kPa",
            "total_brake_power_kW",
            "volumetric_efficiency",
            "peak_pressure_kPa",
            "peak_temperature_K",
            "wall_heat_J_per_cylinder",
            "energy_residual_J",
        ],
        value=[
            Float64(cycles),
            indicated_work,
            friction_work,
            brake_work,
            released_energy,
            indicated_work / released_energy,
            brake_work / released_energy,
            indicated_work / unit_displacement / 1e3,
            brake_work / unit_displacement / 1e3,
            brake_work * cycle_frequency * cfg.geometry.cylinders / 1e3,
            data.intake_mass_g[end] / (reference_mass * 1e3),
            maximum(data.pressure_kPa),
            maximum(data.temperature_K),
            wall_heat,
            energy_residual,
        ],
        unit=["-", "J", "J", "J", "J", "-", "-", "kPa", "kPa", "kW",
            "-", "kPa", "K", "J", "J"],
    )
end

function simulate_advanced_cycle(cfg::AdvancedEngineConfig=default_advanced_config();
    initial_state=nothing, step_deg=0.25, reltol=1e-8)
    initial = isnothing(initial_state) ? advanced_initial_state(cfg) : initial_state
    solution = solve_advanced_cycle(cfg, initial; step_deg=step_deg, reltol=reltol)
    final = AdvancedState(solution.u[end][1:NSPECIES], solution.u[end][I_U])
    data = advanced_dataframe(solution, cfg)
    convergence = DataFrame(cycle=[1], mass_g=[mixture_mass(cfg, final.moles) * 1e3],
        temperature_K=[mixture_temperature(cfg, final.moles, final.internal_energy)],
        relative_change=[NaN])
    summary = advanced_summary(data, cfg, initial, final, 1)
    return AdvancedCycleResult(data, summary, initial, final, 1, convergence)
end

function simulate_advanced_to_periodic(cfg::AdvancedEngineConfig=default_advanced_config();
    step_deg=0.25, reltol=1e-8, periodic_tolerance=1e-6,
    min_cycles=3, max_cycles=100)
    cfg.physics.gas_exchange_model == :closed && throw(ArgumentError(
        "o trecho fechado termina antes da rejeição de calor; use simulate_advanced_cycle"))
    state = advanced_initial_state(cfg)
    history = DataFrame(cycle=Int[], mass_g=Float64[], temperature_K=Float64[],
        relative_change=Float64[])
    final_solution = nothing
    cycle_initial = state

    for cycle in 1:max_cycles
        cycle_initial = state
        solution = solve_advanced_cycle(cfg, cycle_initial; step_deg=step_deg, reltol=reltol)
        state = AdvancedState(solution.u[end][1:NSPECIES], solution.u[end][I_U])
        old_temperature = mixture_temperature(cfg, cycle_initial.moles,
            cycle_initial.internal_energy)
        new_temperature = mixture_temperature(cfg, state.moles, state.internal_energy)
        composition_change = maximum(abs.(state.moles .- cycle_initial.moles)) /
            max(sum(state.moles), 1e-12)
        mass_change = abs(mixture_mass(cfg, state.moles) -
            mixture_mass(cfg, cycle_initial.moles)) / max(mixture_mass(cfg, state.moles), 1e-12)
        temperature_change = abs(new_temperature - old_temperature) / max(new_temperature, 1e-12)
        relative_change = max(composition_change, mass_change, temperature_change)
        push!(history, (cycle, mixture_mass(cfg, state.moles) * 1e3,
            new_temperature, relative_change))
        final_solution = solution

        if cycle >= min_cycles && relative_change <= periodic_tolerance
            data = advanced_dataframe(solution, cfg)
            summary = advanced_summary(data, cfg, cycle_initial, state, cycle)
            return AdvancedCycleResult(data, summary, cycle_initial, state,
                cycle, history)
        end
    end

    @warn "O modelo avançado não atingiu periodicidade" max_cycles periodic_tolerance
    data = advanced_dataframe(final_solution, cfg)
    summary = advanced_summary(data, cfg, cycle_initial, state, max_cycles)
    return AdvancedCycleResult(data, summary, cycle_initial, state,
        max_cycles, history)
end

advanced_cycle_summary(result::AdvancedCycleResult) = result.summary

function plot_advanced_cycle(result::AdvancedCycleResult; output_path=nothing)
    d = result.data
    p1 = plot(d.angle_deg, d.pressure_kPa; xlabel="Ângulo [graus]",
        ylabel="Pressão [kPa]", label=false, linewidth=2, title="Pressão")
    p2 = plot(d.angle_deg, d.temperature_K; xlabel="Ângulo [graus]",
        ylabel="Temperatura [K]", label=false, linewidth=2, title="Temperatura")
    p3 = plot(d.volume_cm3, d.pressure_kPa; xlabel="Volume [cm³]",
        ylabel="Pressão [kPa]", label=false, linewidth=2, title="Diagrama p-V")
    p4 = plot(d.angle_deg, d.indicated_work_J; xlabel="Ângulo [graus]",
        ylabel="Trabalho acumulado [J]", label="Indicado", linewidth=2,
        title="Trabalho e atrito")
    plot!(p4, d.angle_deg, d.brake_work_J; label="Freio", linewidth=2)
    panel = plot(p1, p2, p3, p4; layout=(2, 2), size=(1100, 760))
    !isnothing(output_path) && savefig(panel, output_path)
    return panel
end

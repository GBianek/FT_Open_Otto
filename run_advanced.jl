using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using OpenOtto
using PrettyTables

# Caso solicitado: combustão ar-combustível do TCC sem cinética, cp(T),
# transferência de calor, válvulas e atrito linear F = μ v_p.
physics = PhysicsOptions(
    combustion_model=:tcc_air_fuel,
    friction_model=:linear_velocity,
    heat_transfer_model=:constant_h,
    gas_exchange_model=:compressible_valves,
    property_model=:nasa,
)

cfg = default_advanced_config(
    physics=physics,
    equivalence_ratio=1.0,
    friction_coefficient=16.0,
)

result = simulate_advanced_to_periodic(cfg; step_deg=0.25)

mkpath(joinpath(@__DIR__, "..", "output"))
figure_path = joinpath(@__DIR__, "..", "output", "ciclo_avancado.png")
plot_advanced_cycle(result; output_path=figure_path)

pretty_table(result.summary)
println("\nGráfico salvo em: $figure_path")

# Exemplo de estudo isolado: desligar apenas o atrito, mantendo todo o resto.
cfg_without_friction = with_physics(cfg, PhysicsOptions(
    combustion_model=:tcc_air_fuel,
    friction_model=:none,
    heat_transfer_model=:constant_h,
    gas_exchange_model=:compressible_valves,
    property_model=:nasa,
))

result_without_friction = simulate_advanced_cycle(
    cfg_without_friction;
    initial_state=result.initial_state,
    step_deg=0.25,
)

println("\nMesmo ciclo, sem atrito:")
pretty_table(result_without_friction.summary)

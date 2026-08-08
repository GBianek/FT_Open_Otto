using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using OpenOtto
using PrettyTables

config = default_config()

result = simulate_to_periodic(
    config;
    initial_pressure=config.exhaust.pressure,
    initial_temperature=750.0,
    step_deg=0.25,
    periodic_tolerance=1e-6,
    max_cycles=100,
)

pretty_table(result.summary)
pretty_table(last(result.convergence, min(10, size(result.convergence, 1))))

mkpath(joinpath(@__DIR__, "..", "output"))
plot_cycle(result; output_path=joinpath(@__DIR__, "..", "output", "ciclo_aberto.png"))

println("\nGráfico salvo em output/ciclo_aberto.png")

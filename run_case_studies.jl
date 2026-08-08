using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using OpenOtto
using PrettyTables

base = default_advanced_config(
    physics=PhysicsOptions(
        combustion_model=:tcc_air_fuel,
        friction_model=:linear_velocity,
        heat_transfer_model=:constant_h,
        gas_exchange_model=:compressible_valves,
        property_model=:nasa,
    ),
    rpm=2000.0,
    equivalence_ratio=1.0,
    friction_coefficient=16.0,
)

studies = run_default_case_studies(
    base;
    step_deg=0.5,
    periodic_tolerance=2e-6,
)

output_dir = joinpath(@__DIR__, "..", "output", "estudos_parametricos")
mkpath(output_dir)

for (slug, study) in pairs(studies)
    println("\n", study.name)
    pretty_table(parametric_summary(study))
    plot_case_study(
        study;
        output_path=joinpath(output_dir, "$(slug)_pv_palpha_tv.png"),
    )
end

println("\nPainéis salvos em: ", output_dir)

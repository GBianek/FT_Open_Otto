using Test
using OpenOtto

@testset "Geometria biela-manivela" begin
    cfg = default_config()
    g = cfg.geometry
    @test cylinder_volume(g, 0.0) ≈ cylinder_volume(g, 2π)
    @test cylinder_volume(g, π) / cylinder_volume(g, 0.0) ≈ g.compression_ratio
    @test dvolume_dangle(g, 0.0) ≈ 0.0 atol=1e-14
    @test dvolume_dangle(g, π) ≈ 0.0 atol=1e-14
end

@testset "Termoquímica e composição avançadas" begin
    cfg = default_advanced_config()
    nu = reaction_stoichiometry(cfg)
    reaction_mass = sum(nu[i] * species_molar_mass(cfg, i) for i in eachindex(nu))

    @test reaction_mass ≈ 0.0 atol=1e-12
    @test reaction_heat_per_mol(cfg) / 1e3 ≈ 802.56 rtol=2e-4
    @test species_cp(cfg, 4, 2000.0) > species_cp(cfg, 4, 300.0)
    @test sum(intake_mole_fractions(cfg)) ≈ 1.0
    @test sum(exhaust_mole_fractions(cfg)) ≈ 1.0
end

@testset "Ciclo avançado, atrito e primeira lei" begin
    cfg = default_advanced_config()
    result = simulate_advanced_cycle(cfg; step_deg=1.0, reltol=1e-7)
    summary = Dict(result.summary.metric .=> result.summary.value)

    @test all(result.data.mass_g .> 0)
    @test all(result.data.temperature_K .> 0)
    @test minimum(result.data.fuel_umol) >= -1e-3
    @test summary["released_energy_J_per_cylinder"] > 0
    @test summary["friction_work_J_per_cylinder"] > 0
    @test summary["brake_work_J_per_cylinder"] ≈
        summary["indicated_work_J_per_cylinder"] - summary["friction_work_J_per_cylinder"]
    @test abs(summary["energy_residual_J"]) <=
        1e-5 * summary["released_energy_J_per_cylinder"]
end

@testset "Configurações dos estudos paramétricos" begin
    base = default_advanced_config()
    changed = reconfigure_advanced(base;
        rpm=6000.0,
        equivalence_ratio=0.70,
        compression_ratio=16.0,
        property_model=:constant,
    )
    @test changed.rpm == 6000.0
    @test changed.equivalence_ratio == 0.70
    @test changed.geometry.compression_ratio == 16.0
    @test changed.physics.property_model == :constant
    @test changed.geometry.bore == base.geometry.bore
    @test changed.intake === base.intake

    definitions = default_case_studies(base)
    @test length(definitions.rotation.cases) == 4
    @test length(definitions.mixture.cases) == 3
    @test length(definitions.compression.cases) == 4
    @test length(definitions.properties.cases) == 2
end

@testset "Referência Otto ideal pelos switches" begin
    cfg = default_advanced_config(physics=ideal_physics())
    result = simulate_advanced_cycle(cfg; step_deg=0.5, reltol=1e-8)
    summary = Dict(result.summary.metric .=> result.summary.value)
    gamma = result.data.gamma[1]
    eta_otto = 1 - cfg.geometry.compression_ratio^(-(gamma - 1))

    @test summary["indicated_efficiency"] ≈ eta_otto rtol=2e-5
    @test summary["friction_work_J_per_cylinder"] == 0.0
    @test summary["wall_heat_J_per_cylinder"] == 0.0
end

@testset "Lei de Wiebe" begin
    w = WiebeLaw(-15.0, 50.0, 5.0, 2.0)
    @test wiebe_fraction(w, deg2rad(-20.0)) == 0.0
    @test wiebe_fraction(w, deg2rad(40.0)) == 1.0
    angles = range(deg2rad(-15.0), deg2rad(35.0); length=20_001)
    integral = sum((wiebe_derivative(w, angles[i]) + wiebe_derivative(w, angles[i + 1])) / 2 *
                   (angles[i + 1] - angles[i]) for i in 1:length(angles)-1)
    @test integral ≈ 1.0 rtol=2e-5
end

@testset "Vazão compressível" begin
    gas = GasProperties()
    mdot = orifice_mass_flow(1e-4, 0.7, 200e3, 300.0, 50e3, gas)
    @test isfinite(mdot)
    @test mdot > 0
    @test orifice_mass_flow(1e-4, 0.7, 100e3, 300.0, 100e3, gas) == 0.0
end

@testset "Balanço de um ciclo" begin
    result = simulate_cycle(default_config(); step_deg=1.0, reltol=1e-7)
    @test all(result.data.mass_g .> 0)
    @test all(result.data.temperature_K .> 0)
    residual = result.summary[result.summary.metric .== "energy_residual_J", :value][1]
    qin = result.summary[result.summary.metric .== "gross_heat_release_J_per_cylinder", :value][1]
    @test abs(residual) <= 1e-6 * max(abs(qin), 1.0)
end

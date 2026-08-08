export StudyCase,
       ParametricStudyResult,
       default_case_studies,
       run_case_study,
       run_default_case_studies,
       parametric_summary,
       plot_case_study

"""Um caso nomeado de uma bateria paramétrica."""
struct StudyCase
    label::String
    config::AdvancedEngineConfig
end

"""Resultados de uma bateria na mesma ordem dos casos de entrada."""
struct ParametricStudyResult
    name::String
    parameter::String
    cases::Vector{StudyCase}
    results::Vector{AdvancedCycleResult}

    function ParametricStudyResult(name, parameter, cases, results)
        length(cases) == length(results) ||
            throw(ArgumentError("cada caso deve possuir um resultado"))
        new(String(name), String(parameter), collect(cases), collect(results))
    end
end

"""
Casos padronizados para rotação, riqueza, razão de compressão e propriedades.

Cada bateria mantém os demais parâmetros no valor do `base`, evitando misturar
efeitos físicos distintos em uma mesma comparação.
"""
function default_case_studies(base::AdvancedEngineConfig=default_advanced_config())
    rotations = [StudyCase(label, reconfigure_advanced(base; rpm=value))
        for (value, label) in (
            (2000.0, "2000 rpm"),
            (4000.0, "4000 rpm"),
            (6000.0, "6000 rpm"),
            (8000.0, "8000 rpm"),
        )]
    mixtures = [StudyCase(label,
        reconfigure_advanced(base; equivalence_ratio=value))
        for (value, label) in (
            (0.70, "ϕ = 0,70"),
            (0.85, "ϕ = 0,85"),
            (1.00, "ϕ = 1,00"),
        )]
    compressions = [StudyCase(label,
        reconfigure_advanced(base; compression_ratio=value))
        for (value, label) in (
            (8.0, "r = 8"),
            (10.0, "r = 10"),
            (13.0, "r = 13"),
            (16.0, "r = 16"),
        )]
    properties = [
        StudyCase("cp constante", reconfigure_advanced(base; property_model=:constant)),
        StudyCase("cp(T), NASA", reconfigure_advanced(base; property_model=:nasa)),
    ]
    return (
        rotation=(name="Rotação", parameter="rpm", cases=rotations),
        mixture=(name="Composição da mistura", parameter="equivalence_ratio", cases=mixtures),
        compression=(name="Razão de compressão", parameter="compression_ratio", cases=compressions),
        properties=(name="Modelo de propriedades", parameter="property_model", cases=properties),
    )
end

function run_case_study(name, parameter, cases;
    step_deg=0.5, reltol=3e-8, periodic_tolerance=2e-6,
    min_cycles=3, max_cycles=60)
    results = AdvancedCycleResult[]
    for case in cases
        push!(results, simulate_advanced_to_periodic(
            case.config;
            step_deg=step_deg,
            reltol=reltol,
            periodic_tolerance=periodic_tolerance,
            min_cycles=min_cycles,
            max_cycles=max_cycles,
        ))
    end
    return ParametricStudyResult(name, parameter, cases, results)
end

function run_default_case_studies(base::AdvancedEngineConfig=default_advanced_config(); kwargs...)
    definitions = default_case_studies(base)
    return (
        rotation=run_case_study(definitions.rotation.name,
            definitions.rotation.parameter, definitions.rotation.cases; kwargs...),
        mixture=run_case_study(definitions.mixture.name,
            definitions.mixture.parameter, definitions.mixture.cases; kwargs...),
        compression=run_case_study(definitions.compression.name,
            definitions.compression.parameter, definitions.compression.cases; kwargs...),
        properties=run_case_study(definitions.properties.name,
            definitions.properties.parameter, definitions.properties.cases; kwargs...),
    )
end

function study_summary_value(result::AdvancedCycleResult, metric)
    row = findfirst(==(metric), result.summary.metric)
    isnothing(row) && throw(KeyError(metric))
    return result.summary.value[row]
end

"""Tabela compacta das grandezas mais úteis para interpretar cada bateria."""
function parametric_summary(study::ParametricStudyResult)
    table = DataFrame(
        study=String[],
        case=String[],
        cycles=Int[],
        pmax_MPa=Float64[],
        Tmax_K=Float64[],
        imep_kPa=Float64[],
        bmep_kPa=Float64[],
        eta_i=Float64[],
        eta_b=Float64[],
        volumetric_efficiency=Float64[],
        energy_residual_J=Float64[],
    )
    for (case, result) in zip(study.cases, study.results)
        push!(table, (
            study.name,
            case.label,
            result.cycles,
            study_summary_value(result, "peak_pressure_kPa") / 1e3,
            study_summary_value(result, "peak_temperature_K"),
            study_summary_value(result, "imep_kPa"),
            study_summary_value(result, "bmep_kPa"),
            study_summary_value(result, "indicated_efficiency"),
            study_summary_value(result, "brake_efficiency"),
            study_summary_value(result, "volumetric_efficiency"),
            study_summary_value(result, "energy_residual_J"),
        ))
    end
    return table
end

function add_study_curves!(chart, study, xcolumn, ycolumn;
    angle_window=nothing, yscale=1.0)
    for (case, result) in zip(study.cases, study.results)
        data = result.data
        mask = if isnothing(angle_window)
            trues(nrow(data))
        else
            (data.angle_deg .>= angle_window[1]) .&
                (data.angle_deg .<= angle_window[2])
        end
        x = data[mask, xcolumn]
        y = data[mask, ycolumn] .* yscale
        plot!(chart, x, y; label=case.label, linewidth=2)
    end
    return chart
end

"""Gera um painel horizontal com P×V, P×α e T×V para uma bateria."""
function plot_case_study(study::ParametricStudyResult; output_path=nothing)
    pv = plot(; xlabel="Volume V [cm³]", ylabel="Pressão p [MPa]",
        title="P × V", legend=:topright, gridalpha=0.25)
    pa = plot(; xlabel="Ângulo α [graus]", ylabel="Pressão p [MPa]",
        title="P × α", legend=:topright, gridalpha=0.25, xlims=(-180, 180))
    tv = plot(; xlabel="Volume V [cm³]", ylabel="Temperatura T [K]",
        title="T × V", legend=:topright, gridalpha=0.25)

    add_study_curves!(pv, study, :volume_cm3, :pressure_kPa; yscale=1e-3)
    add_study_curves!(pa, study, :angle_deg, :pressure_kPa;
        angle_window=(-180.0, 180.0), yscale=1e-3)
    add_study_curves!(tv, study, :volume_cm3, :temperature_K)

    panel = plot(pv, pa, tv; layout=(1, 3), size=(1650, 480),
        plot_title=study.name)
    !isnothing(output_path) && savefig(panel, output_path)
    return panel
end

# OpenOtto.jl

Modelo aberto, zero-dimensional e de tempo finito para um motor Otto de quatro
tempos. O projeto prossegue a formulação geométrica e angular do TCC de Gabriel
N. Bianek, mas substitui o sistema fechado com cinética química por balanços de
massa e energia em um volume de controle aberto.

## O que esta versão resolve

- geometria exata pistão-biela-manivela usada no TCC;
- ciclo completo de 720 graus;
- admissão e escape com válvulas de área variável;
- vazão compressível, inclusive bloqueamento e refluxo;
- primeira lei extensiva para sistema aberto;
- liberação de calor de tempo finito por Wiebe normalizada;
- troca de calor com as paredes;
- convergência para um ciclo periódico;
- DataFrame de resultados, tabela-resumo e quatro gráficos.

A versão 0.2 acrescentou um segundo solver, mantendo o solver 0.1 intacto:

- combustão ar-combustível do TCC com fração queimada prescrita, sem cinética;
- transporte de CH4, O2, N2, CO2, H2O, CO e H2;
- `cp(T)`, `cv(T)`, entalpia e energia interna por polinômios NASA;
- atrito viscoso do pistão `F = mu*v_p`, trabalho de freio e BMEP;
- *switches* independentes para combustão, atrito, parede, válvulas e propriedades.

A versão 0.3 acrescenta baterias paramétricas reproduzíveis para:

- rotações de 2000 a 8000 rpm;
- misturas metano-ar com `phi = 0.70`, `0.85` e `1.00`;
- razões de compressão `r = 8`, `10`, `13` e `16`;
- comparação isolada entre `cp` constante e `cp(T)` por NASA;
- painéis comparativos `P-V`, `P-alpha` e `T-V`.

Leia [MODELO_MATEMATICO.md](MODELO_MATEMATICO.md) para o núcleo 0.1,
[MODELO_AVANCADO.md](MODELO_AVANCADO.md) para as novas equações e
[VERIFICACAO.md](VERIFICACAO.md) para o benchmark numérico independente.

## Instalação e execução

No terminal, dentro da pasta do projeto:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/run_example.jl
```

O exemplo imprime a tabela-resumo e salva `output/ciclo_aberto.png`.

Para o modelo avançado solicitado:

```bash
julia --project=. examples/run_advanced.jl
```

Ele salva `output/ciclo_avancado.png`.

Para executar todos os estudos paramétricos:

```bash
julia --project=. examples/run_case_studies.jl
```

Os quatro painéis são salvos em `output/estudos_parametricos/` e as tabelas
resumidas são impressas no terminal.

Para usar no REPL ou em um notebook:

```julia
using Pkg
Pkg.activate("caminho/para/OpenOtto")
Pkg.instantiate()

using OpenOtto
using PrettyTables

cfg = default_config()
result = simulate_to_periodic(cfg)

pretty_table(result.summary)
first(result.data, 10)
plot_cycle(result)
```

## Modelo avançado e switches

```julia
physics = PhysicsOptions(
    combustion_model=:tcc_air_fuel,       # ou :isochoric, :wiebe_source
    friction_model=:linear_velocity,      # ou :none
    heat_transfer_model=:constant_h,      # ou :adiabatic
    gas_exchange_model=:compressible_valves, # ou :closed
    property_model=:nasa,                 # ou :constant
)

cfg = default_advanced_config(
    physics=physics,
    equivalence_ratio=1.0,
    friction_coefficient=16.0,
)

result = simulate_advanced_to_periodic(cfg)
plot_advanced_cycle(result)
```

Uma configuração existente pode ser alterada sem reconstruí-la manualmente:

```julia
cfg_6000 = reconfigure_advanced(cfg; rpm=6000.0)
cfg_lean = reconfigure_advanced(cfg; equivalence_ratio=0.70)
cfg_r16 = reconfigure_advanced(cfg; compression_ratio=16.0)
cfg_cp = reconfigure_advanced(cfg; property_model=:constant)
```

Todos os switches idealizados podem ser reunidos e comparados com a expressão
analítica do ciclo Otto:

```julia
ideal = simulate_advanced_cycle(default_advanced_config(physics=ideal_physics()))
```

Nesse modo fechado, a integração percorre o trecho compressão–combustão–expansão
de `-180°` a `+180°`; não use a rotina de convergência periódica.

O caso-base avançado usa metano e aceita `0 < equivalence_ratio <= 1`. A
composição rica e a dissociação ficam para uma etapa posterior.

## Parâmetros de referência

`default_config()` preserva o caso-base do TCC sempre que ele é aplicável:

| Parâmetro | Valor |
|---|---:|
| Cilindrada total | 1667 cm³ |
| Cilindros | 6 |
| Motor | quadrado (`D=S`) |
| Razão `L/R` | 3,7 |
| Razão de compressão | 13 |
| Rotação | 2000 rpm |
| Admissão | 101,35 kPa; 303,15 K |

Os eventos de válvulas, a contrapressão de escape, a parede e a Wiebe são
parâmetros iniciais plausíveis, não dados calibrados de um motor específico.

## Personalização

Os tipos são imutáveis para tornar cada caso explícito e reproduzível. Exemplo:

```julia
g = EngineGeometry(0.086, 0.086, 0.1591, 10.5, 4)
gas = GasProperties(R=287.0, gamma=1.35)

cfg = EngineConfig(
    g,
    gas,
    Reservoir(95e3, 300.0),
    Reservoir(115e3, 800.0),
    Valve(-375.0, -140.0, 3.2e-4, 0.70),
    Valve(135.0, 375.0, 2.8e-4, 0.72),
    WiebeLaw(-18.0, 48.0, 5.0, 2.0);
    rpm=2500.0,
    wall_temperature=450.0,
    heat_transfer_coefficient=400.0,
    equivalence_ratio=0.9,
)

result = simulate_to_periodic(cfg; step_deg=0.25)
```

## Convenção angular

O PMS de combustão é `0°`. O intervalo simulado é:

| Intervalo idealizado | Tempo predominante |
|---|---|
| `-360°` a `-180°` | admissão |
| `-180°` a `0°` | compressão |
| `0°` a `180°` | expansão |
| `180°` a `360°` | escape |

Os eventos reais das válvulas podem ultrapassar esses limites e gerar overlap.
A área é repetida com período de 720°, de modo que uma janela como `-380°` a
`-130°` continua corretamente pelo contorno `+360°/-360°`.

## Testes

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Os testes verificam a geometria, a integral da Wiebe, a vazão compressível e o
fechamento numérico da primeira lei.

## Limites físicos

O código ainda não é um modelo preditivo calibrado de um motor real. O solver
0.1 conserva `cp`, `cv` e composição constantes; o solver 0.2 remove essas
hipóteses, mas ainda não trata mistura rica, dissociação, NOx, blow-by, crevices
ou dinâmica dos coletores. A transferência de calor e o atrito são modelos
concentrados. É uma base conservativa e extensível, não um substituto de
validação experimental.

## Relação com o FTOtto

O repositório anterior discretizava processos de um sistema fechado e corrigia
um expoente politrópico para satisfazer a primeira lei. Aqui, `m` e `U` são
integrados diretamente. Isso elimina a necessidade do expoente politrópico e
permite que as válvulas alterem simultaneamente massa, entalpia e pressão.

Referências de base: Brunetti (2012), Naaktgeboren (2017) e o TCC de Bianek
(2024). A formulação do escoamento em orifício segue as relações isentrópicas
usuais para gases perfeitos.

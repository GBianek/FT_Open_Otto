# Verificação de referência

Como o ambiente usado para preparar esta versão não continha um executável
Julia, foram feitas duas verificações separadas:

1. análise sintática de todos os arquivos `.jl` com uma gramática Julia;
2. reprodução independente das mesmas equações em Python/SciPy, sem reutilizar
   o código Julia.

Para `default_config()`, integração máxima de 0,25 grau e tolerância relativa de
`1e-9`, a reprodução independente forneceu:

| Grandeza | Resultado de referência |
|---|---:|
| Ciclos até periodicidade | 4 |
| Mudança relativa final | `1,05e-8` |
| Calor liberado por cilindro | 933,485 J |
| Trabalho indicado líquido por cilindro | 480,444 J |
| Rendimento indicado | 0,51468 |
| IMEP | 1729,25 kPa |
| Pressão máxima | 9368,11 kPa |
| Temperatura máxima | 3747,99 K |
| Resíduo da primeira lei | `-4,70e-13 J` |

A temperatura máxima elevada é uma consequência esperada da aproximação de
calores específicos constantes e da ausência de dissociação. Ela não deve ser
interpretada como previsão calibrada. Ao trocar por propriedades dependentes de
temperatura e uma composição transportada, esse pico deverá diminuir.

Ao executar o projeto em Julia, pequenas diferenças são normais por causa do
controlador de passo do `Tsit5`. Diferenças grandes devem ser investigadas com
`julia --project=. -e 'using Pkg; Pkg.test()'`.

## Verificação independente da versão 0.2

O mesmo procedimento foi repetido para o solver avançado, incluindo os
polinômios NASA, os sete balanços de espécie, a inversão `U -> T`, o escoamento
de entalpia química e o atrito `F=mu*v_p`. O caso é o de
`default_advanced_config()`, com CH4, `phi=1`, `mu=16 kg/s`, passo máximo de
`0,5°` e tolerância relativa de `3e-8` na reprodução Python/SciPy.

| Grandeza | Resultado de referência v0.2 |
|---|---:|
| Ciclos até periodicidade | 5 |
| Mudança relativa final | `1,64e-7` |
| Energia química liberada por cilindro | 750,681 J |
| Trabalho indicado por cilindro | 329,814 J |
| Perda por atrito por cilindro | 26,828 J |
| Trabalho de freio por cilindro | 302,986 J |
| Rendimento indicado | 0,43935 |
| Rendimento de freio | 0,40361 |
| IMEP | 1187,09 kPa |
| BMEP | 1090,53 kPa |
| Potência de freio, seis cilindros | 30,299 kW |
| Pressão máxima | 6660,52 kPa |
| Temperatura máxima | 2647,06 K |
| Resíduo da primeira lei | `-4,90e-13 J` |

O balanço atômico da reação de metano apresentou resíduo de massa de
`1,39e-17 kg/mol`, no nível de arredondamento de ponto flutuante. A entalpia de
reação calculada a 298,15 K foi `802,557 kJ/mol`, correspondente à água no
estado gasoso.

## Retorno ao ciclo Otto ideal

Com `ideal_physics()`, o trecho fechado usa mistura fresca, propriedades
constantes, parede adiabática, calor isocórico e atrito nulo. Para a mistura
CH4/ar do caso-base, `gamma=1,387039` e `r=13`. A referência analítica é

```math
eta_{Otto}=1-r^{-(gamma-1)}=0,6294384.
```

O teste automatizado exige que o rendimento obtido por integração coincida com
essa expressão dentro de `2e-5` de erro relativo. Há ainda testes separados para
conservação de massa da reação, entalpia de combustão, tendência de `cp(T)`,
positividade dos estados, identidade `W_b=W_i-W_f` e fechamento da primeira lei.

Como não havia executável Julia no ambiente de empacotamento, os arquivos Julia
foram analisados com uma gramática da linguagem, e a validação numérica foi feita
por implementação independente. A suíte `Pkg.test()` permanece como verificação
final no ambiente Julia do usuário.

## Estudos paramétricos da versão 0.3

O espelho independente foi generalizado para aceitar rotação, riqueza da
mistura, razão de compressão e modelo de propriedades. Todos os casos atingiram
periodicidade, mantiveram estados positivos e fecharam a primeira lei entre
`1e-14` e `2e-12 J`.

| Estudo | Caso | pmax [MPa] | Tmax [K] | IMEP [kPa] | BMEP [kPa] | eta_i | eta_b | eta_v |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Rotação | 2000 rpm | 6,660 | 2647 | 1187 | 1091 | 0,439 | 0,404 | 0,897 |
| Rotação | 4000 rpm | 5,670 | 2662 | 919 | 725 | 0,405 | 0,320 | 0,740 |
| Rotação | 6000 rpm | 3,990 | 2597 | 555 | 265 | 0,357 | 0,171 | 0,506 |
| Rotação | 8000 rpm | 3,016 | 2522 | 363 | -24 | 0,321 | -0,021 | 0,368 |
| Mistura | phi=0,70 | 5,668 | 2185 | 880 | 783 | 0,454 | 0,404 | 0,894 |
| Mistura | phi=0,85 | 6,182 | 2425 | 1039 | 942 | 0,447 | 0,405 | 0,896 |
| Mistura | phi=1,00 | 6,660 | 2647 | 1187 | 1091 | 0,439 | 0,404 | 0,897 |
| Compressão | r=8 | 4,377 | 2642 | 1037 | 941 | 0,379 | 0,344 | 0,907 |
| Compressão | r=10 | 5,320 | 2648 | 1109 | 1013 | 0,408 | 0,372 | 0,902 |
| Compressão | r=13 | 6,660 | 2647 | 1187 | 1091 | 0,439 | 0,404 | 0,897 |
| Compressão | r=16 | 7,950 | 2642 | 1243 | 1147 | 0,462 | 0,426 | 0,894 |
| Propriedades | cp constante | 8,463 | 3313 | 1423 | 1326 | 0,528 | 0,492 | 0,895 |
| Propriedades | cp(T), NASA | 6,660 | 2647 | 1187 | 1091 | 0,439 | 0,404 | 0,897 |

O resultado a 8000 rpm é uma extrapolação deliberada. A combinação de válvulas
fixas, enchimento reduzido e atrito linear produz trabalho de freio negativo;
isso identifica o limite do conjunto de parâmetros, não uma rotação máxima
universal. A comparação de propriedades mostra por que `cp` constante não deve
ser usado quantitativamente em altas temperaturas: ele superestima tanto a
temperatura quanto a pressão e o trabalho do ciclo.

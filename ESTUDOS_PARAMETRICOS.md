# Estudos paramétricos

Esta bateria usa como referência o OpenOtto avançado a 2000 rpm, `phi=1`,
`r=13`, propriedades NASA, válvulas compressíveis, parede com coeficiente
constante e atrito `F=mu*v_p`. Em cada estudo somente o parâmetro indicado é
alterado.

## Casos

| Bateria | Valores |
|---|---|
| Rotação | 2000, 4000, 6000 e 8000 rpm |
| Composição metano-ar | `phi=0.70`, `0.85` e `1.00` |
| Razão de compressão | `r=8`, `10`, `13` e `16` |
| Propriedades | `cp` constante e `cp(T)` por NASA |

O comando de execução é:

```bash
julia --project=. examples/run_case_studies.jl
```

## Leitura física

### Rotação

Com áreas e fases de válvula fixas, o menor tempo disponível reduz o enchimento:
a eficiência volumétrica cai de 0,897 a 2000 rpm para 0,368 a 8000 rpm. Ao mesmo
tempo, o trabalho de atrito por ciclo cresce linearmente com a rotação no modelo
`F=mu*v_p`. Por isso o BMEP torna-se negativo no caso extremo de 8000 rpm. Este
caso deve ser interpretado como teste de estresse do modelo concentrado.

### Composição

A mistura mais pobre contém menos combustível por quantidade de carga. O pico
de temperatura cai de 2647 K em `phi=1` para 2185 K em `phi=0.70`, enquanto o
IMEP diminui de 1187 para 880 kPa. O rendimento indicado cresce ligeiramente
porque as perdas térmicas diminuem, mas o trabalho disponível é menor.

### Razão de compressão

Entre `r=8` e `r=16`, o pico de pressão cresce de 4,38 para 7,95 MPa e o
rendimento indicado cresce de 0,379 para 0,462. A temperatura máxima varia
pouco neste conjunto aberto porque a massa admitida, os gases residuais, a
troca térmica e a combustão prescrita se reajustam em cada ciclo periódico.

### Propriedades

O caso com `cp` constante prevê 3313 K e 8,46 MPa, contra 2647 K e 6,66 MPa com
`cp(T)`. Como o calor específico real aumenta com a temperatura, o modelo NASA
absorve mais energia sensível por kelvin e reduz os picos. A hipótese constante
também superestima o rendimento indicado: 0,528 contra 0,439.

## Limitações

- não há calibração experimental das válvulas, da parede ou do atrito;
- o perfil de Wiebe é mantido em graus para todas as rotações;
- não há mistura rica, dissociação, detonação ou limite de propagação de chama;
- os coletores permanecem como reservatórios de pressão e temperatura fixas;
- 8000 rpm está fora da região em que o atrito linear deve ser considerado
  quantitativamente confiável.

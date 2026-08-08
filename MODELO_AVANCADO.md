# Modelo avançado: combustão prescrita, composição, `cp(T)` e atrito

## 1. O que mudou na versão 0.2

O solver original da versão 0.1 foi mantido para reprodutibilidade. O solver
avançado integra os números de mols de sete espécies, a energia interna total e
estados auxiliares de trabalho e transferência. Nesta etapa, o combustível de
referência é o metano e a mistura deve ser pobre ou estequiométrica
(`0 < phi <= 1`). Não há equação de cinética química.

As opções físicas são independentes:

| Hipótese ou fenômeno | Opção idealizada | Opção removida/realista |
|---|---|---|
| Combustão | `:isochoric` | `:wiebe_source` ou `:tcc_air_fuel` |
| Atrito | `:none` | `:linear_velocity` |
| Parede adiabática | `:adiabatic` | `:constant_h` |
| Sistema fechado | `:closed` | `:compressible_valves` |
| Calores específicos constantes | `:constant` | `:nasa` |

Cada linha é um *switch* de `PhysicsOptions`; portanto, um estudo pode alterar
uma hipótese por vez sem mudar as demais.

Com `gas_exchange_model=:closed`, o solver percorre somente o trecho de alta
pressão, de `-180°` a `+180°`, iniciado no PMI com mistura fresca. Essa opção é
útil para a referência Otto ar-padrão. Como a rejeição de calor isocórica não
faz parte desse trecho, use `simulate_advanced_cycle`; a rotina periódica rejeita
deliberadamente essa combinação.

## 2. Estado e primeira lei

O vetor de composição é

```math
\mathbf n=(n_F,n_{O_2},n_{N_2},n_{CO_2},n_{H_2O},n_{CO},n_{H_2}).
```

Para gases ideais, a energia e a pressão são

```math
U(T,\mathbf n)=\sum_j n_j\,[\bar h_j(T)-R_uT],
\qquad
p=\frac{n_{tot}R_uT}{V}.
```

`U` inclui energia sensível e entalpias de formação. Em cada avaliação do lado
direito da EDO, o programa inverte numericamente `U(T,n)` para obter `T`. A
primeira lei angular do volume de controle é

```math
\frac{dU}{d\alpha}=
\frac{dQ_{ext}}{d\alpha}
-\frac{\dot Q_w}{\omega}
-p\frac{dV}{d\alpha}
+\frac{1}{\omega}\sum_k \dot N_k\bar h_k.
```

Na combustão ar-combustível, `Q_ext=0`: a elevação de temperatura decorre da
mudança de composição a energia total conservada. Isso evita somar uma fonte de
calor ao mesmo tempo em que a entalpia química dos reagentes desaparece.

## 3. Combustão do TCC sem cinética

Para um combustível genérico `C_c H_h O_o N_n`, em mistura pobre ou
estequiométrica, usa-se

```math
C_cH_hO_oN_n + \nu_{O_2}O_2 \rightarrow
cCO_2+\frac{h}{2}H_2O+\frac{n}{2}N_2,
\qquad
\nu_{O_2}=c+\frac{h}{4}-\frac{o}{2}.
```

O ar admitido possui `O2 + 3.76 N2`, e a razão de equivalência controla o
excesso de oxigênio. A quantidade que pode reagir é capturada no instante de
ignição e limitada pelo combustível e pelo oxigênio aprisionados. A evolução é
prescrita pela fração queimada normalizada do TCC,

```math
\frac{d\mathbf n}{d\alpha}=\boldsymbol\nu\,n_{F,0}
\frac{dy}{d\alpha}.
```

Logo, a duração angular ainda é um dado da lei de Wiebe, mas não há constante
de Arrhenius, energia de ativação ou mecanismo cinético. O calor químico
acumulado mostrado nos resultados é somente um diagnóstico:

```math
Q_{rel}=n_{F,queimado}\left[-\sum_j\nu_j\bar h_j^0\right].
```

Os gases residuais não precisam da correlação fechada `zeta(P,r_v)` do TCC:
eles permanecem no cilindro porque admissão e escape são resolvidos como um
sistema aberto e o estado é repetido até a periodicidade.

## 4. Propriedades dependentes da temperatura

Para `property_model=:nasa`, cada espécie usa o polinômio NASA de sete
coeficientes em dois intervalos:

```math
\frac{\bar c_p}{R_u}=a_1+a_2T+a_3T^2+a_4T^3+a_5T^4,
```

```math
\frac{\bar h}{R_uT}=a_1+\frac{a_2T}{2}+\frac{a_3T^2}{3}
+\frac{a_4T^3}{4}+\frac{a_5T^4}{5}+\frac{a_6}{T}.
```

As propriedades da mistura são médias molares instantâneas. Assim, `cp`, `cv`
e `gamma` mudam tanto com `T` quanto com a composição. Os coeficientes
embutidos são os do GRI-Mech 3.0 para CH4, O2, N2, CO2, H2O, CO e H2.

## 5. Atrito linear

O modelo solicitado é aplicado ao pistão:

```math
F_f=\mu v_p,\qquad
v_p=\omega\frac{dx}{d\alpha}.
```

Como o atrito sempre dissipa energia, o trabalho perdido por ângulo é

```math
\frac{dW_f}{d\alpha}=|F_f\,dx/d\alpha|
=\mu\omega\left(\frac{dx}{d\alpha}\right)^2.
```

Com rotação prescrita, esse termo não altera a trajetória termodinâmica do gás:
ele converte trabalho indicado em perda mecânica. Portanto,

```math
W_b=W_i-W_f,
\qquad BMEP=\frac{W_b}{V_d}.
```

O valor padrão `mu=16 kg/s` é apenas um caso-base e deve ser calibrado com
dados de motorização ou pressão no cilindro.

## 6. Limites desta etapa

- somente misturas pobres/estequiométricas; ainda não há CO e H2 de mistura rica;
- sem dissociação, NOx, blow-by, crevices ou dinâmica dos coletores;
- metano é o combustível fornecido de fábrica; outros combustíveis exigem seus
  coeficientes NASA;
- a duração de combustão é prescrita, não predita;
- `h` da parede e `mu` são parâmetros concentrados que requerem calibração.

Esses limites são deliberados: isolam primeiro os efeitos de combustão finita,
composição, `cp(T)` e atrito antes da volta da cinética química.

## 7. Fontes técnicas

- GRI-Mech 3.0, arquivos termoquímicos e definição dos polinômios NASA:
  <https://combustion.berkeley.edu/gri-mech/version30/text30.html>;
- P. L. Curto-Risso, A. Medina e A. Calvo Hernández, modelo Otto irreversível
  com atrito, transferência de calor e perdas internas:
  <https://doi.org/10.1063/1.2986214>;
- G. N. Bianek (2024), TCC fornecido pelo autor, especialmente as Equações 9 a
  14, 24 a 30 e 61 a 68.

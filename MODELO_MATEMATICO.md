# Modelo matemático do motor Otto aberto de tempo finito

## 1. Escopo da versão 0.1

O cilindro é tratado como um volume de controle zero-dimensional, homogêneo e
quase-estacionário. O modelo percorre 720 graus de virabrequim e resolve
admissão, compressão, liberação finita de calor, expansão e escape. A cinética
química do TCC foi retirada: a taxa de liberação de calor é prescrita por uma
função de Wiebe simples.

Nesta primeira abordagem:

- o fluido é um gás ideal com `R`, `cp` e `cv` constantes;
- os coletores são reservatórios de `p` e `T` constantes;
- cada válvula é um orifício compressível quase-estacionário;
- o interior do cilindro possui uma única pressão e uma única temperatura;
- a transferência de calor usa um coeficiente médio constante;
- não há atrito mecânico, blow-by, dinâmica dos coletores nem espécies químicas.

## 2. Geometria do conjunto pistão-biela-manivela

Com `R_c = S/2` como raio da manivela, `L` como comprimento da biela e `alpha`
medido a partir do PMS, a posição do pistão é

```math
x(\alpha)=R_c(1-\cos\alpha)+L-\sqrt{L^2-R_c^2\sin^2\alpha}.
```

Essa é a mesma geometria das Equações 5 a 8 do TCC. Para `A_p = pi D^2/4`,

```math
V(\alpha)=V_c+A_p x(\alpha),
\qquad
V_c=\frac{V_d}{r_v-1},
\qquad
V_d=A_pS.
```

A derivada usada diretamente no trabalho de fronteira é

```math
\frac{dV}{d\alpha}=A_p\left[
R_c\sin\alpha+
\frac{R_c^2\sin\alpha\cos\alpha}
{\sqrt{L^2-R_c^2\sin^2\alpha}}
\right].
```

## 3. Balanços no volume de controle

As variáveis conservativas são a massa `m` e a energia interna extensiva `U`.
Com velocidade angular constante `omega`, usa-se o ângulo como variável
independente (`dt = d alpha / omega`). O balanço de massa é

```math
\frac{dm}{d\alpha}=\frac{\dot m_i+\dot m_e}{\omega},
```

onde cada vazão é algébrica e positiva quando entra no cilindro. A primeira lei
para o volume de controle é

```math
\frac{dU}{d\alpha}=
\frac{dQ_{comb}}{d\alpha}
-\frac{\dot Q_w}{\omega}
-p\frac{dV}{d\alpha}
+\frac{\dot H_i+\dot H_e}{\omega}.
```

As propriedades instantâneas são recuperadas por

```math
T=\frac{U}{m c_v},
\qquad
p=\frac{mRT}{V}.
```

Para uma porta com vazão algébrica `mdot`, o transporte de entalpia é

```math
\dot H =
\begin{cases}
\dot m c_p T_{res}, & \dot m\ge 0,\\
\dot m c_p T,       & \dot m<0.
\end{cases}
```

Assim, refluxos também são tratados sem trocar manualmente a equação.

## 4. Vazão nas válvulas

A área instantânea é uma primeira aproximação suave:

```math
A_v(\alpha)=A_{max}\sin^2\left[
\pi\frac{\alpha-\alpha_{abre}}
{\alpha_{fecha}-\alpha_{abre}}
\right]
```

dentro da janela de abertura e zero fora dela. A vazão de montante para jusante
é calculada por escoamento isentrópico em orifício. Para `r_p=p_d/p_u`,
As janelas são periódicas em 720 graus, preservando o overlap no contorno entre
o fim e o início do ciclo numérico.

```math
\dot m=C_d A_v\frac{p_u}{\sqrt{RT_u}}\Phi(r_p),
```

com

```math
\Phi=\sqrt{\gamma\left(\frac{2}{\gamma+1}\right)^{
(\gamma+1)/(\gamma-1)}}
```

no regime bloqueado e

```math
\Phi=\sqrt{\frac{2\gamma}{\gamma-1}
\left(r_p^{2/\gamma}-r_p^{(\gamma+1)/\gamma}\right)}
```

no regime não bloqueado. O código identifica o lado de montante pela pressão e
atribui o sinal correto à vazão no cilindro.

## 5. Liberação finita de calor

Se `z=(alpha-alpha_0)/Delta alpha`, a fração queimada é

```math
x_b=\frac{1-\exp[-a z^{m_w+1}]}{1-\exp(-a)},
\qquad 0<z<1.
```

A normalização, ausente na forma convencional assintótica, garante `x_b=1` ao
fim da duração prescrita. O termo da primeira lei é

```math
\frac{dQ_{comb}}{d\alpha}=Q_{ciclo}\frac{dx_b}{d\alpha}.
```

Sem espécies, `Q_ciclo` é estimado a partir da carga nominal no PMI nas
condições do coletor de admissão:

```math
m_{carga,ref}=\frac{p_i V_{PMI}}{RT_i},
\qquad
m_{comb,ref}=m_{carga,ref}\frac{\phi}{AFR_{st}+\phi},
\qquad
Q_{ciclo}=\eta_c PCI\,m_{comb,ref}.
```

Essa escolha é deliberadamente simples. A próxima extensão natural é transportar
massas de ar, combustível e gases residuais separadamente e calcular o calor com
a massa de combustível efetivamente aprisionada no fechamento da admissão.

## 6. Transferência de calor e métricas

A área molhada aproximada é

```math
A_w(\alpha)=2A_p+\pi D x(\alpha),
```

e

```math
\dot Q_w=h_w A_w(T-T_w).
```

O trabalho indicado é integrado como estado auxiliar,

```math
W_i=\oint p\,dV,
```

e então

```math
IMEP=\frac{W_i}{V_d},
\qquad
P_i=W_i\frac{N}{120}z.
```

O fator `N/120` é a frequência de ciclos de um motor de quatro tempos.

## 7. Fechamento periódico

Um único ciclo depende do estado residual arbitrado no PMS de admissão. Por
isso, a rotina principal aplica repetidamente o mapa de 720 graus até que massa
e temperatura no fim do ciclo coincidam com as do início dentro da tolerância.
O último ciclo é então usado para os diagramas e métricas.

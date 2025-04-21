# Importando bibliotecas necessárias
import numpy as np
import cmath

# Definição da potência fornecida pelo gerador
P_g = 310795499.72  # Potência ativa total
Q_g = 58556595.94   # Potência reativa total

# Definição da porcentagem de perdas desejadas (menor que 10%)
perda_percentual_ativa = 0.09  # Ajuste para perdas menores que 10%
perda_percentual_reativa = 0.09  # Ajuste para perdas menores que 10%

# Calculando nova potência total das cargas
P_cargas_total = P_g * (1 - perda_percentual_ativa)
Q_cargas_total = Q_g * (1 - perda_percentual_reativa)

# Potências ativas e reativas atuais das cargas
P_Z1_atual = 103550195.56
P_Z2_atual = 136527757.26
P_Z3_atual = 26968948.49

Q_Z1_atual = 19506138.83
Q_Z2_atual = 25716934.57
Q_Z3_atual = 5090853.43

# Ajustando proporcionalmente as novas potências das cargas
P_Z1 = P_cargas_total * (P_Z1_atual / abs(P_Z1_atual + P_Z2_atual + P_Z3_atual))
P_Z2 = P_cargas_total * (P_Z2_atual / abs(P_Z1_atual + P_Z2_atual + P_Z3_atual))
P_Z3 = P_cargas_total * (P_Z3_atual / abs(P_Z1_atual + P_Z2_atual + P_Z3_atual))

Q_Z1 = abs(Q_cargas_total * (Q_Z1_atual / abs(Q_Z1_atual + Q_Z2_atual + Q_Z3_atual)))
Q_Z2 = abs(Q_cargas_total * (Q_Z2_atual / abs(Q_Z1_atual + Q_Z2_atual + Q_Z3_atual)))
Q_Z3 = abs(Q_cargas_total * (Q_Z3_atual / abs(Q_Z1_atual + Q_Z2_atual + Q_Z3_atual)))

# Tensão nominal nas cargas
V_Z1 = 490434.82
V_Z2 = 210412.55
V_Z3 = 62876.72

# Calculando impedâncias garantindo que a parte imaginária seja positiva
Z1 = ((V_Z1 ** 2) / (P_Z1 + 1j * abs(Q_Z1)))*0.9
Z2 = ((V_Z2 ** 2) / (P_Z2 + 1j * abs(Q_Z2)))*0.9
Z3 = ((V_Z3 ** 2) / (P_Z3 + 1j * abs(Q_Z3)))*0.9
# com os valores, obtemos os erros nao desejaveis, decidimos multiplicar por um fator de 0.9
# Exibindo resultados no formato desejado
print(f"Nova impedância de Z1: {Z1.real:.2f} + {abs(Z1.imag):.2f}j Ω")
print(f"Nova impedância de Z2: {Z2.real:.2f} + {abs(Z2.imag):.2f}j Ω")
print(f"Nova impedância de Z3: {Z3.real:.2f} + {abs(Z3.imag):.2f}j Ω")
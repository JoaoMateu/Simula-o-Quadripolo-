# Importando bibliotecas necessárias
import numpy as np
import cmath

#------------------Declarando as funções---------------------------------#

# MATRIZ PI(LINHA DE TRANSMISSÃO):
def M_Linha(Linha):
 # Admitância Y da linha (Y1=Y2)
 Y = 1/(1/(1j*w*(Linha[2]/2)))
 # Impedância Z da linha
 Z = Linha[0] + (1j*w*Linha[1])
 # Calculando matriz ABCD
 A = 1 + (Y*Z)
 B = Z
 C = (2*Y) + (Y*Y*Z)
 D = 1 + (Y*Z)
 matrizLinha = np.array([[A, B], [C, D]])
 return matrizLinha

# Função Impedância em série
def Carga(Z):
 A = 1
 B = Z
 C = 0
 D = 1
 matriz = np.array([[A, B],[C, D]])
 return matriz

# Função Matriz Admitância
def CargaDerivada(Z):
 A = 1
 B = 0
 C = 1/Z
 D = 1
 matriz = np.array([[A, B],[C, D]])
 return matriz

# Função MATRIZ T (TRANSFORMADOR):
def Transformador(Transfp):
 Z1 = 7.6e-3 + 1j*3.8e-3 # Indicando as impedâncias Z1 e Z2 do circuito equivalente do transformador
 Z2 = 33.9e-3 + 1j*0.85e-3
 V1 = Transfp[0] # Recebendo os valores de entrada
 V2 = Transfp[1]
 Rm = Transfp[2]
 Xm = Transfp[3]
 Y = (Rm + 1j*Xm) / (Rm*1j*Xm) # Indicando a admitância Y do circuito equivalente do transformador

 N = V1 / V2 # Relação V1/V2
 A = N * (1 + (Y*Z1)) # Construção da Matriz de transmissão
 B = (Z1 + Z2 + Y*Z1*Z2) / N
 C = N * Y
 D = (1 + Y*Z2) / N
 matriz = np.array([[A, B],[C, D],])
 return matriz

# FUNÇÃO MATRIZ CASCATA:
def Associar_Matriz_em_Cascata(matriz1, matriz2):
 bloco = np.dot(matriz1, matriz2)
 return bloco

# FUNÇÃO MATRIZ PARALELO:
def Matriz_Paralelo(matriz1, matriz2):
 Aa = matriz1[0][0]
 Ba = matriz1[0][1]
 Ca = matriz1[1][0]
 Da = matriz1[1][1]
 Ab = matriz2[0][0]
 Bb = matriz2[0][1]
 Cb = matriz2[1][0]
 Db = matriz2[1][1]
 A = (Aa*Bb + Ab*Ba)/(Ba + Bb)
 B = (Ba*Bb)/(Ba + Bb)
 C = Ca + Cb + ((Aa - Ab)*(Db - Da))/(Ba + Bb)
 D = (Bb*Da + Ba*Db)/(Ba + Bb)
 matriz = np.array([[A, B],[C, D],])
 return matriz

#------------------Declarando as Variáveis---------------------------------#

w = 2 * np.pi * 60 #--> Frequência do sistema (rad/s)

# Impedância série Thevenin
Zf = 2 + 1j*0.38

#distâncias das linhas de transmissão (em km)
DistLT1 = DistLT2 = 100
DistLT3 = 80
DistLT4 = 120

# parâmetros da LINHAS:
LT1_LT2 = [DistLT1*0.182, DistLT1*2.28e-3, DistLT1*0.0140e-6]
LT3_ = [DistLT3*0.182, DistLT3*2.28e-3, DistLT3*0.0140e-6]
LT4_ = [DistLT4*0.182, DistLT4*2.28e-3, DistLT4*0.0140e-6]

#transformadores [v1, v2, Rm, Xm]
Transformador1 = [69000, 500000, 4320, 5050]
Transformador2 = [500000, 230000, 432000, 505000]
Transformador3 = [230000, 69000, 402000, 607000]

# Definição das cargas Z1, Z2 e Z3

Carga1 = (840 +(1j*np.pi*2*60*4.6))
Carga2 = (266.13 +(1j*np.pi*2*60*0.643))
Carga3 = (120.28 +(1j*np.pi*2*60*0.29))

#--->
#Carga1 = (1906.263 +(1j*359.091))

#Carga2 = (266.13 +(1j*50.13))

#Carga3 = (120.28 +(1j*22.70))

#<----


#CRIANDO OBJETOS DOS COMPONENTES DA LINHA
Zth = Carga(Zf) #Impedância série do começo da linha
T1 = Transformador(Transformador1) #Primeiro Transformador
T2 = Transformador(Transformador2) #Segundo Transformador
T3 = Transformador(Transformador3) #Terceiro Transformador
LT1 = M_Linha(LT1_LT2) #Primeira Linha
LT2 = M_Linha(LT1_LT2) #Segunda Linha
LT3 = M_Linha(LT3_) #Terceira Linha
LT4 = M_Linha(LT4_) #Quarta Linha
Z1 = CargaDerivada(Carga1) #Cargas 1, 2 e 3
Z2 = CargaDerivada(Carga2)
Z3 = CargaDerivada(Carga3)

#CRIANDO BLOCOS DAS LINHAS POR ASSOCIAÇÃO DE MATRIZES
matriz1 = Associar_Matriz_em_Cascata (Zth,T1)
parallel_LT1_2 = Matriz_Paralelo (LT1,LT2)
matriz2 = Associar_Matriz_em_Cascata (matriz1, parallel_LT1_2) #MatrizTransmissão até a carga 1

matriz3 = Associar_Matriz_em_Cascata (matriz2, Z1)
matriz4 = Associar_Matriz_em_Cascata (matriz3, LT3)
matriz5 = Associar_Matriz_em_Cascata (matriz4, T2) #Matriz Transmissão até a carga 2

matriz6 = Associar_Matriz_em_Cascata (matriz5, Z2)
matriz7 = Associar_Matriz_em_Cascata (matriz6, LT4)
matriz8 = Associar_Matriz_em_Cascata (matriz7, T3) #Matriz Transmissão atéa carga 3

Fim_da_Linha = Associar_Matriz_em_Cascata (matriz8, Z3)

print(Fim_da_Linha)

#----------------------Item 3------------------------------#
#Encontrando Vc e Ic em Z1
print("Tensão e corrente em rms Z1:")
Vg = 69e+3
Vc = (Vg*Carga1)/(matriz2[0][0]*Carga1 + matriz2[0][1])
Ic = Vc/ Carga1
#Potencia ativa e reativa
Ic_cong1 = np.conjugate(Ic)
s1 = Vc* Ic_cong1
P1 = s1.real
Q1 = s1.imag
absVc, angVc = cmath.polar(Vc)
absIc, angIc = cmath.polar(Ic)
angVc = np.rad2deg(angVc)
angIc = np.rad2deg(angIc)
print(f"Vcarga1 = {absVc:.2f} ∠ {angVc:.2f}º V")
print(f"Icarga1 = {absIc:.2f} ∠ {angIc:.2f}º A\n")

#Encontrando Vc e IC em Z2
print("Tensão e corrente em rms Z2:")
Vg = 69e+3
Vc = (Vg*Carga2)/(matriz5[0][0]*Carga2 + matriz5[0][1])
Ic = Vc/Carga2
#Potencia ativa e reativa
Ic_cong2 = np.conjugate(Ic)
s2 = Vc* Ic_cong2
P2 = s2.real
Q2 = s2.imag
absVc, angVc = cmath.polar(Vc)
absIc, angIc = cmath.polar(Ic)
angVc = np.rad2deg(angVc)
angIc = np.rad2deg(angIc)
print(f"Vcarga2 = {absVc:.2f} ∠ {angVc:.2f}º V")
print(f"Icarga2 = {absIc:.2f} ∠ {angIc:.2f}º A\n")

#----------------------Item 1------------------------------#
#Encontrando Vc e IC em Z3
print("Tensão e corrente em rms Z3:")
Vg = 69e+3
Vc = (Vg*Carga3)/(matriz8[0][0]*Carga3 + matriz8[0][1]) #Tensão na carga 3
Ic = Vc/Carga3 #Corrente na carga 3
#Potencia ativa e reativa
Ic_cong3 = np.conjugate(Ic)
s3 = Vc* Ic_cong3
P3 = s3.real
Q3 = s3.imag
absVc, angVc = cmath.polar(Vc)
absIc, angIc = cmath.polar(Ic)
angVc = np.rad2deg(angVc)
angIc = np.rad2deg(angIc)
print(f"Vcarga3 = {absVc:.2f} ∠ {angVc:.2f}º V")
print(f"Icarga3 = {absIc:.2f} ∠ {angIc:.2f}º A\n")

#----------------------Item 4------------------------------#
#Potencia ativa e reativa no gerador
Ic_g = Vg / matriz1[0][1]
# Potência aparente (S) no gerador
S_g = Vg * np.conjugate(Ic_g) # Produto da tensão do gerador e o conjugado da corrente

# Separando a potência ativa (P) e a potência reativa (Q) no gerador
P_g = S_g.real # Parte real de S_g é a potência ativa
Q_g = S_g.imag # Parte imaginária de S_g é a potência reativa

# Imprime os resultados
print(f"Potência ativa e reativa do gerador:")
print(f"Potência ativa (P_g) = {P_g:.2f} W")
print(f"Potência reativa (Q_g) = {Q_g:.2f} var\n")

#----------------------Item 5------------------------------#
print(f"Potência ativa e reativa de Z3: \nP_Z3 = {P3:.2f} W e Q_Z3 = {Q3:.2f} var\n")
print(f"Potência ativa e reativa de Z2: \nP_Z2 = {P2:.2f} W e Q_Z2 = {Q2:.2f} var\n")
print(f"Potência ativa e reativa de Z1: \nP_Z1 = {P1:.2f} W e Q_Z1 = {Q1:.2f} var\n")

#----------------------Item 6------------------------------#
#Calculo das perdas
print(f"Perda de potência ativa e reativa do sistema:")
P_cargas_total = P1+P2+P3
Perdas_ativa = P_g - P_cargas_total
print(f"Perda da potência ativa = {Perdas_ativa:.2f} W")
Q_cagas_total = Q1+Q2+Q3
Perdas_reativa = Q_g - Q_cagas_total
print(f"Perda da potência Reativa = {Perdas_reativa:.2f} Var")
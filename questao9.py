import numpy as np
import cmath

# Definição das distâncias das linhas de transmissão (em km)
DISTANCIA_LT1 = DISTANCIA_LT2 = 100  # Distância das linhas 1 e 2
DISTANCIA_LT3 = 80                   # Distância da linha 3
DISTANCIA_LT4 = 120                  # Distância da linha 4

# Definição dos transformadores [tensão_primária, tensão_secundária, resistência_magnetização, reatância_magnetização]
TRANSFORMADOR_1 = [69000, 502311, 4320, 13.3955]    # Transformador 1
TRANSFORMADOR_2 = [500000, 230000, 432000, 1339.55] # Transformador 2
TRANSFORMADOR_3 = [230000, 69000, 402000, 1610.12]   # Transformador 3

# Impedância série Thevenin (em ohms)
Z_THEVENIN = 2 + 1j * 0.38  # Impedância equivalente do sistema

# Frequência angular do sistema (rad/s)
FREQUENCIA_ANGULAR = 2 * np.pi * 60  # Frequência de 60 Hz convertida para rad/s

# Definição das impedâncias das cargas (em ohms)
CARGA_1 = 840 + 1j * FREQUENCIA_ANGULAR * 4.6    # Impedância da carga 1
CARGA_2 = 117.55 + 1j * FREQUENCIA_ANGULAR * 0.643  # Impedância da carga 2
CARGA_3 = 52.9 + 1j * FREQUENCIA_ANGULAR * 0.29    # Impedância da carga 3

# Parâmetros das linhas de transmissão [resistência, indutância, capacitância] por km
PARAMETROS_LT1_LT2 = [DISTANCIA_LT1 * 0.182, DISTANCIA_LT1 * 2.28e-3, DISTANCIA_LT1 * 0.0140e-6]  # Linhas 1 e 2
PARAMETROS_LT3 = [DISTANCIA_LT3 * 0.182, DISTANCIA_LT3 * 2.28e-3, DISTANCIA_LT3 * 0.0140e-6]      # Linha 3
PARAMETROS_LT4 = [DISTANCIA_LT4 * 0.182, DISTANCIA_LT4 * 2.28e-3, DISTANCIA_LT4 * 0.0140e-6]      # Linha 4

# Função para calcular a matriz de transmissão de uma linha
def matriz_transmissao_linha(parametros):
    resistencia, indutancia, capacitancia = parametros
    admitancia = 1j * FREQUENCIA_ANGULAR * (capacitancia / 2)
    impedancia = resistencia + 1j * FREQUENCIA_ANGULAR * indutancia
    A = 1 + (admitancia * impedancia)
    B = impedancia
    C = (2 * admitancia) + (admitancia * admitancia * impedancia)
    D = 1 + (admitancia * impedancia)
    return np.array([[A, B], [C, D]], dtype=complex)

# Função para calcular a matriz de transmissão de uma carga em série
def matriz_carga_serie(impedancia):
    return np.array([[1, impedancia], [0, 1]], dtype=complex)

# Função para calcular a matriz de transmissão de uma carga em derivação
def matriz_carga_derivacao(impedancia):
    return np.array([[1, 0], [1 / impedancia, 1]], dtype=complex)

# Função para calcular a matriz de transmissão de um transformador
def matriz_transmissao_transformador(transformador):
    tensao_primaria, tensao_secundaria, resistencia_mag, reatancia_mag = transformador
    z1 = 7.6e-3 + 1j * 3.8e-3
    z2 = 33.9e-3 + 1j * 0.85e-3
    admitancia = (resistencia_mag + 1j * reatancia_mag) / (resistencia_mag * 1j * reatancia_mag)
    relacao_transformacao = tensao_primaria / tensao_secundaria
    A = relacao_transformacao * (1 + (admitancia * z1))
    B = (z1 + z2 + admitancia * z1 * z2) / relacao_transformacao
    C = relacao_transformacao * admitancia
    D = (1 + admitancia * z2) / relacao_transformacao
    return np.array([[A, B], [C, D]], dtype=complex)

# Função para multiplicar matrizes em cascata
def multiplicar_matrizes_cascata(matriz1, matriz2):
    return np.dot(matriz1, matriz2)

# Função para calcular a matriz de transmissão de elementos em paralelo
def matriz_transmissao_paralelo(matriz1, matriz2):
    Aa, Ba, Ca, Da = matriz1[0][0], matriz1[0][1], matriz1[1][0], matriz1[1][1]
    Ab, Bb, Cb, Db = matriz2[0][0], matriz2[0][1], matriz2[1][0], matriz2[1][1]
    try:
        A = (Aa * Bb + Ab * Ba) / (Ba + Bb)
        B = (Ba * Bb) / (Ba + Bb)
        C = Ca + Cb + ((Aa - Ab) * (Db - Da)) / (Ba + Bb)
        D = (Bb * Da + Ba * Db) / (Ba + Bb)
    except ZeroDivisionError:
        raise ValueError("Divisão por zero na matriz de transmissão em paralelo.")
    return np.array([[A, B], [C, D]], dtype=complex)

# Criando objetos dos componentes do sistema
Z_TH = matriz_carga_serie(Z_THEVENIN)
T1 = matriz_transmissao_transformador(TRANSFORMADOR_1)
T2 = matriz_transmissao_transformador(TRANSFORMADOR_2)
T3 = matriz_transmissao_transformador(TRANSFORMADOR_3)
LT1 = matriz_transmissao_linha(PARAMETROS_LT1_LT2)
LT2 = matriz_transmissao_linha(PARAMETROS_LT1_LT2)
LT3 = matriz_transmissao_linha(PARAMETROS_LT3)
LT4 = matriz_transmissao_linha(PARAMETROS_LT4)
Z1 = matriz_carga_derivacao(CARGA_1)
Z2 = matriz_carga_derivacao(CARGA_2)
Z3 = matriz_carga_derivacao(CARGA_3)

# Calculando blocos das linhas por associação de matrizes
MATRIZ_1 = multiplicar_matrizes_cascata(Z_TH, T1)  # Até o transformador 1
PARALELO_LT1_LT2 = matriz_transmissao_paralelo(LT1, LT2)  # LT1 e LT2 em paralelo
SERIE_PARALELO_LT3 = multiplicar_matrizes_cascata(PARALELO_LT1_LT2, LT3)  # (LT1 || LT2) em série com LT3
MATRIZ_2 = multiplicar_matrizes_cascata(MATRIZ_1, SERIE_PARALELO_LT3)  # Até antes da carga 1
MATRIZ_3 = multiplicar_matrizes_cascata(MATRIZ_2, Z1)  # Incluindo a carga 1
MATRIZ_4 = multiplicar_matrizes_cascata(MATRIZ_3, LT4)  # Até a linha 4
MATRIZ_5 = multiplicar_matrizes_cascata(MATRIZ_4, T2)  # Até antes da carga 2
MATRIZ_6 = multiplicar_matrizes_cascata(MATRIZ_5, Z2)  # Incluindo a carga 2
MATRIZ_7 = multiplicar_matrizes_cascata(MATRIZ_6, T3)  # Até antes da carga 3
MATRIZ_FINAL = multiplicar_matrizes_cascata(MATRIZ_7, Z3)  # Sistema completo

# Função para formatar valores com unidades apropriadas
def formatar_valor(valor, unidade_base, precisao=2):
    prefixos = [(1e9, 'G'), (1e6, 'M'), (1e3, 'k'), (1, ''), (1e-3, 'm'), (1e-6, 'µ'), (1e-9, 'n')]
    for fator, prefixo in prefixos:
        if abs(valor) >= fator or fator == 1:
            valor_formatado = valor / fator
            return f"{valor_formatado:.{precisao}f} {prefixo}{unidade_base}"
    return f"{valor:.{precisao}f} {unidade_base}"

# Função para calcular tensão e corrente usando matriz ABCD
def calcular_tensao_corrente(matriz, tensao_gerador, impedancia_carga):
    A, B = matriz[0][0], matriz[0][1]
    try:
        tensao_carga = tensao_gerador / (A + B / impedancia_carga)
        corrente_carga = tensao_carga / impedancia_carga
    except ZeroDivisionError:
        raise ValueError("Divisão por zero no cálculo de tensão ou corrente.")
    return tensao_carga, corrente_carga

# Função para calcular potências
def calcular_potencias(tensao, corrente):
    corrente_conjugada = np.conjugate(corrente)
    potencia_aparente = tensao * corrente_conjugada
    return potencia_aparente.real, potencia_aparente.imag

# Função para exibir resultados
def exibir_resultados(tensao, corrente, potencia_ativa, potencia_reativa, nome_carga):
    print(f"{'='*50}")
    print(f"Resultados para a carga {nome_carga}:")
    print(f"{'-'*50}")
    print(f"Potência ativa: {formatar_valor(potencia_ativa, 'W')}")
    print(f"Potência reativa: {formatar_valor(potencia_reativa, 'var')}")
    print("Formato retangular:")
    print(f"  Tensão: {formatar_valor(tensao.real, 'V')} + j{formatar_valor(tensao.imag, 'V')}")
    print(f"  Corrente: {formatar_valor(corrente.real, 'A')} + j{formatar_valor(corrente.imag, 'A')}")
    print("Formato polar:")
    abs_tensao, ang_tensao = cmath.polar(tensao)
    abs_corrente, ang_corrente = cmath.polar(corrente)
    ang_tensao, ang_corrente = np.rad2deg(ang_tensao), np.rad2deg(ang_corrente)
    print(f"  Tensão: {formatar_valor(abs_tensao, 'V')} ∠ {ang_tensao:.2f}°")
    print(f"  Corrente: {formatar_valor(abs_corrente, 'A')} ∠ {ang_corrente:.2f}°")
    print(f"{'='*50}\n")

# Tensão do gerador
TENSAO_GERADOR = 69e3  # 69 kV RMS

# Calculando para as cargas
tensao_z1, corrente_z1 = calcular_tensao_corrente(MATRIZ_2, TENSAO_GERADOR, CARGA_1)
potencia_ativa_z1, potencia_reativa_z1 = calcular_potencias(tensao_z1, corrente_z1)
exibir_resultados(tensao_z1, corrente_z1, potencia_ativa_z1, potencia_reativa_z1, "Z1")

tensao_z2, corrente_z2 = calcular_tensao_corrente(MATRIZ_5, TENSAO_GERADOR, CARGA_2)
potencia_ativa_z2, potencia_reativa_z2 = calcular_potencias(tensao_z2, corrente_z2)
exibir_resultados(tensao_z2, corrente_z2, potencia_ativa_z2, potencia_reativa_z2, "Z2")

tensao_z3, corrente_z3 = calcular_tensao_corrente(MATRIZ_7, TENSAO_GERADOR, CARGA_3)
potencia_ativa_z3, potencia_reativa_z3 = calcular_potencias(tensao_z3, corrente_z3)
exibir_resultados(tensao_z3, corrente_z3, potencia_ativa_z3, potencia_reativa_z3, "Z3")

# Calculando potência no gerador
A1, B1 = MATRIZ_1[0][0], MATRIZ_1[0][1]
try:
    corrente_gerador = TENSAO_GERADOR / B1
except ZeroDivisionError:
    raise ValueError("Divisão por zero no cálculo da corrente do gerador.")
potencia_aparente_gerador = TENSAO_GERADOR * np.conjugate(corrente_gerador)
potencia_ativa_gerador = potencia_aparente_gerador.real
potencia_reativa_gerador = potencia_aparente_gerador.imag
print(f"{'='*50}")
print("Potência ativa e reativa do gerador:")
print(f"{'-'*50}")
print(f"Potência ativa: {formatar_valor(potencia_ativa_gerador, 'W')}")
print(f"Potência reativa: {formatar_valor(potencia_reativa_gerador, 'var')}")
print(f"{'='*50}\n")

# Calculando perdas do sistema
print(f"{'='*50}")
print("Perda de potência ativa e reativa do sistema:")
print(f"{'-'*50}")
potencia_ativa_total_cargas = potencia_ativa_z1 + potencia_ativa_z2 + potencia_ativa_z3
perdas_ativa = potencia_ativa_gerador - potencia_ativa_total_cargas
potencia_reativa_total_cargas = potencia_reativa_z1 + potencia_reativa_z2 + potencia_reativa_z3
perdas_reativa = potencia_reativa_gerador - potencia_reativa_total_cargas
print(f"Perda de potência ativa: {formatar_valor(perdas_ativa, 'W')}")
print(f"Perda de potência reativa: {formatar_valor(perdas_reativa, 'var')}")
print(f"{'='*50}")
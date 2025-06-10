import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# --- 1. PARÁMETROS DE ENTRADA (MODIFICA ESTOS VALORES SEGÚN TU POZO) ---

# Geometría del Pozo
ID_CASING_IN = 6.625
OD_TUBING_IN = 4.5
H_FUGA_M = 2307
RUGOSIDAD_TUBERIA_M = 0.000045 # Rugosidad de la tubería (acero comercial) en metros

# --- Parámetros del Casing ---
PROFUNDIDAD_GLC_INICIAL_CASING_M = 2370 # 3000*50/100
RHO_LIQUIDO_CASING_KG_M3 = 1074
RHO_GAS_CASING_KG_M3 = 0.725
VISCOSIDAD_LIQUIDO_PA_S = 0.01 # Viscosidad del líquido en Pascal-segundo (10 cP)
PROFUNDIDAD_FLUIDO_ANULAR_INICIAL_M = 0 # m

# --- Propiedades del Fluido Anular ---
RHO_ANULAR_KG_M3 = 720
V_GAS_ANULAR_INICIAL_M3 = 1

# --- Propiedades del Gas (para Factor Z) ---
GRAVEDAD_ESPECIFICA_GAS = 0.591
TEMPERATURA_SISTEMA_K = 353.15 # Temperatura promedio del sistema en Kelvin

# --- Condiciones Iniciales ---
P_SUP_CASING_BAR = 39
P_SUP_ANULAR_BAR = 0

# --- Parámetros de la Fuga ---
DIAMETRO_FUGA_IN = 3/25.4
CD = 0.6

# --- Parámetros de la Simulación ---
DELTA_T_S = 5
MAX_T_MIN = 10000

# --- CONSTANTES ---
G = 9.81
M_PER_IN = 0.0254
PA_PER_BAR = 100000
R = 8.314

# --- 2. FUNCIONES AUXILIARES AVANZADAS (Sin cambios) ---
def calcular_propiedades_pseudocriticas(sg):
    ppc_pa = (4.892 - 0.4048 * sg) * 1e6
    tpc_k = (94.72 + 170.75 * sg)
    return ppc_pa, tpc_k

def calcular_z_factor_papay(ppr, tpr):
    z = 1 - (3.53 * ppr) / (10**(0.9813 * tpr)) + (0.274 * ppr**2) / (10**(0.8157 * tpr))
    return z

def calcular_factor_friccion_churchill(re, d_h, rugosidad):
    if re <= 2300: return 64 / re
    A = (-2.457 * np.log((7/re)**0.9 + 0.27 * (rugosidad/d_h)))**16
    B = (37530/re)**16
    f = 8 * ((8/re)**12 + 1/((A + B)**1.5))**(1/12)
    return f

# --- 3. CONVERSIONES Y CÁLCULOS INICIALES ---
id_casing_m = ID_CASING_IN * M_PER_IN
od_tubing_m = OD_TUBING_IN * M_PER_IN
diametro_fuga_m = DIAMETRO_FUGA_IN * M_PER_IN
p_sup_casing_inicial_pa = P_SUP_CASING_BAR * PA_PER_BAR
p_sup_anular_inicial_pa = P_SUP_ANULAR_BAR * PA_PER_BAR
a_anular_m2 = (np.pi / 4) * (id_casing_m**2 - od_tubing_m**2)
a_casing_m2 = (np.pi / 4) * id_casing_m**2
a_fuga_m2 = (np.pi / 4) * (diametro_fuga_m**2)
d_hidraulico_anular_m = id_casing_m - od_tubing_m
ppc_anular_pa, tpc_anular_k = calcular_propiedades_pseudocriticas(GRAVEDAD_ESPECIFICA_GAS)
tpr_anular = TEMPERATURA_SISTEMA_K / tpc_anular_k

# Presión Casing
h_gas_casing_ini = PROFUNDIDAD_GLC_INICIAL_CASING_M
h_liq_casing_ini = H_FUGA_M - h_gas_casing_ini
p_hidro_gas_casing_ini = RHO_GAS_CASING_KG_M3 * G * h_gas_casing_ini
p_hidro_liq_casing_ini = RHO_LIQUIDO_CASING_KG_M3 * G * h_liq_casing_ini
p_casing_inicial_total = p_sup_casing_inicial_pa + p_hidro_gas_casing_ini + p_hidro_liq_casing_ini
print(f"Cálculo de Presión Casing a {H_FUGA_M} m:")
print(f"  P. Superficie: {p_sup_casing_inicial_pa/PA_PER_BAR:.2f} bar")
print(f"  P. Hidro. Gas ({h_gas_casing_ini} m): {p_hidro_gas_casing_ini/PA_PER_BAR:.2f} bar")
print(f"  P. Hidro. Líquido ({h_liq_casing_ini} m): {p_hidro_liq_casing_ini/PA_PER_BAR:.2f} bar")
print(f"  ------------------------------------------")
print(f"  PRESIÓN TOTAL CASING: {p_casing_inicial_total/PA_PER_BAR:.2f} bar\n")

# Presión Anular
h_col_anular_ini = max(0, H_FUGA_M - PROFUNDIDAD_FLUIDO_ANULAR_INICIAL_M)
p_hidro_anular_ini = RHO_ANULAR_KG_M3 * G * h_col_anular_ini
p_anular_inicial_total = p_sup_anular_inicial_pa + p_hidro_anular_ini
print(f"Cálculo de Presión Anular a {H_FUGA_M} m:")
print(f"  P. Superficie: {p_sup_anular_inicial_pa/PA_PER_BAR:.2f} bar")
print(f"  P. Hidro. Anular ({h_col_anular_ini} m): {p_hidro_anular_ini/PA_PER_BAR:.2f} bar")
print(f"  ------------------------------------------")
print(f"  PRESIÓN TOTAL ANULAR: {p_anular_inicial_total/PA_PER_BAR:.2f} bar\n")

# Conclusión
delta_p_inicial = p_casing_inicial_total - p_anular_inicial_total
print("--- CONCLUSIÓN INICIAL ---")
print(f"Delta P (Casing - Anular) = {delta_p_inicial/PA_PER_BAR:.2f} bar")
if delta_p_inicial <= 0:
    print("El Delta P inicial es negativo o cero. La fuga del casing al anular es físicamente imposible.")
    print("La simulación no se iniciará.")
    sys.exit() # Termina el script
print("El Delta P inicial es positivo. Iniciando simulación dinámica...")
print("-" * 30)

# <<< NUEVO: Validación de la altura del fluido anular >>>
if PROFUNDIDAD_FLUIDO_ANULAR_INICIAL_M > H_FUGA_M:
    print(f"ADVERTENCIA: La fuga ({H_FUGA_M} m) está por encima del nivel de fluido anular inicial ({PROFUNDIDAD_FLUIDO_ANULAR_INICIAL_M} m).")

# --- 4. INICIALIZACIÓN DE LA SIMULACIÓN ---
resultados = pd.DataFrame(columns=[
    'Tiempo (s)', 'Delta_P (Pa)', 'Caudal (m³/s)', 'Velocidad (m/s)', 
    'Distancia Propagada (m)', 'P_sup_anular (bar)', 'P_casing_fuga (bar)'
])
t = 0.0
v_fugado_total_m3 = 0.0
max_t_s = MAX_T_MIN * 60
p_sup_anular_actual_pa = p_sup_anular_inicial_pa
profundidad_glc_casing_actual_m = PROFUNDIDAD_GLC_INICIAL_CASING_M

print("Iniciando simulación UNIVERSAL...")
# Informar al usuario qué modelo de densidad se está utilizando
if RHO_LIQUIDO_CASING_KG_M3 < RHO_ANULAR_KG_M3:
    print("-> Modelo activado: Fluido de fuga LIGERO (flota).")
else:
    print("-> Modelo activado: Fluido de fuga DENSO (se hunde).")
print("-" * 30)


# --- 5. BUCLE DE SIMULACIÓN ITERATIVA ---
while t < max_t_s:
    # --- A. CÁLCULO DE PRESIÓN DEL CASING (CON DEPLECIÓN) ---
    h_gas_casing_m = profundidad_glc_casing_actual_m
    h_liquido_casing_m = H_FUGA_M - h_gas_casing_m
    if h_liquido_casing_m < 0: h_liquido_casing_m = 0
    p_hidro_gas_casing = RHO_GAS_CASING_KG_M3 * G * h_gas_casing_m
    p_hidro_liquido_casing = RHO_LIQUIDO_CASING_KG_M3 * G * h_liquido_casing_m
    p_casing_a_fuga_pa = p_sup_casing_inicial_pa + p_hidro_gas_casing + p_hidro_liquido_casing
    
    # --- B. CÁLCULO DE PRESIÓN DEL ANULAR ---
    # B1. Presión en superficie con Z-Factor
    if V_GAS_ANULAR_INICIAL_M3 > 0:
        vol_gas_actual = V_GAS_ANULAR_INICIAL_M3 - v_fugado_total_m3
        if vol_gas_actual > 1e-9:
            p_sup_guess = p_sup_anular_actual_pa
            for _ in range(5):
                ppr_guess = p_sup_guess / ppc_anular_pa
                z_guess = calcular_z_factor_papay(ppr_guess, tpr_anular)
                ppr_inicial = p_sup_anular_inicial_pa / ppc_anular_pa
                z_inicial = calcular_z_factor_papay(ppr_inicial, tpr_anular)
                nRT = (p_sup_anular_inicial_pa * V_GAS_ANULAR_INICIAL_M3) / z_inicial
                p_sup_guess = (z_guess * nRT) / vol_gas_actual
            p_sup_anular_actual_pa = p_sup_guess
        else:
            print("Volumen de gas anular completamente comprimido. Deteniendo."); break

    # B2. Presión hidrostática (fuerzas corporales)
    h_fugado_m = v_fugado_total_m3 / a_anular_m2

   # <<< INICIO DE LÓGICA CONDICIONAL DE DENSIDAD ACTUALIZADA >>>
    if RHO_LIQUIDO_CASING_KG_M3 < RHO_ANULAR_KG_M3:
        # CASO 1: Fluido de fuga LIGERO (flota)
        h_columna_anular_original = max(0, H_FUGA_M - PROFUNDIDAD_FLUIDO_ANULAR_INICIAL_M)
        p_hidro_anular_pa = (RHO_ANULAR_KG_M3 * G * h_columna_anular_original) + (RHO_LIQUIDO_CASING_KG_M3 * G * h_fugado_m)
    else:
        # CASO 2: Fluido de fuga DENSO (se hunde)
        h_columna_anular_original = max(0, H_FUGA_M - PROFUNDIDAD_FLUIDO_ANULAR_INICIAL_M)
        p_hidro_anular_pa = RHO_ANULAR_KG_M3 * G * h_columna_anular_original
    # <<< FIN DE LÓGICA CONDICIONAL >>>


    # B3. Caída de presión por fricción
    delta_p_friccion = 0
    if not resultados.empty:
        velocidad_anterior = resultados.iloc[-1]['Velocidad (m/s)']
        if velocidad_anterior > 1e-6:
             reynolds = (RHO_LIQUIDO_CASING_KG_M3 * velocidad_anterior * d_hidraulico_anular_m) / VISCOSIDAD_LIQUIDO_PA_S
             f_factor = calcular_factor_friccion_churchill(reynolds, d_hidraulico_anular_m, RUGOSIDAD_TUBERIA_M)
             delta_p_friccion = f_factor * (h_fugado_m / d_hidraulico_anular_m) * (RHO_LIQUIDO_CASING_KG_M3 * velocidad_anterior**2 / 2)
    
    p_anular_a_fuga_pa = p_sup_anular_actual_pa + p_hidro_anular_pa + delta_p_friccion
    
    # --- C. CÁLCULO DEL CAUDAL Y ACTUALIZACIÓN ---
    delta_p = p_casing_a_fuga_pa - p_anular_a_fuga_pa
    if delta_p <= 0:
        print(f"\nEquilibrio alcanzado en t = {t/60:.2f} minutos."); break
        
    caudal_m3_s = CD * a_fuga_m2 * np.sqrt(max(0, 2 * delta_p / RHO_LIQUIDO_CASING_KG_M3))
    velocidad_m_s = caudal_m3_s / a_anular_m2
    delta_v_m3 = caudal_m3_s * DELTA_T_S
    
    # Almacenar resultados
    nueva_fila = {
        'Tiempo (s)': t, 'Delta_P (Pa)': delta_p, 'Caudal (m³/s)': caudal_m3_s,
        'Velocidad (m/s)': velocidad_m_s, 'Distancia Propagada (m)': h_fugado_m,
        'P_sup_anular (bar)': p_sup_anular_actual_pa / PA_PER_BAR,
        'P_casing_fuga (bar)': p_casing_a_fuga_pa / PA_PER_BAR
    }
    resultados = pd.concat([resultados, pd.DataFrame([nueva_fila])], ignore_index=True)
    
    # --- D. ACTUALIZACIÓN DEL ESTADO PARA LA SIGUIENTE ITERACIÓN ---
    v_fugado_total_m3 += delta_v_m3
    delta_h_casing = delta_v_m3 / a_casing_m2
    profundidad_glc_casing_actual_m += delta_h_casing
    t += DELTA_T_S

if t >= max_t_s: print(f"\nTiempo máximo de simulación ({MAX_T_MIN} min) alcanzado.")
print("\n--- Resultados Finales ---")
if not resultados.empty:
    print(f"Volumen total fugado: {v_fugado_total_m3:.2f} m³ ({v_fugado_total_m3 * 6.2898:.2f} barriles)")
    print(f"Presión final en sup. anular: {resultados.iloc[-1]['P_sup_anular (bar)']:.2f} bar")
    print(f"Altura final del fluido fugado en anular: {resultados.iloc[-1]['Distancia Propagada (m)']:.2f} m")

# --- 6. VISUALIZACIÓN DE RESULTADOS (Sin cambios) ---
if not resultados.empty:
    fig, axs = plt.subplots(4, 1, figsize=(12, 18), sharex=True)
    fig.suptitle('Evolución de la Fuga (Modelo Universal Avanzado)', fontsize=16)
    axs[0].plot(resultados['Tiempo (s)'], resultados['Caudal (m³/s)'], 'b-')
    axs[0].set_ylabel('Caudal (m³/s)'); axs[0].set_title('Caudal de Fuga'); axs[0].grid(True, ls='--', alpha=0.6)
    axs[1].plot(resultados['Tiempo (s)'], resultados['Velocidad (m/s)'], 'g-')
    axs[1].set_ylabel('Velocidad (m/s)'); axs[1].set_title('Velocidad de Propagación'); axs[1].grid(True, ls='--', alpha=0.6)
    axs[2].plot(resultados['Tiempo (s)'], resultados['Distancia Propagada (m)'], 'r-')
    axs[2].set_ylabel('Distancia (m)'); axs[2].set_title('Distancia de Propagación'); axs[2].grid(True, ls='--', alpha=0.6)
    axs[3].plot(resultados['Tiempo (s)'], resultados['P_casing_fuga (bar)'], 'k-', label='Presión Casing en Fuga')
    axs[3].plot(resultados['Tiempo (s)'], resultados['P_sup_anular (bar)'], 'm--', label='Presión Sup. Anular')
    axs[3].set_xlabel('Tiempo (s)'); axs[3].set_ylabel('Presión (bar)'); axs[3].set_title('Evolución de Presiones'); axs[3].grid(True, ls='--', alpha=0.6)
    axs[3].legend()
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.show()

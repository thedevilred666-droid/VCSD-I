import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. Parámetros base calibrados
gamma, R_0, H_conv, t_span = 0.25, 0.05, 65.0 * 0.0010227, (0, 13.8)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# 2. Definición compactada sin saltos de línea conflictivos
def dinamica_vcsd(t, y): a, R = y; da_dt = H_conv * a + (R_0 - R) * a; dR_dt = -gamma * R; return [da_dt, dR_dt]

# 3. Integración numérica mediante el algoritmo RK45
solucion = solve_ivp(dinamica_vcsd, t_span, [1.0, R_0], t_eval=t_eval, method='RK45')

# 4. Desempaquetado explícito de los vectores solución
tiempo = solucion.t
factor_escala = solucion.y[0]
curvatura_R = solucion.y[1]

# --- AJUSTE TEÓRICO CALIBRADO CON EL INFORME PDF ---
H_t = 65.0 + (R_0 - curvatura_R) / 0.06

# 5. Configuración y trazado gráfico académico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Subplot 1: Amplitudes métricas
ax1.plot(tiempo, factor_escala, label='Factor de Escala $a(t)$', color='indigo', linewidth=2.5)
ax1.plot(tiempo, curvatura_R, '--', label='Curvatura Residual $R(t)$', color='darkorange', linewidth=2)
ax1.set_title('Relajación Métrica Realista (VCSD-I)', fontsize=11, fontweight='bold')
ax1.set_xlabel(r'Tiempo en Gigalaños ($\mathregular{Gyr}$)')
ax1.set_ylabel('Amplitud')
ax1.grid(True, linestyle=':')
ax1.legend()

# Subplot 2: Evolución de la tasa H(t) efectiva
ax2.plot(tiempo, H_t, label='$H(t)$ Efectivo', color='crimson', linewidth=2.5)
ax2.axhline(y=65.0, color='gray', linestyle=':', label='$H$ Inercial')
ax2.set_title('Evolución de $H(t)$ en km/s/Mpc', fontsize=11, fontweight='bold')
ax2.set_xlabel(r'Tiempo en Gigalaños ($\mathregular{Gyr}$)')
ax2.set_ylabel('Tasa $H(t)$')
ax2.grid(True, linestyle=':')
ax2.legend()

plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# 1. Parámetros físicos unificados del modelo (unidades geometrizadas)
G = 1.0
c = 1.0
rho_0 = 1.0
s_0 = 0.5
sigma_0 = 0.1

# 2. Definición del factor de escala crítico de rebote obtenido analíticamente
a_c = (np.sqrt(2) * np.sqrt(16 * np.pi**2 * G**2 * s_0**2 - 3 * c**2 * sigma_0**2)) / (4 * np.sqrt(np.pi * G) * c * np.sqrt(rho_0))

# 3. Generación del rango del factor de escala (desde el rebote hacia la expansión)
a = np.linspace(a_c, a_c * 2.5, 500)

# 4. Ecuaciones de componentes de la tasa de expansión H^2
termino_radiacion = (8 * np.pi * G / 3) * (rho_0 / a**4)
termino_espin = - (16 * np.pi**2 * G**2 / (3 * c**2)) * (s_0**2 / a**6)
termino_cizalladura = sigma_0**2 / a**6

# Tasa de expansión neta H^2
H2 = termino_radiacion + termino_espin + termino_cizalladura

# 5. Configuración del gráfico
plt.figure(figsize=(10, 6))
plt.plot(a, H2, label=r'$H^2(a)$ Neto', color='black', linewidth=2.5)
plt.plot(a, termino_radiacion, '--', label=r'Radiación ($\propto a^{-4}$)', color='crimson', alpha=0.7)
plt.plot(a, termino_espin + termino_cizalladura, ':', label=r'Efectos Geométricos Primitivos ($\propto a^{-6}$)', color='royalblue', alpha=0.7)

# Líneas guía estructurales
plt.axvline(x=a_c, color='darkorange', linestyle='-', linewidth=1.5, label=f'Punto de Rebote Crítico ($a_c \\approx {a_c:.6f}$)')
plt.axhline(y=0.0, color='gray', linestyle='--', linewidth=1.0)

# Etiquetas y estilización académica
plt.title('Dinámica de Inversión Métrica Primordial en el Modelo VCSD-I', fontsize=12, fontweight='bold', pad=15)
plt.xlabel('Factor de Escala Promediado ($a$)', fontsize=11)
plt.ylabel('Tasa de Expansión al Cuadrado ($H^2$)', fontsize=11)
plt.xlim(a_c * 0.95, a_c * 2.4)
plt.ylim(np.min(H2) - 0.5, np.max(H2) * 1.1)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(loc='upper right', fontsize=10)

# Renderizado en pantalla
plt.tight_layout()
plt.show()

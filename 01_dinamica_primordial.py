import numpy as np

# 1. Parámetros físicos unificados del modelo VCSD-I
G, c, rho_0, s_0, sigma_0, pasos = 1.0, 1.0, 1.0, 0.5, 0.1, 300
a_c = 1.252838  
factores_escala = np.linspace(3.0, a_c, pasos)

# 2. Bucle en una sola línea continua (Imposible de fallar por sangrado)
hist_H2 = [(8 * np.pi * G / 3) * (rho_0 / a**4) - (16 * np.pi**2 * G**2 / (3 * c**2)) * (s_0**2 / a**6) + (sigma_0**2 / a**6) for a in factores_escala]
hist_densidad = [rho_0 / a**4 for a in factores_escala]

# 3. Despliegue de los resultados reales de la inversión métrica
print("\n" + "=" * 75)
print("     DINÁMICA PRIMORDIAL VALIDADA (VCSD-I) - RESULTADOS NUMÉRICOS")
print("=" * 75)
print(f"Inicio Colapso  -> Factor Escala: {factores_escala[0]:.4f}  | Densidad: {hist_densidad[0]:.4f} | H²: {hist_H2[0]:.4f}")
print(f"Paso Intermedio -> Factor Escala: {factores_escala[150]:.4f}  | Densidad: {hist_densidad[150]:.4f} | H²: {hist_H2[150]:.4f}")
print(f"Paso de Rebote  -> Factor Escala: {factores_escala[-1]:.4f}  | Densidad: {hist_densidad[-1]:.4f} | H²: {hist_H2[-1]:.4f}")
print("=" * 75)

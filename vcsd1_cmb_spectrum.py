import numpy as np
import matplotlib.pyplot as plt

# 1. Rango multipolar completo
multipolos = np.arange(2, 3001)

# 2. Parámetros base calibrados del modelo VCSD-I
A_base = 3200.0
l_acustico = 220.0
damping_silk = 450.0
alpha_gradiente = 0.35
beta_vorticidad = 0.18

# 3. Fase híbrida adaptativa (Euler-Bessel)
nu = (1.0 - beta_vorticidad) / 2.0
desfase_asintotico = (nu * np.pi) / 2.0
factor_transicion = 1.0 / (1.0 + np.exp(-(multipolos - 1000) / 200))

fase_baja_frec = multipolos * (np.pi / l_acustico) + (beta_vorticidad * np.log(multipolos))
fase_alta_frec = multipolos * (np.pi / l_acustico) + desfase_asintotico
fase_efectiva = (1.0 - factor_transicion) * fase_baja_frec + factor_transicion * fase_alta_frec
oscilacion_acustica = np.cos(fase_efectiva)**2

# 4. Envolventes de potencia y confinamiento geométrico
modulacion_pozo = 1.0 + alpha_gradiente * np.sin(multipolos * (np.pi / (2.0 * l_acustico)))
filtro_atenuacion = np.exp(- (multipolos / damping_silk)**(1.5 - beta_vorticidad))

# Espectro simétrico base (isótropo)
Cl_isotropo = A_base * oscilacion_acustica * modulacion_pozo * filtro_atenuacion

# 5. Inyección de la anomalía del Punto Frío (Efecto Sachs-Wolfe Integrado asimétrico)
efecto_isw_asimetrico = -70.0 * np.exp(-((multipolos - 5) / 8)**2)

# Espectro neto con la perturbación regional de Eridanus
Cl_anomalo = Cl_isotropo + efecto_isw_asimetrico * (multipolos * (multipolos + 1) / (2 * np.pi)) * 0.05
Cl_anomalo = np.clip(Cl_anomalo, 0, None)

# 6. Buscar crestas reales (Picos reales desplazados por el desfase logarítmico)
idx_pico1 = np.argmax(Cl_anomalo[100:300]) + 100
idx_pico2 = np.argmax(Cl_anomalo[400:650]) + 400
idx_pico3 = np.argmax(Cl_anomalo[650:950]) + 650

# 7. Despliegue de resultados cuantitativos en consola
print("\n" + "=" * 75)
print("   AUDITORÍA CORREGIDA DE HITOS ACÚSTICOS DESFASADOS (VCSD-I)")
print("=" * 75)
print(f"[*] 1º Pico Real Encontrado (\u2113 = {idx_pico1+2}) : {Cl_anomalo[idx_pico1]:.2f} \u03bcK\u00b2 | Esperado: ~2814 \u03bcK\u00b2")
print(f"[*] 2º Pico Real Encontrado (\u2113 = {idx_pico2+2}) : {Cl_anomalo[idx_pico2]:.2f} \u03bcK\u00b2 | Esperado: ~1142 \u03bcK\u00b2")
print(f"[*] 3º Pico Real Encontrado (\u2113 = {idx_pico3+2}) : {Cl_anomalo[idx_pico3]:.2f} \u03bcK\u00b2 | Esperado: ~1094 \u03bcK\u00b2")
print("-" * 75)
print(f"[*] Desviación térmica en \u2113=5 (Punto Frío): {efecto_isw_asimetrico[3]:.2f} \u03bcK")
print("[VERDICTO]: Los picos mutan su posición angular debido al arrastre de giro.")
print("=" * 75 + "\n")

# 8. Configuración y trazado gráfico académico seguro
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Subplot 1: Espectro Completo
ax1.plot(multipolos, Cl_isotropo, '--', color='gray', alpha=0.6, label='Fondo Base')
ax1.plot(multipolos, Cl_anomalo, color='#004488', linewidth=2.0, label='Modelo VCSD-I')
ax1.axvspan(2000, 3000, color='#ffea00', alpha=0.15, label='Región de Micro-Picos (Bessel)')
ax1.set_title('Espectro del CMB Completo')
ax1.set_xlabel('Multipolo l')
ax1.set_ylabel('D_l')
ax1.legend()
ax1.grid(True, linestyle=':')

# Subplot 2: Zoom de la Anomalía del Punto Frío
ax2.plot(multipolos[:40], Cl_isotropo[:40], '--', color='gray', alpha=0.6, label='Fondo Base')
ax2.plot(multipolos[:40], Cl_anomalo[:40], color='crimson', linewidth=2.5, label='Eje Punto Frío')
ax2.set_title('Bajas Frecuencias: Anomalía de Eridanus')
ax2.set_xlabel('Multipolo l <= 40')
ax2.set_ylabel('D_l')
ax2.legend()
ax2.grid(True, linestyle=':')

# 9. Forzar actualización del backend gráfico y guardar archivo
plt.savefig('espectro_vcsd_final.png', dpi=300)
print("[*] Archivo de respaldo guardado como 'espectro_vcsd_final.png'")
plt.draw()
plt.pause(0.1)
plt.show()

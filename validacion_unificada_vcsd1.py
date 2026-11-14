import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import jv

# ==============================================================================
# 1. EXPANSIÓN TARDÍA POR ALISAMIENTO MÉTRICO DE VACÍOS (CINEMÁTICA EXACTA)
# ==============================================================================
c_si = 2.99792458e8
gamma, H_inercial, Delta_H_max, conv = 0.25, 65.0, 0.86, 0.0010227
edo_tardia = lambda t, y: [((H_inercial + Delta_H_max * (1.0 - y.item(1))) * conv) * y.item(0), -gamma * y.item(1)]
t_eval = np.linspace(0, 13.8, 500)
sol = solve_ivp(edo_tardia, (0, 13.8), [1.0, 1.0], t_eval=t_eval)
H_t = H_inercial + Delta_H_max * (1.0 - sol.y[1])
H_final = H_t.item(-1)

# ==============================================================================
# 2. CURVAS DE ROTACIÓN GALÁCTICAS (UMBRAL COSMOLÓGICO a0 = c*H0 / 2pi)
# ==============================================================================
G_astro = 4.3009e-6
a0_astro = ((c_si * (H_final * 1000 / 3.0857e22)) / (2 * np.pi)) * (3.0857e19 / 1e6)
radios = np.linspace(0.2, 40.0, 300)
masa_r = lambda M_tot, r_d: M_tot * (radios**2 / (radios + r_d)**2)
vel_rot = lambda M_tot, r_d: np.sqrt(np.sqrt((G_astro * masa_r(M_tot, r_d) / radios**2)**2 + (G_astro * masa_r(M_tot, r_d) / radios**2) * a0_astro) * radios)

# ==============================================================================
# 3. ESPECTRO ANGULAR DEL CMB (EULER-BESSEL + REFRACCIÓN DEL POZO REGIONAL)
# ==============================================================================
multipolos = np.arange(2, 1600)
beta, l_acustico, damping_silk = 0.18, 282.5, 1480.0
nu = (1.0 - beta) / 2.0
x = (multipolos * (np.pi / l_acustico)) * (1.0 + 0.565 * np.exp(-multipolos / 320.0))

Psi_pozo = -0.45 * np.exp(-multipolos / 350.0)
modo_temp = (x**nu) * jv(-nu, x) + Psi_pozo
modo_dopp = (x**nu) * jv(1.0 - nu, x)

entrada_horizonte = (x**4) / (1.8**4 + x**4)
meseta_sw = 0.12 / (1.0 + (x / 1.5)**4)
potencia_bessel = meseta_sw + entrada_horizonte * (modo_temp**2 + 0.580**2 * modo_dopp**2)
filtro_silk = np.exp(-(multipolos / damping_silk)**(1.5 - beta))
Dl = 2848.0 * (potencia_bessel * filtro_silk) / np.max(potencia_bessel * filtro_silk)

# ==============================================================================
# 4. GENERACIÓN DE GRÁFICA CONSOLIDADA (FIGURA 1 DEL DOCUMENTO MAESTRO)
# ==============================================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

axes[0].plot(t_eval, H_t, color='crimson', linewidth=2.2, label='H(t) Efectivo')
axes[0].axhline(H_inercial, color='gray', linestyle=':', label='H Inercial (65.0)')
axes[0].set_title('Expansión Tardía por Alisamiento')
axes[0].set_xlabel('Tiempo Cósmico [Gyr]')
axes[0].set_ylabel('H(t) [km/s/Mpc]')
axes[0].grid(True, linestyle=':')
axes[0].legend()

axes[1].plot(radios, vel_rot(5.0e8, 1.2), linewidth=2.0, label='Galaxia Enana')
axes[1].plot(radios, vel_rot(6.0e10, 3.0), linewidth=2.0, label='Vía Láctea')
axes[1].plot(radios, vel_rot(3.0e11, 4.5), linewidth=2.0, label='Superespiral')
axes[1].set_title('Curvas de Rotación (Umbral a_0)')
axes[1].set_xlabel('Radio Galáctico [kpc]')
axes[1].set_ylabel('Velocidad Tangencial [km/s]')
axes[1].set_ylim(0, 320)
axes[1].grid(True, linestyle=':')
axes[1].legend()

axes[2].plot(multipolos, Dl, color='#004488', linewidth=1.8, label='VCSD-I (Euler-Bessel)')
axes[2].set_title('Espectro Angular del CMB (Convención 4π)')
axes[2].set_xlabel('Multipolo Angular (l)')
axes[2].set_ylabel('D_l [μK²]')
axes[2].set_xlim(2, 1500)
axes[2].grid(True, linestyle=':')
axes[2].legend()

plt.tight_layout()
plt.savefig('validacion_unificada_vcsd1.png', dpi=300)
plt.show()

plt.show()

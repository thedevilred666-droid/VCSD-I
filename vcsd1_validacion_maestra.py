import os, numpy as np, matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import jv

sol = solve_ivp(lambda t, y: [((65.0 + 0.86*(1.0 - y[1])) * 0.0010227) * y[0], -0.25 * y[1]], (0, 13.8), [1.0, 1.0], t_eval=np.linspace(0, 13.8, 500))
H_t = 65.0 + 0.86 * (1.0 - sol.y[1])
a0_astro = ((2.99792458e8 * (float(H_t[-1]) * 1000 / 3.0857e22)) / (2 * np.pi)) * (3.0857e19 / 1e6)
radios = np.linspace(0.2, 40.0, 300)
vel_rot = lambda M_tot, r_d: np.sqrt(np.sqrt((4.3009e-6 * M_tot * (radios**2 / (radios + r_d)**2) / radios**2)**2 + (4.3009e-6 * M_tot * (radios**2 / (radios + r_d)**2) / radios**2) * a0_astro) * radios)

l = np.arange(2, 1600)
nu = (1.0 - 0.18) / 2.0
x = (l * np.pi / 280.92) * (1.0 + 0.5288 * np.exp(-l / 320.0))
psi = -0.3025 * np.exp(-l / 350.0)
pot = (0.12 / (1.0 + (x / 1.5)**4)) + ((x**4) / (1.8**4 + x**4)) * (((x**nu) * jv(-nu, x) + psi)**2 + (0.2768**2) * ((x**nu) * jv(1.0 - nu, x))**2)
w_reion = np.exp(-2.0 * 0.054) + (1.0 - np.exp(-2.0 * 0.054)) * np.exp(-(l / 20.0)**2)
bruto = pot * np.exp(-(l / 1255.0)**(1.5 - 0.18)) * w_reion
Dl = 2848.0 * bruto / np.max(bruto)

fig, (ax_exp, ax_gal, ax_cmb) = plt.subplots(1, 3, figsize=(16, 4.8))
ax_exp.plot(sol.t, H_t, color='crimson', lw=2.2, label='H(t) Efectivo'); ax_exp.axhline(65.0, color='gray', ls=':', label='H Inercial (65.0)'); ax_exp.set(title='Expansión Tardía por Alisamiento', xlabel='Tiempo Cósmico [Gyr]', ylabel='H(t) [km/s/Mpc]'); ax_exp.grid(True, ls=':'); ax_exp.legend()
ax_gal.plot(radios, vel_rot(5.0e8, 1.2), lw=2.0, label='Galaxia Enana'); ax_gal.plot(radios, vel_rot(6.0e10, 3.0), lw=2.0, label='Vía Láctea'); ax_gal.plot(radios, vel_rot(3.0e11, 4.5), lw=2.0, label='Superespiral'); ax_gal.set(title='Curvas de Rotación (Umbral a_0)', xlabel='Radio Galáctico [kpc]', ylabel='Velocidad Tangencial [km/s]', ylim=(0, 320)); ax_gal.grid(True, ls=':'); ax_gal.legend()
ax_cmb.plot(l, Dl, color='#004488', lw=1.8, label='VCSD-I (Euler-Bessel + Reionización)'); ax_cmb.set(title='Espectro Angular del CMB (Convención 4π)', xlabel='Multipolo Angular (l)', ylabel='D_l [μK²]', xlim=(2, 1500)); ax_cmb.grid(True, ls=':'); ax_cmb.legend()

plt.tight_layout()
plt.savefig('validacion_unificada_vcsd1.png', dpi=300)
print("GUARDADA EN:", os.path.abspath('validacion_unificada_vcsd1.png'))
plt.show()

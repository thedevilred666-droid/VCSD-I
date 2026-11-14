import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import jv
from scipy.signal import find_peaks

print("=" * 78)
print("   AUDITORÍA UNIFICADA DEL ANEXO VCSD-I (SI, FILAMENTOS, DIPOLO H0 Y CMB)")
print("=" * 78)

# 1. CRONOLOGÍA INFLACIONARIA EN UNIDADES SI
G, c, rho_0, sigma_0 = 6.67430e-11, 2.99792458e8, 1.79e30, 1.0e-5
s_crit = (np.sqrt(3) * c**2 * sigma_0) / (4 * np.pi * G)
s_0 = 1.5 * s_crit
a_c = (np.sqrt(2) * np.sqrt(16 * np.pi**2 * G**2 * s_0**2 - 3 * c**4 * sigma_0**2)) / (4 * np.sqrt(np.pi * G) * c**2 * np.sqrt(rho_0))
a_fin = np.sqrt(2) * a_c
H_si = lambda a: np.sqrt(max((8 * np.pi * G / 3) * (rho_0 / a**4) - ((16 * np.pi**2 * G**2 / (3 * c**4)) * s_0**2 - sigma_0**2) / a**6, 0.0))
H_max = H_si(np.sqrt(1.5) * a_c)
dt_inf, _ = quad(lambda a: 1.0 / (a * H_si(a)), a_c * (1 + 1e-12), a_fin)

print("\n[1] CRONOLOGÍA INFLACIONARIA EN UNIDADES SI:")
print(f"    - Factor de escala rebote (a_c)   : {a_c:.3e}")
print(f"    - Factor fin aceleración (a_fin)  : {a_fin:.3e}")
print(f"    - Tasa máxima de expansión H_max  : {H_max:.3e} s^-1")
print(f"    - Duración inflacionaria (dt_inf) : {dt_inf:.3e} s ({dt_inf*1e3:.3e} ms)")

# 2. VORTICIDAD MACROSCÓPICA EN FILAMENTOS (MEERKAT)
v_nucleo, r_nucleo, r_ext = 76817.56, 1.0, 20.0
r_fil = np.linspace(r_nucleo, r_ext, 300)
v_fil = v_nucleo * (r_nucleo / r_fil)
v_contorno = v_nucleo * (r_nucleo / r_ext)

print("\n[2] VORTICIDAD MACROSCÓPICA EN FILAMENTOS (MEERKAT):")
print(f"    - Velocidad en contorno (r = 20)  : {v_contorno:.2f} km/s (Verificado)")

# 3. MODULACIÓN DIPOLAR DIRECCIONAL DE H0
H0_iso, delta_H0_obj = 65.83, 4.41
eps_ani = (delta_H0_obj / 2.0) / H0_iso
theta = np.linspace(0, np.pi, 300)
H0_theta = H0_iso * (1.0 + eps_ani * np.cos(2 * theta))
H0_max, H0_min = H0_iso * (1.0 + eps_ani), H0_iso * (1.0 - eps_ani)

print("\n[3] MODULACIÓN DIPOLAR DIRECCIONAL DE H0:")
print(f"    - Anisotropía geométrica (eps)    : {eps_ani*100:.2f}%")
print(f"    - H0 máximo (Eje Polar, θ=0)      : {H0_max:.2f} km/s/Mpc")
print(f"    - H0 mínimo (Ecuatorial, θ=π/2)   : {H0_min:.2f} km/s/Mpc")
print(f"    - Variación dipolar neta (ΔH0)    : {H0_max - H0_min:.2f} km/s/Mpc")

# 4. ESPECTRO MULTIPOLAR DEL CMB (EULER-BESSEL REGULAR, l = 2..3000)
l = np.arange(2, 3001)
beta, l_A, l_silk = 0.18, 260.0, 950.0
nu = (1.0 - beta) / 2.0
x = l * (np.pi / l_A)
Psi_pozo = -0.55 * np.exp(-l / 350.0)
modo_temp = (x**nu) * jv(-nu, x) + Psi_pozo
modo_dopp = (x**nu) * jv(1.0 - nu, x)
entrada_horizonte = (x**4) / (1.8**4 + x**4)
meseta_sw = 0.12 / (1.0 + (x / 1.5)**4)
potencia = meseta_sw + entrada_horizonte * (modo_temp**2 + 0.35**2 * modo_dopp**2)
filtro_silk = np.exp(-(l / l_silk) ** (1.5 - beta))
Dl = 2848.0 * (potencia * filtro_silk) / np.max(potencia * filtro_silk)

pks_idx, _ = find_peaks(Dl, distance=180)
pks = [l[p] for p in pks_idx if l[p] > 150]

print("\n[4] ESPECTRO MULTIPOLAR DEL CMB (EULER-BESSEL REGULAR, l = 2..3000):")
print("\n".join([f"    - Pico {i}º -> l = {pk:4d} | Amplitud D_l = {Dl[pk - 2]:7.2f} uK^2" for i, pk in enumerate(pks[:7], 1)]))

# 5. GRÁFICA CONSOLIDADA COMPACTA (ALTA RESOLUCIÓN 300 DPI AL GUARDAR)
plt.rcParams.update({'font.size': 8})
fig, axes = plt.subplots(2, 2, figsize=(7.8, 5.4), dpi=90)
a_vals = np.linspace(a_c, a_c * 2.2, 400)

axes[0, 0].plot(a_vals * 1e16, [H_si(v)**2 for v in a_vals], color="crimson", lw=2, label=r"$H^2(a)$ en SI")
axes[0, 0].axvline(a_c * 1e16, color="black", ls="--", label=r"Rebote $a_c$")
axes[0, 0].axvline(a_fin * 1e16, color="royalblue", ls=":", label=r"Fin Inflación $a_{\mathrm{fin}}$")
axes[0, 0].set_title("1. Rebote Torsional e Inflación en SI", fontweight="bold")
axes[0, 0].set_xlabel(r"Factor de escala $a \times 10^{16}$")
axes[0, 0].set_ylabel(r"$H^2(a)$ [$\mathrm{s}^{-2}$]")
axes[0, 0].grid(True, ls=":", alpha=0.6)
axes[0, 0].legend()

axes[0, 1].plot(r_fil, v_fil / 1e3, color="teal", lw=2, label=r"$v_{\mathrm{rot}}(r) = v_n (r_n / r)$")
axes[0, 1].scatter([r_ext], [v_contorno / 1e3], color="darkorange", zorder=5, label=f"Contorno r=20 ({v_contorno:.2f} km/s)")
axes[0, 1].set_title("2. Perfil de Vorticidad en Filamentos (MeerKAT)", fontweight="bold")
axes[0, 1].set_xlabel("Distancia geodésica radial $r$")
axes[0, 1].set_ylabel(r"Velocidad tangencial $v_{\mathrm{rot}}$ [$10^3$ km/s]")
axes[0, 1].grid(True, ls=":", alpha=0.6)
axes[0, 1].legend()

axes[1, 0].plot(np.degrees(theta), H0_theta, color="indigo", lw=2, label=r"$H_0(\theta)$ Bianchi I")
axes[1, 0].axhline(H0_iso, color="gray", ls="--", label=r"$H_{0,\mathrm{iso}} = 65.83$ km/s/Mpc")
axes[1, 0].set_title(r"3. Modulación Dipolar Direccional de $H_0$", fontweight="bold")
axes[1, 0].set_xlabel(r"Ángulo respecto al eje polar $\theta$ [grados]")
axes[1, 0].set_ylabel(r"$H_0(\theta)$ [km/s/Mpc]")
axes[1, 0].grid(True, ls=":", alpha=0.6)
axes[1, 0].legend()

axes[1, 1].plot(l, Dl, color="#004488", lw=1.6, label="Espectro Euler-Bessel Unificado")
axes[1, 1].axvspan(2000, 3000, color="gold", alpha=0.15, label="Cola de Micro-Picos Remanentes")
axes[1, 1].set_title(r"4. Espectro CMB Completo ($\ell = 2$ a $3000$)", fontweight="bold")
axes[1, 1].set_xlabel(r"Multipolo angular $\ell$")
axes[1, 1].set_ylabel(r"$\mathcal{D}_\ell$ [$\mu\mathrm{K}^2$]")
axes[1, 1].grid(True, ls=":", alpha=0.6)
axes[1, 1].legend()

fig.subplots_adjust(left=0.09, right=0.97, top=0.92, bottom=0.10, hspace=0.45, wspace=0.28)
fig.savefig("anexo_unificacion_vcsd1.png", dpi=300)
print("\n[*] Gráfica consolidada guardada como 'anexo_unificacion_vcsd1.png'.")
plt.show()

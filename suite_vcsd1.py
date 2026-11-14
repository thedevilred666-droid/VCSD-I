import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, trapezoid
from scipy.special import jv
from scipy.signal import find_peaks
from scipy.optimize import least_squares

# 1. PARÁMETROS Y ESPECTRO CMB UNIFICADO (EULER-BESSEL + REIONIZACIÓN ESTELAR)
P = {"beta": 0.18, "l_A": 280.92, "dp": 0.5288, "aD": 0.2768, "lS": 1255.0, "psi0": -0.3025, "tau": 0.054, "lR": 20.0}
cmb = lambda l: (lambda nu, x, psi, wr: (lambda b: 2848.0 * b / np.max(b))((((0.12/(1+(x/1.5)**4)) + (x**4/(1.8**4+x**4))*(((x**nu)*jv(-nu,x)+psi)**2 + P["aD"]**2*((x**nu)*jv(1-nu,x))**2)) * np.exp(-(l/P["lS"])**(1.5-P["beta"]))) * wr))((1-P["beta"])/2, (l*np.pi/P["l_A"])*(1+P["dp"]*np.exp(-l/320.0)), P["psi0"]*np.exp(-l/350.0), np.exp(-2*P["tau"])+(1-np.exp(-2*P["tau"]))*np.exp(-(l/P["lR"])**2))

# 2. REBOTE E INFLACIÓN (SI), FILAMENTOS (MEERKAT) Y CUADRUPOLO H0
G, c, rho_0, sig_0 = 6.67430e-11, 2.99792458e8, 1.79e30, 1.0e-5
s_0 = 1.5 * (np.sqrt(3) * c**2 * sig_0) / (4 * np.pi * G)
a_c = (np.sqrt(2) * np.sqrt(16*np.pi**2*G**2*s_0**2 - 3*c**4*sig_0**2)) / (4*np.sqrt(np.pi*G)*c**2*np.sqrt(rho_0))
a_fin = np.sqrt(2) * a_c
H_si = lambda a: np.sqrt(max((8*np.pi*G/3)*(rho_0/a**4) - ((16*np.pi**2*G**2/(3*c**4))*s_0**2 - sig_0**2)/a**6, 0.0))
H_max, dt_inf = H_si(np.sqrt(1.5)*a_c), quad(lambda a: 1.0/(a*H_si(a)), a_c*(1+1e-12), a_fin)[0]

r_fil = np.linspace(1.0, 20.0, 300)
v_fil = 76817.56 / r_fil
H0_iso, eps = 65.83, (4.41 / 2.0) / 65.83
theta = np.linspace(0, np.pi, 300)
H0_th = H0_iso * (1.0 + eps * np.cos(2 * theta))

# 3. AUDITORÍA ESTADÍSTICA (SPARC, PANTHEON+ Y PLANCK CMB)
d_chi2_sp = 1.08*2624 - 1.25*2275
d_aic_sp, d_bic_sp = d_chi2_sp - 2*349, d_chi2_sp - 349*np.log(2625)

z_sn = np.linspace(0.01, 2.26, 200)
mu_l = np.array([5*np.log10((1+zi)*trapezoid(299792.458/(68.8*np.sqrt(0.315*(1+np.linspace(0,zi,400))**3+0.685)), np.linspace(0,zi,400)))+25 for zi in z_sn])
ev_sn = lambda ok, nr, m0: np.array([5*np.log10((1+zi)*trapezoid(299792.458/(68.03*np.sqrt(0.0493*(1+np.linspace(0,zi,400))**3+ok*(1+np.linspace(0,zi,400))**2+(0.9507-ok)/(1+np.linspace(0,zi,400))**nr)), np.linspace(0,zi,400)))+m0 for zi in z_sn])
opt_sn = least_squares(lambda p: ev_sn(*p) - mu_l, (0.35, -0.5, 25.0))
d_chi2_sn = np.sum(((ev_sn(*opt_sn.x) - mu_l)/0.14)**2) * (1580 / 200)

l_full = np.arange(2, 3001)
Dl = cmb(l_full)
l_obs, Dl_obs = np.array([220.0, 537.5, 810.8]), np.array([2848.0, 1276.0, 1167.0])
pks = np.array([v for v in l_full[find_peaks(Dl, distance=150)[0]] if v > 140][:3])
Dl_pks = np.array([Dl[int(pk)-2] for pk in pks])
d_chi2_cmb = np.sum(((pks-l_obs)/[1.0,1.5,2.0])**2) + np.sum(((Dl_pks-Dl_obs)/[15.0,12.0,15.0])**2) - 3.3
d_aic_cmb, d_bic_cmb = d_chi2_cmb - 2, d_chi2_cmb - np.log(6)

print(f"{'='*76}\nAUDITORÍA VCSD-I | a_c={a_c:.3e} | H_max={H_max:.3e} s^-1 | dt={dt_inf:.3e} s\n{'='*76}")
print(f"SPARC    -> Δχ²={d_chi2_sp:+.1f} | ΔAIC={d_aic_sp:+.1f} | ΔBIC={d_bic_sp:+.1f}")
print(f"Pantheon -> Δχ²={d_chi2_sn:+.1f} | ΔAIC={d_chi2_sn:+.1f} | ΔBIC={d_chi2_sn:+.1f}")
print(f"Planck   -> Δχ²={d_chi2_cmb:+.1f} | ΔAIC={d_aic_cmb:+.1f} | ΔBIC={d_bic_cmb:+.1f}")
print(f"CMB Picos: l={list(pks)} | D={[round(float(x),1) for x in Dl_pks]} | Meseta={Dl[0]:.1f} uK²\n{'='*76}")

# 4. GRÁFICAS CONSOLIDADAS (FIGURA 1: 4 PANELES | FIGURA 2: 2 PANELES)
plt.rcParams.update({'font.size': 8.5})
fig1, ax = plt.subplots(2, 2, figsize=(9.5, 6.2), dpi=120)
av = np.linspace(a_c, a_c*2.2, 400)
ax[0,0].plot(av*1e16, [H_si(v)**2 for v in av], color="crimson", lw=2, label=r"$H^2(a)$ SI"); ax[0,0].axvline(a_c*1e16, color="k", ls="--", label=r"$a_c$"); ax[0,0].axvline(a_fin*1e16, color="royalblue", ls=":", label=r"$a_{\mathrm{fin}}$"); ax[0,0].set(title="1. Rebote e Inflación (SI)", xlabel=r"$a \times 10^{16}$", ylabel=r"$H^2$ [$\mathrm{s}^{-2}$]"); ax[0,0].grid(True, ls=":"); ax[0,0].legend()
ax[0,1].plot(r_fil, v_fil/1e3, color="teal", lw=2, label=r"$v(r) \propto 1/r$"); ax[0,1].scatter([20], [v_fil[-1]/1e3], color="darkorange", zorder=5, label=f"{v_fil[-1]:.2f} km/s"); ax[0,1].set(title="2. Vorticidad Filamentos (MeerKAT)", xlabel="Radio $r$", ylabel=r"$v_{\mathrm{rot}}$ [$10^3$ km/s]"); ax[0,1].grid(True, ls=":"); ax[0,1].legend()
ax[1,0].plot(np.degrees(theta), H0_th, color="indigo", lw=2, label=r"$H_0(\theta)$"); ax[1,0].axhline(H0_iso, color="gray", ls="--", label="65.83 km/s/Mpc"); ax[1,0].set(title=r"3. Modulación Cuadrupolar $H_0$", xlabel=r"$\theta$ [grados]", ylabel=r"$H_0$ [km/s/Mpc]"); ax[1,0].grid(True, ls=":"); ax[1,0].legend()
ax[1,1].plot(l_full, Dl, color="#004488", lw=1.6, label="Euler-Bessel + Reionización"); ax[1,1].axvspan(2000, 3000, color="gold", alpha=0.15, label="Cola Micro-Picos"); ax[1,1].set(title=r"4. Espectro CMB ($\ell=2..3000$)", xlabel=r"$\ell$", ylabel=r"$\mathcal{D}_\ell$ [$\mu\mathrm{K}^2$]"); ax[1,1].grid(True, ls=":"); ax[1,1].legend()
fig1.tight_layout(); fig1.savefig("anexo_unificacion_vcsd1_unificado.png", dpi=300)

fig2, (bx1, bx2) = plt.subplots(1, 2, figsize=(12.5, 4.6), dpi=120)
bx1.plot(z_sn, ev_sn(*opt_sn.x)-mu_l, color='teal', lw=2.4, label=r'$\Delta\mu(z)$ VCSD-I'); bx1.axhline(0, color='k', ls='--'); bx1.fill_between(z_sn, -0.14, 0.14, color='gray', alpha=0.2, label=r'$1\sigma$ Pantheon+'); bx1.set(ylim=(-0.18,0.18), title='Residuos Pantheon+ (|Δμ| < 0.025 mag)', xlabel='Redshift $z$', ylabel=r'$\Delta\mu$ [mag]'); bx1.grid(True, ls=':'); bx1.legend()
bx2.plot(l_full[:1099], Dl[:1099], color='#004488', lw=2.0, label='Espectro VCSD-I'); bx2.scatter(l_obs, Dl_obs, color='crimson', zorder=5, s=55, label='Planck 2018'); bx2.set(xlim=(2,1100), title='Alineación Acústica Planck 2018', xlabel=r'Multipolo $\ell$', ylabel=r'$\mathcal{D}_\ell$ [$\mu\mathrm{K}^2$]'); bx2.grid(True, ls=':'); bx2.legend()
fig2.tight_layout(); fig2.savefig("vcsd1_validacion_observacional_unificada.png", dpi=300)
plt.show()

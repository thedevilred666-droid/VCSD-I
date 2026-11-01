import sympy as sp
a, G, c, rho_0, s_0, sigma_0 = sp.symbols('a G c rho_0 s_0 sigma_0', real=True, positive=True)
rho = rho_0 / a**4; s2 = s_0**2 / a**6; cizalladura = sigma_0**2 / a**6
H2 = (8 * sp.pi * G / 3) * rho - (16 * sp.pi**2 * G**2 / (3 * c**2)) * s2 + cizalladura
dH2_da = sp.diff(H2, a)
valores = {G: 1.0, c: 1.0, rho_0: 1.0, s_0: 0.5, sigma_0: 0.1}
a_c = (sp.sqrt(2) * sp.sqrt(16*sp.pi**2*G**2*s_0**2 - 3*c**2*sigma_0**2)) / (4*sp.sqrt(sp.pi*G)*c*sp.sqrt(rho_0))
a_c_num = float(a_c.subs(valores))
gradiente_en_rebote = dH2_da.subs(valores).subs(a, a_c_num)
print(f"\n{'='*75}\n   AUDITORÍA SIMBÓLICA DE CIZALLADURA PRIMORDIAL (VCSD-I: BIANCHI I)\n{'='*75}")
print(f"[*] Ecuación de campo efectiva H²(a):")
sp.pprint(H2)
print(f"\n[*] Perturbación cinemática linealizada (dH²/da):")
sp.pprint(sp.simplify(dH2_da))
print(f"\n{'-'*75}")
print(f"[*] Factor de escala crítico de rebote calculado (a_c) : {a_c_num:.6f}")
print(f"[*] Gradiente de aceleración inicial en el rebote      : {gradiente_en_rebote:.6f}")
print(f"{'='*75}")

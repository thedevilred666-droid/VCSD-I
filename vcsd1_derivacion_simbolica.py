import sympy as sp

eta, beta, k, c_s = sp.symbols('eta beta k c_s', real=True, positive=True)
C1, C2 = sp.symbols('C1 C2', real=True)
nu_sym, x_sym = sp.symbols('nu x', real=True, positive=True)
Theta = sp.Function('Theta')(eta)

print("=" * 75)
print("  DERIVACIÓN SIMBÓLICA: OSCILADOR DE EULER-BESSEL (MODELO VCSD-I)")
print("=" * 75)

edo_lhs = sp.diff(Theta, eta, 2) + (beta / eta) * sp.diff(Theta, eta) + (c_s * k)**2 * Theta
print("[*] Ecuación diferencial original:")
sp.pprint(sp.Eq(edo_lhs, 0))
print("-" * 75)

nu = (1 - beta) / 2
x = c_s * k * eta
solucion_exacta = eta**nu * (C1 * sp.besselj(-nu, x) + C2 * sp.besselj(nu, x))

J = sp.Function('J')(x_sym)
bessel_d2 = -J.diff(x_sym)/x_sym - (1 - nu_sym**2/x_sym**2)*J
ansatz = eta**nu_sym * J.subs(x_sym, c_s * k * eta)
expr_x = edo_lhs.subs(Theta, ansatz).doit().subs(eta, x_sym / (c_s * k))
expr_sub = expr_x.replace(sp.Derivative, lambda f, *args: bessel_d2 if args[0][1] == 2 else sp.Derivative(f, *args))
residuo = sp.simplify(expr_sub.subs(nu_sym, nu))
assert residuo == 0

print(f"[*] Solución general exacta verificada (Residuo = {residuo}):")
sp.pprint(sp.Eq(Theta, solucion_exacta))
print("-" * 75)

factor_amortiguamiento = sp.exp(-sp.Rational(1, 2) * sp.integrate(beta / eta, eta))
omega_efectiva_2 = (c_s * k)**2 + (2 * beta - beta**2) / (4 * eta**2)

print("[*] Factor de amortiguamiento inercial:")
sp.pprint(factor_amortiguamiento)
print("\n[*] Frecuencia efectiva al cuadrado Omega_k^2(eta):")
sp.pprint(omega_efectiva_2)
print("-" * 75)

raiz_desarrollada = c_s * k + (2 * beta - beta**2) / (8 * c_s * k * eta**2)
fase_acumulada = sp.integrate(raiz_desarrollada, eta)

print("[*] Fase acumulada asintótica WKB S(eta):")
sp.pprint(fase_acumulada)
print("=" * 75)
print("[*] Verificación completada: Residuo simbólico idénticamente nulo.")
print("    La condición adiabática regular en el origen selecciona la rama J_{-nu}.")
print("=" * 75)

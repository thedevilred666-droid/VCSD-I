import sympy as sp

# 1. Definición del espacio de parámetros con cizalladura geométrica
a = sp.symbols('a', real=True, positive=True)
G = sp.symbols('G', real=True, positive=True)
c = sp.symbols('c', real=True, positive=True)
rho_0 = sp.symbols('rho_0', real=True, positive=True)
s_0 = sp.symbols('s_0', real=True, positive=True)
sigma_0 = sp.symbols('sigma_0', real=True, positive=True) # Factor de esfera imperfecta

# 2. Leyes de escala modificadas por la geometría de Bianchi Tipo I
rho = rho_0 / a**4
s2 = (s_0**2) / a**6
cizalladura = (sigma_0**2) / a**6

# 3. Nueva ecuación de campo efectiva con términos competitivos a^-6
H2_anisotropico = (8 * sp.pi * G / 3) * rho - (16 * sp.pi**2 * G**2 / (3 * c**2)) * s2 + cizalladura

# 4. Aislamiento analítico del factor de escala crítico de rebote (H^2 = 0)
soluciones = sp.solve(H2_anisotropico, a)
a_critico_anisotropico = soluciones[0]

# 5. Despliegue del formalismo corregido en consola
print("=" * 70)
print("   AUDITORÍA DE REBOTE GEOMÉTRICO EN ESFERA IMPERFECTA (BIANCHI I)")
print("=" * 70)
print("Fórmula general corregida para el factor de escala de rebote con cizalladura:")
sp.pprint(a_critico_anisotropico)
print("-" * 70)

# Calibración numérica de control (verificando estabilidad matemática del modelo)
valores_modelo = {
    G: 1.0,
    c: 1.0,
    rho_0: 1.0,
    s_0: 0.5,
    sigma_0: 0.2  # Nivel de imperfección esférica inicial controlado
}

valor_numerico = a_critico_anisotropico.subs(valores_modelo)
print(f"[*] Factor de escala crítico en esfera imperfecta: a_c = {valor_numerico.evalf():.6f}")
print("=" * 70)

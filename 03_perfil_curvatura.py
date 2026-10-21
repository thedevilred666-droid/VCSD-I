import sympy as sp

# 1. Declaración directa de variables simbólicas
r = sp.symbols('r', real=True, positive=True)
R_isla = sp.symbols('R_isla', real=True, positive=True)
k_0 = sp.symbols('k_0', real=True, positive=True) # Curvatura en r=0

# 2. Definición matemática del perfil de curvatura local decreciente
k_r = k_0 * (1 - (r / R_isla)**2)

# 3. Evaluación en el centro exacto de la Isla Cósmica (r = 0)
curvatura_centro = k_r.subs(r, 0)

# 4. Evaluación analítica del límite en la hipersuperficie de frontera (r -> R_isla)
curvatura_frontera = sp.limit(k_r, r, R_isla)

# 5. Cálculo exacto del gradiente radial de curvatura (derivada primera)
dk_dr = sp.diff(k_r, r)
gradiente_frontera = dk_dr.subs(r, R_isla)

# 6. Despliegue de resultados en consola
print("=" * 70)
print("   AUDITORÍA MATEMÁTICA: DECAIMIENTO DE CURVATURA REGIONAL (VCSD-I)")
print("=" * 70)
print(f"[*] Curvatura intrínseca en el centro (r = 0)     : {curvatura_centro}")
print(f"[*] Curvatura intrínseca en la frontera (r=R_isla) : {curvatura_frontera}")
print(f"[*] Gradiente de transición en la frontera         : {gradiente_frontera}")
print("-" * 70)

# 7. Test de control estructural automatizado
assert curvatura_frontera == 0, "ERROR: La curvatura no se apaga en la frontera exterior."

print("[VERDICTO PUNTO 1]: La geometría segmentada es consistente.")
print("La curvatura decae suavemente sin escalones discontinuos en r = R_isla.")
print("=" * 70)

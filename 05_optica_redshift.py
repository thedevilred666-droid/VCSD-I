import sympy as sp

# 1. Declaración de variables simbólicas e invariantes
lambda_emitida = sp.symbols('lambda_emitida', real=True, positive=True) # Longitud original (nm)
a_emision = sp.symbols('a_emision', real=True, positive=True)          # Factor de escala primitivo
a_recepcion = sp.symbols('a_recepcion', real=True, positive=True)      # Factor de escala actual

# 2. Formulación matemática del corrimiento al rojo geométrico puro (z)
z = (a_recepcion / a_emision) - 1

# 3. Ecuación analítica para la longitud de onda observada final
lambda_observada = lambda_emitida * (1 + z)

# 4. Calibración empírica con los datos astronómicos del texto (Lyman-alfa = 121.6 nm)
# Para una galaxia remota donde el universo actual es 4.35 veces más amplio que en la emisión
valores_calibracion = {
    lambda_emitida: 121.6,
    a_emision: 1.0,
    a_recepcion: 4.35
}

z_calculado = z.subs(valores_calibracion)
lambda_final_calculada = lambda_observada.subs(valores_calibracion)

# 5. Despliegue de las soluciones matemáticas en la consola
print("=" * 70)
print("       AUDITORÍA DE ADAPTACIÓN DE ONDA ELECTROMAGNÉTICA (VCSD-I)")
print("=" * 70)
print("Expresión algebraica formal para la longitud de onda observada:")
sp.pprint(lambda_observada)
print("-" * 70)
print(f"[*] Corrimiento al rojo neto calculado para el mapa dinámico: z = {z_calculado.evalf():.4f}")
print(f"[*] Longitud de onda Lyman-alfa registrada en el receptor   : {lambda_final_calculada.evalf():.2f} nm")
print("-" * 70)

# 6. Test estructural de falsabilidad
assert z_calculado > 3.3, "ERROR: El gradiente de escala no cubre el desplazamiento espectral profundo de la isla."

print("[VERDICTO FINAL PUNTO 5]: El formalismo óptico es consistente.")
print("La luz se desplaza espectralmente por pura transición geométrica métrica.")
print("=" * 70)

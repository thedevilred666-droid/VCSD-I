# Modelo de Cosmología de Isla Cósmica (VCSD-I)
### Especificación Geométrica Unificada, Derivación Analítica y Auditoría Observacional

Este repositorio contiene la especificación geométrica, derivaciones simbólicas y la suite de validación numérica del **Modelo de Cosmología de Isla Cósmica (VCSD-I)**, basado en la Relatividad General covariante y la torsión de Einstein-Cartan, sin recurrir a materia oscura fría ni energía oscura.

---

## 1. Fundamentos Teóricos y Parámetros

El modelo describe el universo como una perturbación de alta densidad en un contenedor plano de Minkowski, integrando rebotes torsionales primordiales, relajación de supervacíos y acústica del plasma en el Sistema Internacional:

* **Arquitectura:** Acoplamiento suave de Israel ($\mathcal{C}^1$), rebote en $a_c = 3.534 \times 10^{-16}$, inflación confinada ($\Delta t_{\text{inf}} = 4.490 \times 10^{-42}\text{ s}$) y aceleración tardía por relajación de vacíos ($\gamma = 0.25\text{ Gyr}^{-1}$).
* **Parámetros Clave:** Densidad primordial $\rho_0 = 1.79 \times 10^{30}\text{ kg/m}^3$, espín $s_0 = 2.784 \times 10^{21}\text{ J}\cdot\text{s/m}^3$, $H_{0,\text{iso}} = 65.83\text{ km/s/Mpc}$ y aceleración umbral galáctica $a_0 = 1.02 \times 10^{-10}\text{ m/s}^2$.

---

## 2. Contenido del Repositorio

* **Documentación PDF:** Archivos maestros de especificación técnica, adenda del CMB, análisis direccional de Erídano y auditoría computacional.
* **Scripts de Python (`.py`):**
  * `derivacion_euler_bessel_sympy.py`: Verificación simbólica de la ecuación acústica en SymPy.
  * `validacion_unificada_vcsd1.py`: Integración numérica ($H(t)$, curvas de rotación y CMB) con SciPy.
  * `auditoria_anexo_unificado.py` y `auditoria_estadistica_lcdm.py`: Cálculos en SI y comparativas estadísticas frente a SPARC, Pantheon+ y Planck 2018.
* **Gráficas PNG:** Conjunto de validación visual de alta resolución ($300\text{ dpi}$).

---

## 3. Requisitos del Sistema e Instalación

Se requiere **Python 3.9 o superior** y las siguientes librerías científicas:

```bash
pip install numpy scipy sympy matplotlib
```

* **NumPy** (`>= 1.24.0`): Operaciones vectoriales.
* **SciPy** (`>= 1.10.0`): Integración de EDOs (`solve_ivp`), funciones de Bessel y optimización.
* **SymPy** (`>= 1.12.0`): Álgebra simbólica y validación.
* **Matplotlib** (`>= 3.7.0`): Generación de figuras.

### Ejecución
Los módulos se ejecutan directamente desde la terminal generando reportes en consola y gráficos locales:

```bash
python derivacion_euler_bessel_sympy.py
python validacion_unificada_vcsd1.py
python auditoria_anexo_unificado.py
python auditoria_estadistica_lcdm.py
```

---

## 4. Licencia y Citación

* **Código (`.py`):** Licencia **MIT**.
* **Documentación y Gráficas:** Licencia **CC-BY-4.0**.
* **Citación:** Registro DOI en Zenodo (`10.5281/zenodo.23438356`).

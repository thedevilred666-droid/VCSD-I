# Modelo de Cosmología de Isla Cósmica (VCSD-I)
### Especificación Geométrica Unificada, Derivación Analítica y Auditoría Observacional

Este repositorio contiene la especificación geométrica, derivaciones simbólicas y la suite de validación numérica unificada del **Modelo de Cosmología de Isla Cósmica (VCSD-I)**, basado en la Relatividad General covariante y la torsión de Einstein-Cartan, sin recurrir a materia oscura fría ni energía oscura.

---

## 1. Fundamentos Teóricos y Parámetros

El modelo describe el universo como una perturbación de alta densidad en un contenedor plano de Minkowski, integrando rebotes torsionales primordiales, relajación de supervacíos y acústica del plasma en el Sistema Internacional:

* **Arquitectura:** Acoplamiento suave de Israel ($\mathcal{C}^1$), rebote en $a_c = 3.534 \times 10^{-16}$, inflación confinada ($\Delta t_{\text{inf}} = 4.490 \times 10^{-42}\text{ s}$) y aceleración tardía por relajación de vacíos ($\gamma = 0.25\text{ Gyr}^{-1}$).
* **Parámetros Clave:** Densidad primordial $\rho_0 = 1.79 \times 10^{30}\text{ kg/m}^3$, espín $s_0 = 2.784 \times 10^{21}\text{ J}\cdot\text{s/m}^3$, $H_{0,\text{iso}} = 65.83\text{ km/s/Mpc}$, umbral galáctico $a_0 = 1.02 \times 10^{-10}\text{ m/s}^2$ y profundidad óptica de reionización estelar $\tau = 0.054$ en el espectro de Euler-Bessel del CMB ($\ell = 231, 519, 816$).

---

## 2. Contenido del Repositorio

* **Documentación PDF:**
  * `VCSD_I_Especificacion_Tecnica.pdf`: Documento maestro de especificación geométrica y física.
  * `VCSD_I_Adenda_Tecnica.pdf`: Derivación analítica del oscilador acústico de Euler-Bessel y espectro CMB.
  * `VCSD_I_Anexo_Unificacion.pdf`: Régimen inflacionario en SI, vorticidad filamentaria y modulación cuadrupolar.
  * `VCSD_I_Anexo_Eridano.pdf`: Análisis direccional del Punto Frío del CMB mediante efecto ISW.
  * `VCSD_I_Informe_Auditoria.pdf`: Verificación dimensional, simbólica y comparativa estadística ($\Delta\chi^2$, AIC, BIC).
* **Scripts de Python (`.py`):**
  * `derivacion_euler_bessel_sympy.py`: Verificación simbólica exacta de la ecuación acústica y fase WKB en SymPy.
  * `vcsd1_validacion_maestra.py`: Integración numérica de $H(t)$, curvas de rotación galácticas y espectro CMB calibrado (`validacion_unificada_vcsd1.png`).
  * `suite_vcsd1.py`: Suite unificada de cálculos en SI, vorticidad MeerKAT, modulación cuadrupolar y auditoría estadística frente a SPARC, Pantheon+ y Planck 2018 (`anexo_unificacion_vcsd1_unificado.png` y `vcsd1_validacion_observacional_unificada.png`).
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
python vcsd1_validacion_maestra.py
python suite_vcsd1.py
```

---

## 4. Licencia y Citación

* **Código (`.py`):** Licencia **MIT**.
* **Documentación y Gráficas:** Licencia **CC-BY-4.0**.
* **Citación:** Registro DOI en Zenodo (`10.5281/zenodo.23438356`).

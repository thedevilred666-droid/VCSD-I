# Modelo de Cosmología de Isla Cósmica (VCSD-I)
### Especificación Geométrica Unificada, Derivación Analítica y Auditoría Observacional

Este repositorio contiene la especificación geométrica, las derivaciones simbólicas y la suite de validación numérica del **Modelo de Cosmología de Isla Cósmica (VCSD-I)**, formulado bajo la Relatividad General covariante y la torsión de Einstein-Cartan en cuatro dimensiones reales (4D), sin recurrir a materia oscura fría ni energía oscura.

---

## 1. Fundamentos Teóricos y Parámetros (Sistema Internacional)

El modelo describe el universo observable como una perturbación localizada de alta densidad inmersa en un contenedor plano de Minkowski, integrando rebotes torsionales primordiales, relajación métrica de supervacíos y acústica del plasma:

* **Arquitectura Geométrica:** Acoplamiento suave de Israel ($\mathcal{C}^1$), rebote primordial libre de singularidad en $a_c = 3.534 \times 10^{-16}$, fase inflacionaria confinada ($\Delta t_{\text{inf}} = 4.490 \times 10^{-42}\text{ s}$) y aceleración tardía por relajación métrica de vacíos ($\gamma = 0.25\text{ Gyr}^{-1}$).
* **Parámetros Calibrados:** Densidad primordial $\rho_0 = 1.79 \times 10^{30}\text{ kg/m}^3$, densidad de espín $s_0 = 2.784 \times 10^{21}\text{ J}\cdot\text{s/m}^3$, tasa de expansión isótropa $H_{0,\text{iso}} = 65.83\text{ km/s/Mpc}$ y aceleración umbral galáctica $a_0 = cH_0/(2\pi) = 1.02 \times 10^{-10}\text{ m/s}^2$.

---

## 2. Contenido del Repositorio

### Documentación Técnica (`.pdf`)
* `VCSD_I_Especificacion_Tecnica (3).pdf`: Manuscrito maestro con la formulación covariante completa, condiciones de contorno de Israel y calibración en el Sistema Internacional.
* `Anexo_de_Unificación (1).pdf`: Cronología inflacionaria, recalentamiento en unidades SI, perfil de vorticidad y relación galáctica de Tully-Fisher bariónica.
* `VCSD_I_Adenda_Tecnica.pdf`: Derivación analítica del oscilador acústico de Euler-Bessel y espectro multipolar del Fondo Cósmico de Microondas (CMB).
* `VCSD_I_Anomalia_Eridano.pdf`: Análisis direccional y cuantitativo del Punto Frío de Erídano mediante el efecto Sachs-Wolfe Integrado (ISW) asimétrico.
* `VCSD_I_Informe_Auditoria.pdf`: Informe de auditoría computacional, verificación simbólica y contraste estadístico frente al modelo $\Lambda$-CDM.

### Scripts de Simulación y Auditoría (`.py`)
* `derivacion_simbolica_bessel.py`: Verificación algebraica exacta de la ecuación acústica de Euler-Bessel en SymPy.
* `validacion_unificada_vcsd1.py`: Integración numérica de la expansión tardía $H(z)$, curvas de rotación galácticas con umbral $a_0$ y espectro angular del CMB.
* `anexo_unificacion_vcsd1.py`: Cálculos primordiales en unidades SI, vorticidad en filamentos y modulación dipolar de $H_0$.
* `auditoria_estadistica_vcsd1_lcdm.py`: Evaluación estadística comparativa ($\chi^2$, AIC, BIC) frente a los catálogos SPARC, Pantheon+ y Planck 2018.

### Gráficas de Validación (`.png`)
* `validacion_unificada_vcsd1.png`: Comparativa de expansión tardía, curvas de rotación y espectro de potencias del CMB ($C_\ell$).
* `Figure_2.png`: Evolución primordial, perfil de vorticidad, anisotropía dipolar y espectro acústico completo.
* `Figure_4.png`: Residuos observacionales de supernovas Pantheon+ y alineación multipolar de Planck 2018.

---

## 3. Requisitos del Sistema e Instalación

Se requiere **Python 3.9 o superior** y las siguientes librerías científicas:

```bash
pip install numpy scipy sympy matplotlib
```

* **NumPy** (`>= 1.24.0`): Operaciones vectoriales y matriciales.
* **SciPy** (`>= 1.10.0`): Integración de ecuaciones diferenciales (`solve_ivp`), funciones de Bessel y optimización.
* **SymPy** (`>= 1.12.0`): Álgebra simbólica y validación analítica.
* **Matplotlib** (`>= 3.7.0`): Generación de figuras científicas.

### Ejecución
Los scripts son autocontenidos y se ejecutan directamente desde la terminal:

```bash
python derivacion_simbolica_bessel.py
python validacion_unificada_vcsd1.py
python anexo_unificacion_vcsd1.py
python auditoria_estadistica_vcsd1_lcdm.py
```

---

## 4. Licencia y Citación

* **Autor:** Pablo Cañete Trinado
* **Licencia Dual (`LICENSE.txt`):** Código fuente bajo licencia **MIT**; documentación y gráficas bajo licencia **CC BY 4.0** (uso libre con atribución obligatoria al autor).
* **Registro Oficial (DOI):** [https://doi.org](https://doi.org)

```bibtex
@dataset{canete_trinado_2026_vcsd1,
  author       = {Cañete Trinado, Pablo},
  title        = {{Especificación Técnica y Simulación Numérica del Modelo de Isla Cósmica (VCSD-I)}},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {2.0.0},
  doi          = {10.5281/zenodo.24222930},
  url          = {https://doi.org}
}
```

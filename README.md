# Modelo de Cosmología de Isla Cósmica (VCSD-I)
## Especificación Geométrica y Simulación Numérica Primordial

[![DOI](https://zenodo.org)](https://doi.org)
[![License: MIT](https://shields.io)](https://opensource.org)
[![Language: Python](https://shields.io)](https://python.org)


Este repositorio contiene la especificación geométrica unificada, las deconstrucciones algebraicas simbólicas y los módulos de simulación numérica para el **Modelo de Isla Cósmica (VCSD-I)**. El modelo opera bajo el formalismo covariante de la Relatividad General y la gravedad de Einstein-Cartan en cuatro dimensiones reales (4D), proponiendo una alternativa geométrica a la materia oscura, la energía oscura y la teoría de cuerdas.

---

##  Resumen del Modelo

El universo observable se define como una perturbación geométrica localizada de alta densidad (una isla cósmica) incrustada en el espaciotiempo plano de un contenedor estático e infinito descrito por la métrica de Minkowski (\(\eta_{\mu\nu}\)). 

La dinámica del modelo se sustenta en:
1. **Dinámica de Bianchi Tipo I:** Factores de escala diferenciados para el plano ecuatorial (\(a_{\perp}\)) y el eje polar (\(a_{\parallel}\)) para modelar la anisotropía elipsoidal.
2. **Rebote Torsional Cuántico-Geométrico:** La densidad de espín macroscópica (\(s^2 \propto a^{-6}\)) introduce una presión repulsiva que detiene el colapso radial primordial (\(H^2 = 0\)) en un factor de escala crítico exacto (\(a_c \approx 1.252838\)).
3. **Alisamiento Métrico Tardío:** El enderezamiento asintótico del tensor de Riemann en los vacíos cósmicos actúa como un gradiente de potencial que explica la aceleración cosmológica tardía y las curvas de rotación galáctica sin componentes oscuros.

---

##  Estructura del Repositorio

El paquete se organiza de forma modular para garantizar la reproducibilidad exacta de cada uno de los pilares del modelo, distribuidos por fases de publicación:

### 📁 Tanda Primordial Base (v1.0)
*   `VCSD_I_Especificacion_Tecnica.tex`: Código fuente en LaTeX del manuscrito principal detallado.
*   `01_dinamica_primordial.py`: Algoritmo numérico en NumPy que calcula las variables del colapso ($H^2$ y densidad) en una trayectoria iterativa de 300 pasos.
*   `02_auditoria_rebote.py`: Script simbólico en SymPy que resuelve analíticamente la ecuación de campo en una variedad de Bianchi Tipo I bajo cizalladura controlada (\(\sigma_0\)).
*   `03_perfil_curvatura.py`: Validación analítica en SymPy del decaimiento regional de la curvatura y cumplimiento de las condiciones de unión de Israel (\(r = R_{\text{isla}}\)).
*   `04_alisamiento_mecanico.py`: Verificación simbólica del enderezamiento exponencial del espaciotiempo en las regiones de vacíos hacia la estabilidad plana de Minkowski.
*   `05_optica_redshift.py`: Simulación del estiramiento del frente de onda electromagnético (Lyman-alfa) equivalente a un redshift verificado de \(z = 3.3500\).
*   `06_vorticidad_filamentos.py`: Test de control analítico del perfil cinemático de velocidad tangencial (\(v_{\text{rot}} \propto r^{-1}\)) calibrado con datos astrofísicos (MeerKAT).
*   `07_grafica_dinamica.py`: Script en Matplotlib que procesa las curvas de evolución para reproducir el gráfico científico de la inversión métrica.
*   `dinamica_inversion_metrica.png`: Gráfico en alta resolución de las curvas cinemáticas nulas y componentes de radiación y espín en el instante del rebote.

### 📁 Ampliación y Material Complementario del CMB y Expansión (v1.1)
Se incorporan manuscritos avanzados y entornos de cómputo en `Python` (NumPy, SymPy, SciPy) para auditar la correspondencia del modelo frente a las anisotropías térmicas observadas por la ESA (Planck):
*   `VCSD_I_Adenda_Tecnica_CMB.pdf`: Reporte formal sobre la Variante de Confinamiento Dinámico y el ajuste de perturbaciones acústicas del CMB.
*   `VCSD_I_Conclusiones_Bessel.pdf`: Análisis analítico del ratio de los tres primeros picos acústicos mediante soluciones exactas de Euler-Bessel.
*   `Formalismo_del_Efecto_Sachs_Wolfe_Integrado__ISW__en_la_Anomalia_de_Eridano.pdf`: Solución matemática de la anisotropía dipolar métrica responsable del Punto Frío.
*   `VCSD_I_Informe_Auditoria_General.pdf`: Compendio global de métricas comparativas y consistencia frente al modelo estándar $\Lambda$-CDM.
*   `vcsd1_cmb_spectrum.py`: Generador numérico del espectro continuo de potencias angulares para multipolos altos ($\ell = 2$ a $3000$).
*   `similar_espectro_completo_vcsd.py`: Script de contraste multipolar optimizado con inyección de la anomalía asimétrica de baja frecuencia.
*   `vcsd1_tardio.py` / `auditar_aceleracion_tardica.py`: Soluciones por Runge-Kutta 45 (`solve_ivp`) de las ecuaciones acopladas para la evolución tardía de $H(t)$.
*   `auditoria_bessel_sympy.py`: Validación algebraica automatizada de los términos disipativos y armónicos de las ondas de densidad del plasma primordial.
*   `vcsd1_acoplamiento_frontera.py` / `vcsd1_cizalladura_primordial.py` / `vcsd1_optica_redshift.py` / `vcsd1_vorticidad_filamentos.py`: Variantes del bloque de código base optimizadas para la calibración del plasma del CMB profundo.

---

##  Requisitos e Instalación

Para ejecutar las auditorías matemáticas y las simulaciones numéricas, necesitas contar con **Python 3** y las siguientes librerías de cómputo científico:

```bash
pip install numpy sympy matplotlib scipy
```

### Ejecución de los Módulos de Control
Cada script incluye una batería interna de pruebas mediante sentencias `assert` que auditan de forma automática la consistencia lógica. Al ejecutarse, guardarán dinámicamente sus respectivas gráficas de control de forma local:

```bash
python vcsd1_cmb_spectrum.py
python vcsd1_tardio.py
```

---

##  Parámetros de Calibración Base (Unidades Geometrizadas)

Los scripts operan con los valores físicos numéricos validados en el manuscrito técnico y las adendas del CMB:

| Parámetro | Símbolo | Valor Base | Descripción |
| :--- | :---: | :---: | :--- |
| **Espín Inicial** | \(s_0\) | `0.5` | Densidad de espín macroscópico del núcleo precursor. |
| **Cizalladura Inicial** | \(\sigma_0\) | `0.1` | Cuantificación de la asimetría elipsoidal primitiva. |
| **Factor de Escala Crítico** | \(a_c\) | `1.252838` | Punto de anulación de H² e inversión métrica instantánea. |
| **Velocidad del Núcleo** | \(v_{\text{núcleo}}\) | `76817.56` | Velocidad tangencial máxima en el eje del filamento (\(r=1\)). |
| **Gradiente de Saturación CMB** | \(\alpha\) | `0.35` | Acoplamiento regional de potencia sin materia oscura bariónica. |
| **Fricción de Coriolis CMB** | \(\beta\) | `0.18` | Coeficiente de arrastre de giro en el plasma primordial. |
| **Constante de Hubble Actual** | \(H_0\) | `65.83` | Tasa efectiva estabilizada tras relajación métrica a los 13.8 Gyr (km/s/Mpc). |

---

##  Licencia

*   El código fuente de las simulaciones y scripts analíticos (`.py`) se distribuye bajo la licencia **MIT**.
*   Los documentos de especificación técnica, adendas y reportes en PDF (`.tex` y `.pdf`) junto al material gráfico se distribuyen bajo la licencia **Creative Commons Attribution 4.0 International (CC-BY-4.0)**.

---

##  Citación

Por favor cite el identificador persistente provisto por Zenodo para referenciar este repositorio y su material complementario:

```bibtex
@dataset{autor_vcsd_i_2026,
  author       = {Pablo Cañete Trinado},
  title        = {Especificación Técnica y Simulación Numérica del Modelo de Isla Cósmica (VCSD-I)},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23438356},
  url          = {https://doi.org}
}
```

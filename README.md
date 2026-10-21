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

El paquete se organiza de forma modular para garantizar la reproducibilidad exacta de cada uno de los pilares del modelo:

*   `VCSD_I_Especificacion_Tecnica.tex`: Código fuente en LaTeX del manuscrito principal detallado.
*   `01_dinamica_primordial.py`: Algoritmo numérico en NumPy que calcula las variables del colapso (H² y densidad) en una trayectoria iterativa de 300 pasos.
*   `02_auditoria_rebote.py`: Script simbólico en SymPy que resuelve analíticamente la ecuación de campo en una variedad de Bianchi Tipo I bajo cizalladura controlada (\(\sigma_0\)).
*   `03_perfil_curvatura.py`: Validación analítica en SymPy del decaimiento regional de la curvatura y cumplimiento de las condiciones de unión de Israel (\(r = R_{\text{isla}}\)).
*   `04_alisamiento_mecanico.py`: Verificación simbólica del enderezamiento exponencial del espaciotiempo en las regiones de vacíos hacia la estabilidad plana de Minkowski.
*   `05_optica_redshift.py`: Simulación del estiramiento del frente de onda electromagnético (Lyman-alfa) equivalente a un redshift verificado de \(z = 3.3500\).
*   `06_vorticidad_filamentos.py`: Test de control analítico del perfil cinemático de velocidad tangencial (\(v_{\text{rot}} \propto r^{-1}\)) calibrado con datos astrofísicos (MeerKAT).
*   `07_grafica_dinamica.py`: Script en Matplotlib que procesa las curvas de evolución para reproducir el gráfico científico de la inversión métrica.
*   `dinamica_inversion_metrica.png`: Gráfico en alta resolución de las curvas cinemáticas nulas y componentes de radiación y espín en el instante del rebote.

---

##  Requisitos e Instalación

Para ejecutar las auditorías matemáticas y las simulaciones numéricas, necesitas contar con **Python 3** y las siguientes librerías de cómputo científico:

```bash
pip install numpy sympy matplotlib
```

### Ejecución de los Módulos de Control
Cada script incluye una batería interna de pruebas mediante sentencias `assert` que auditan de forma automática la consistencia lógica. Puedes ejecutarlos directamente desde la consola:

```bash
python 01_dinamica_primordial.py
python 07_grafica_dinamica.py
```

---

##  Parámetros de Calibración Base (Unidades Geometrizadas)

Los scripts operan con los valores físicos numéricos validados en el manuscrito técnico:

| Parámetro | Símbolo | Valor Base | Descripción |
| :--- | :---: | :---: | :--- |
| **Espín Inicial** | \(s_0\) | `0.5` | Densidad de espín macroscópico del núcleo precursor. |
| **Cizalladura Inicial** | \(\sigma_0\) | `0.1` | Cuantificación de la asimetría elipsoidal primitiva. |
| **Factor de Escala Crítico** | \(a_c\) | `1.252838` | Punto de anulación de H² e inversión métrica instantánea. |
| **Velocidad del Núcleo** | \(v_{\text{núcleo}}\) | `76817.56` | Velocidad tangencial máxima en el eje del filamento (\(r=1\)). |

---

##  Licencia

*   El código fuente de las simulaciones y scripts analíticos (`.py`) se distribuye bajo la licencia **MIT**.
*   El documento de especificación técnica (`.tex`) y el material gráfico (`.png`) se distribuyen bajo la licencia **Creative Commons Attribution 4.0 International (CC-BY-4.0)**.

---

##  Citación

Una vez publicado el repositorio, por favor cite el identificador persistente provisto:

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

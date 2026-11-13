# RRuta Latam: Simulación de datos y análisis de potencia

Código del proyecto **RRuta Latam** (semillero MetaCiencia, Universidad El Bosque): diagnóstico del conocimiento, la intención de uso y las barreras percibidas frente a los Reportes Registrados (RR) en Latinoamérica, mediante una encuesta anónima en línea.

Este repositorio contiene el script que se usa en la sección de tamaño de muestra del protocolo (Tabla 1 y Figura 1). Se publica antes de la recolección de datos. El código final de análisis y los datos anonimizados se depositarán después en Zenodo.

## Contenido

```
code/simulacion_analisis_metaciencia.R   Script completo (ver estructura abajo)
results/                                  Salidas del script (CSV, figuras, sessionInfo)
CITATION.cff, .zenodo.json                Metadatos para citar y para Zenodo
LICENSE                                   MIT
```

## Qué hace el script

0. Paquetes y parámetros.
1. Generador de datos (superpoblación simulada, con efectos aleatorios cruzados de país y área).
2. Depuración con las reglas definidas de antemano (duplicados, respuestas apresuradas, inconsistencias).
3. Puntajes de conocimiento e índice de experiencia investigativa.
4. Análisis confirmatorios con corrección de Holm:
   - **C1.** Experiencia y conocimiento de RR (beta-binomial mixto, `glmmTMB`).
   - **C2.** Brecha de familiaridad entre ciencia abierta y RR (logístico mixto pareado, `lme4`).
   - **C3.** Conocimiento de RR e intención de someter un RR (modelo de enlace acumulativo mixto, `ordinal`).
5. Potencia por Monte Carlo para el grid de tamaños de muestra.
6. Muestra aleatoria de la superpoblación y prueba de todos los análisis propuestos (descriptivos, exploratorios, funcionamiento diferencial de ítems, TRI 2PL y sensibilidad).

### Puntajes de conocimiento (correspondencia con el cuestionario)

| Puntaje | Ítems | Preguntas |
|---|---|---|
| Ciencia abierta | 10 | P21 (9 casillas de contenido) y P23 (selección única) |
| Reportes Registrados | 8 | P26 (6 casillas de contenido), P27 y P28 |

"Ninguna de las anteriores" y "No sé" no cuentan como ítems. P22 se reporta aparte (número de modalidades de acceso abierto reconocidas) y no entra en el puntaje.

## Requisitos

R (el protocolo reporta la versión 4.6.2) y los paquetes `tidyverse`, `lme4`, `glmmTMB`, `ordinal`, `mirt`, `parallel` (incluido en R). Opcionales: `brms` (respaldo bayesiano si el ajuste es singular) y `patchwork` (figura combinada).

```r
install.packages(c("tidyverse", "lme4", "glmmTMB", "ordinal", "mirt", "patchwork"))
```

## Cómo ejecutarlo

Desde la raíz del repositorio.

Prueba rápida (3 repeticiones por tamaño de muestra, para detectar errores):

```sh
RRUTA_NSIM=3 Rscript code/simulacion_analisis_metaciencia.R
```

Ejecución completa (1000 repeticiones por tamaño de muestra, valor por defecto):

```sh
Rscript code/simulacion_analisis_metaciencia.R
```

La primera ejecución calibra los efectos verdaderos para que coincidan con los efectos mínimos de interés (C1 y C3: β = 0,15; C2: brecha de 10 puntos porcentuales) y guarda la calibración en `results/calibracion.rds`. Si cambian el número de ítems o esos efectos, la calibración se repite sola.

## Salidas (`results/`)

- `resultados_potencia.csv` y `resultados_potencia_crudo.rds`: potencia por prueba y tamaño de muestra, y resultados de cada repetición.
- `fig_potencia.png`: distribución de los valores p (Figura 1 del protocolo).
- `fig_C1.png`, `fig_C2.png`, `fig_C3.png`, `fig_C1_C2_C3.png`: análisis de la muestra simulada.
- `sessionInfo.txt`: versiones de R y de los paquetes.

## Estado

Versión 1.0.0. El script usa el esquema de puntajes de 10 y 8 ítems y la numeración final del cuestionario. Los archivos de `results/` se generaron con esta versión, con la ejecución completa (1000 repeticiones por tamaño de muestra), y corresponden a la Tabla 1 y la Figura 1 del protocolo.

## Cómo citar

Si usas o adaptas este código, por favor cítalo:

> Leongómez, J. D., Cardona-Gómez, G. P., Cepeda, V., Lemus-Morales, A., Puentes, M., Valencia-Guevara, L., Urrego-Alvarado, D., & Castillo Morales, S. F. (2026). *RRuta Latam: Simulación de datos y análisis de potencia* (Versión 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.24069486

Los metadatos completos están en `CITATION.cff`.

## Licencia

MIT. Ver `LICENSE`. Puedes usar, modificar y redistribuir el código siempre que conserves el aviso de copyright y la licencia.

# Caracterización de la variabilidad genética intra-especie en datos de microbiota intestinal humana derivados de secuenciación con Nanoporos en poblaciones europeas

Trabajo Fin de Máster en Bioinformática — Universitat de València  
Autora: Ana María Parodi Contreras

## Descripción general

Este repositorio contiene el código fuente utilizado en el Trabajo Fin de Máster orientado a caracterizar la variabilidad genética intraespecífica en bacterias de la microbiota intestinal humana, utilizando datos de secuenciación con Oxford Nanopore Technologies del operón ribosomal completo 16S-ITS-23S.

El flujo de trabajo está organizado en tres fases principales y una etapa de análisis estadístico:

- **Fase 0:** preselección de especies bacterianas mediante mapeo global contra la base de datos GROND.
- **Fase 1:** análisis exploratorio por pools, ensamblado de novo del operón ribosomal, selección de referencias, llamada y validación de variantes candidatas.
- **Fase 2:** análisis confirmatorio a nivel individual.
- **Análisis estadístico:** caracterización de la estructura genética entre poblaciones.

## Estructura del repositorio

```text
.
├── Phase_0/
├── Phase_1/
├── Phase_2/
├── Statistical_analysis/
├── config_template.sh
├── .gitignore
└── README.md

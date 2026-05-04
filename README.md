# Caracterización de la variabilidad genética intra-especie en datos de microbiota intestinal humana derivados de secuenciación con Nanoporos en poblaciones europeas
 
Trabajo Fin de Máster en Bioinformática — Universitat de València (ETSE-UV)
 
**Autora:** Ana María Parodi Contreras  
**Tutores:** Dr. Alfonso Benítez-Páez · Dr. Vicente Arnau  
**Mayo, 2026**
 
---
 
## Descripción general
 
Este repositorio contiene el código fuente utilizado en el Trabajo Fin de Máster orientado a caracterizar la variabilidad genética intraespecífica en bacterias de la microbiota intestinal humana, utilizando datos de secuenciación con Oxford Nanopore Technologies del operón ribosomal completo 16S-ITS-23S.
 
El flujo de trabajo está organizado en tres fases principales y una etapa de análisis estadístico:
 
- **Fase 0:** preselección de especies bacterianas mediante mapeo global contra la base de datos GROND.
- **Fase 1:** análisis exploratorio por pools, ensamblado de novo del operón ribosomal, selección de referencias, llamada y validación de variantes candidatas.
- **Fase 2:** análisis confirmatorio a nivel individual.
- **Análisis estadístico:** caracterización de la estructura genética entre poblaciones.
## Cohortes
 
| País | n | Química | Acceso ENA |
|------|---|---------|------------|
| España | 113 | R9.4.1 | [PRJEB102567](https://www.ebi.ac.uk/ena/browser/view/PRJEB102567) |
| Portugal | 130 | R9.4.1 | [PRJEB74092](https://www.ebi.ac.uk/ena/browser/view/PRJEB74092) |
| Noruega | 94 | R10.3 | PRJEB112077* |
 
\* Sin acceso público aún.
 
**Base de datos de referencia:** [GROND v207](https://zenodo.org/records/10889037) (subconjunto refseq207nr_nrRep).
 
## Estructura del repositorio
 
```text
.
├── README.md
├── config.sh
├── validated_variants_blast.tsv
└── scripts/
    ├── Phase_0/
    ├── Phase_1/
    ├── Phase_2/
    └── Statistical_analysis/
```
 
## Phase_0
 
Scripts destinados al preprocesamiento de lecturas, mapeo global contra GROND, filtrado de alineamientos PAF, construcción de la tabla de abundancias y selección inicial de especies.
 
Los pasos principales de esta fase incluyen:
 
- Recorte de los primeros 100 nucleótidos del extremo 5' de cada lectura.
- Mapeo global contra la base de datos GROND.
- Filtrado de alineamientos por longitud, cobertura e identidad.
- Construcción de la tabla OTU.
- Selección de las especies bacterianas más abundantes.
- Construcción de una matriz de lecturas únicas por individuo y especie.
## Phase_1
 
Scripts destinados al análisis exploratorio por pools.
 
Esta fase incluye:
 
- Extracción de lecturas por especie y país.
- Submuestreo de lecturas para ensamblado de novo.
- Ensamblado del operón ribosomal completo mediante Flye.
- Anotación ribosomal con barrnap.
- Validación taxonómica de ensamblados.
- Selección de referencias específicas por especie.
- Construcción de pools por país y especie.
- Mapeo contra referencias específicas.
- Llamada de variantes con bcftools.
- Filtrado de variantes candidatas por calidad y recurrencia.
- Exportación de variantes y frecuencias alélicas.
- Cálculo de densidad de variantes por región funcional del operón ribosomal (16S_rRNA, ITS, 23S_rRNA).
## Phase_2
 
Scripts destinados al análisis confirmatorio a nivel individual.
 
Esta fase incluye:
 
- Selección de pares bacteria-individuo con cobertura suficiente.
- Extracción de lecturas por individuo y especie.
- Submuestreo de lecturas individuales.
- Mapeo individual contra referencias específicas por especie.
- Llamada conjunta de variantes por especie.
- Filtrado por calidad.
- Exportación de genotipos y métricas asociadas en formato TSV.
## Statistical_analysis
 
Scripts en R utilizados para la preparación de datos y análisis estadístico.
 
Esta carpeta incluye scripts para:
 
- Construcción de la matriz de variantes.
- Filtrado a variantes validadas.
- Exclusión de *Ralstonia insidiosa* del análisis principal.
- Cálculo de distancias Jaccard NA-aware.
- Análisis PERMANOVA.
- Análisis de coordenadas principales (PCoA).
- Evaluación de dispersión mediante betadisper.
- Comparaciones por pares entre países.
- Análisis por variante mediante test de Fisher.
- Corrección por FDR (Benjamini-Hochberg).
- Cálculo de frecuencias alélicas medias por país.
- Análisis integrativo mediante vegan::envfit.
- PCA exploratorio basado en frecuencias alélicas medias.
## Resultados principales
 
- **1,771** especies identificadas → **Top 50** → **13** en el análisis final
- **804** variantes candidatas → **418** validadas (BLASTn) → **105** en la matriz individual
- **PERMANOVA:** País R² = 0.076 (p ≤ 0.001); Edad R² = 0.023 (p ≤ 0.001); Sexo no significativo
- **19 variantes concordantes** (Fisher ∩ envfit, FDR < 0.05) en 4 especies: *Alistipes onderdonkii*, *Faecalibacterium prausnitzii*, *Bacteroides ovatus* y *Bacteroides uniformis*
## Software principal
 
Las herramientas principales utilizadas fueron:
 
- Cutadapt 3.5
- minimap2 2.30-r1287
- samtools 1.13
- bcftools 1.13
- Flye 2.9.6-b1802
- barrnap 0.9
- BLASTn 2.11.0+
- SeqKit 2.6.1
- GNU parallel 20210822
- R 4.5.2
Paquetes principales de R: tidyverse, vegan, ape, ggplot2, ggrepel, pheatmap, rstatix.
 
## Archivos no incluidos
 
Este repositorio no incluye: lecturas crudas de secuenciación, ficheros FASTA/FASTQ de gran tamaño, ficheros BAM/SAM, ficheros VCF intermedios, alineamientos PAF, índices de minimap2, bases de datos externas, objetos RDS ni tablas con información individual. Se incluye únicamente `validated_variants_blast.tsv` como tabla auxiliar de variantes validadas utilizada por los scripts de análisis.
 
## Reproducibilidad
 
El archivo `config.sh` contiene las rutas y parámetros del pipeline. Para ejecutar los scripts, se deben modificar las rutas con las del entorno computacional local.
Para ejecutar el flujo de trabajo se requiere: lecturas de secuenciación correspondientes, base de datos GROND, fichero de metadatos de muestras, herramientas bioinformáticas instaladas y rutas locales configuradas en `config.sh`.
 
## Contacto
 
Ana María Parodi Contreras — [apacon90](https://github.com/apacon90)
 

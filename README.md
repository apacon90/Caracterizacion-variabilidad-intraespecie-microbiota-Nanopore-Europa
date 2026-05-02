# CARACTERIZACIÓN DE LA VARIABILIDAD GENÉTICA INTRA-ESPECIE EN DATOS DE MICROBIOTA INTESTINAL HUMANA DERIVADOS DE SECUENCIACIÓN CON NANOPOROS EN POBLACIONES EUROPEAS
## Trabajo Fin de Máster en Bioinformática - Universidad de Valencia

---

##  ÍNDICE
- [Fase 0 — Preselección de especies](#fase-0--preselección-de-especies)
- [Fase 1 — Análisis exploratorio por pools](#fase-1--análisis-exploratorio-por-pools)
- [Fase 2 — Análisis confirmatorio por individuos](#fase-2--análisis-confirmatorio-por-individuos)
- [Archivos de referencia externos](#archivos-de-referencia-externos)
- [Diagrama de flujo resumido](#diagrama-de-flujo-resumido)

---

## FASE 0 - Preselección de especies
<details>
<summary><code>01_trimm.sh</code></summary>

| | Archivo | Descripción |
|---|---|---|
| **INPUT** | `*.fastq.gz` (por muestra) | Lecturas crudas de nanoporos (R9.4.1 / R10.4.1) |
| **OUTPUT** | `*.trimmed.fastq.gz` (por muestra) | Lecturas con los primeros 100 nt del extremo 5' eliminados (Cutadapt) |

</details>


<details>
<summary><code>02_map_filter.sh</code></summary>

| | Archivo | Descripción |
|---|---|---|
| **INPUT** | `*.trimmed.fastq.gz` | Lecturas preprocesadas |
| **INPUT** | `RefSeq207nr` (DB) | Base de datos de referencia |
| **OUTPUT** | `*.paf` (mapeo global) | Alineamientos minimap2 contra RefSeq207nr |
| **OUTPUT** | `*.keep.paf` | Alineamientos filtrados: longitud ≥ 3.000 nt, cobertura ≥ 70%, identidad ≥ 85% |

> Los `*.keep.paf` son el **archivo central** de la Fase 0: los consumen `03_otu_meta_top50.r`, `seleccion_14bacterias.r`, `04_reads_species.sh` (Fase 1) y `04_extract_reads_ge100_n100.sh` (Fase 2).

</details>

<details>
<summary><code>03_otu_meta_top50.r</code></summary>

|  | Archivo | Descripción |
|---|---|---|
| **INPUT** | `otu_table_top50.tsv` | Tabla OTU Top50 |
| **OUTPUT** | `pairs_14bacteria_samples_ge100.tsv` | Pares bacteria–individuo válidos: bacteria, sample_ID, n_reads (solo bacterias con prevalencia >50% a ≥100 reads) |

> `pairs_14bacteria_samples_ge100.tsv` es el **manifiesto de entrada** de toda la Fase 2.

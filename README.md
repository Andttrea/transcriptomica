# Análisis Transcriptómico

Este proyecto documenta un flujo de trabajo completo de análisis transcriptómico (RNA-seq). Incluye procesamiento de lecturas crudas, control de calidad, limpieza, y alineamiento utilizando diferentes herramientas y configuraciones (Single-end vs Paired-end) para comparar su rendimiento.

## Estructura del Proyecto

```
transcriptomica/
├── data/
│   ├── clean/                   # Datos limpios después de procesamiento con fastp
│   ├── fastqc/                  # Reportes de control de calidad inicial
│   │   └── multiqc_report_data/ # Datos consolidados de MultiQC
│	│		└── images			 # Plots de multiqc
│   ├── metadata/               # Metadatos y scripts auxiliares
│   └── SRRs/                    # Secuencias raw descargadas
├── results/
│   ├── hisat2/                  # Alineamientos con HISAT2
│   │   ├── paired_end/          # Alineamientos paired-end
│   │   └── single_end/          # Alineamientos single-end
│   ├── star/                    # Alineamientos con STAR
│   │   ├── paired_end/          # Alineamientos paired-end
│   │   └── single_end/          # Alineamientos single-end
│   └── Reporte.md               # Reporte 
├── src/                         # Scripts de análisis
│   ├── fastp_limpieza.sh        # Script de limpieza de datos
│   ├── hisat2_alineamientos.sh  # Script de alineamiento con HISAT2
│   └── star_alineamientos.sh    # Script de alineamiento con STAR
└── README.md                    # Este archivo

```

## Flujo de Trabajo (Workflow)

1. **Obtención de Datos (`extractSRR.py`)**: Descarga o extracción de los datos de secuencias a partir del SRA toolkit.
2. **Control de Calidad Inicial (FastQC / MultiQC)**: Verificación de la calidad de las lecturas crudas (`data/fastqc`).
3. **Limpieza y Recorte (`fastp_limpieza.sh`)**: Filtro de secuencias de baja calidad y remoción de adaptadores con `fastp`.
4. **Mapeo / Alineamiento (`star_alineamientos.sh`, `hisat2_alineamientos.sh`)**: Se realizaron alineamientos contra un genoma de referencia evaluando dos herramientas (STAR y HISAT2) y dos modalidades (Single-end vs Paired-end). 
5. **Evaluación y Visualización (`graficas.R`)**: Extracción de las tasas de mapeo de los archivos temporales/summary generados por los alineadores para hacer métricas de comparación visuales (p.ej., comparación de % de lecturas mapeadas).

## Dependencias

- SRA Toolkit
- FastQC y MultiQC
- fastp
- STAR
- HISAT2
- R (y paquetes como `ggplot2` o `tidyverse` para generar las gráficas)

## Reporte de Resultados

Para ver los detalles completos de este proyecto, consultar el **[Reporte Analítico](results/Reporte.md)**. 

Autor: Andrea Villarruel García
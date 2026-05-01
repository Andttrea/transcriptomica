## Planteamiento 

Para los análisis bioinformáticos la reproducibilidad y la escalabilidad son factores críticos y necesarios. Como habíamos estado haciendo anteriormente, procesar bastantes scripts individuales suele llevar a errores humanos o nos da dificultades para reanudar procesos interrumpidos.

Para resolver estos problemas se utiliza Nextflow, un marco de trabajo que esta basado en el lenguaje Groovy que permite la creación de pipelines complejos, paralelos y computacionalmente estables. Nexflow utiliza el sistema sistema de DSL2 (Domain Specific Language 2), este permite modularizar los procesos y conectarlos mediante "canales" de datos. Además, gestiona automáticamente el manejo de archivos temporales (directorio work/), permite reanudar ejecuciones desde el último punto exitoso mediante el parámetro -resume y facilita el uso de contenedores (Docker/Singularity) para garantizar que el software se ejecute siempre en la misma versión.

En este contexto, se requiere construir un pipeline reproducible y automatizado en Nextflow que procese datos de RNA-seq en formato FASTQ, integrando  de forma secuencial: un diagnóstico inicial de calidad con `FastQC`, una limpieza profunda con `Fastp`  y y un control de calidad final (volviendo a correr `FastQC`). La importancia de integrar **FastQC y fastp** dentro de un **workflow** único se centra en la necesidad de transformar el análisis de calidad de una tarea estática a un proceso dinámico. 

Todos los datos, scripts y resultados se encuentran ubicados en el servidor chaac de la licenciatura, en el directorio `/export/storage/users/andreavg/transcriptomica/nextflow_tarea`.

## Objetivos

### Objetivo general

Implementar un flujo reproducible en Nextflow para la limipieza y el control de calidad de datos RNA‑seq (FastQC + fastp).
### Objetivo particular

- Emparejar automáticamente archivos R1/R2 y procesarlos por muestra.
- Generar resultados FastQC antes y después del trimming.
- Generar reportes HTML de fastp por muestra

## Estructura del repositorio

```
transcriptomica/
├── nextflow_tarea/                    # Pipeline Nextflow principal
│   ├── src/
│   │   ├── limpieza.nf               # Script principal del workflow Nextflow
│   │   ├── .nextflow/                # Directorio de configuración de Nextflow
│   │   │   ├── history               # Historial de ejecuciones
│   │   │   ├── cache/                # Caché de compilación
│   │   │   └── plr/                  # Archivos de reproducibilidad
│   │   └── work/                     # Directorio de trabajo (archivos intermedios)
│   │       ├── 01/ca3ba312.../       # Tareas ejecutadas
│   │       ├── 02/78673fde.../
│   │       ├── 09/4b3ba9c3.../
│   │       └── ... (múltiples directorios de ejecución)
│   └── results/
│       ├── cleaned/                  # Archivos FASTQ limpios (output de fastp)
│       ├── fastqc_raw/               # Reportes FastQC de archivos crudos (14 muestras)
│       │   ├── SRR9126694_1_fastqc.html
│       │   ├── SRR9126694_2_fastqc.html
│       │   ├── SRR9126930_1_fastqc.html
│       │   ├── SRR9126930_2_fastqc.html
│       │   ├── SRR9126963_1_fastqc.html
│       │   ├── SRR9126963_2_fastqc.html
│       │   ├── SRR9127454_1_fastqc.html
│       │   ├── SRR9127454_2_fastqc.html
│       │   ├── SRR9127455_1_fastqc.html
│       │   ├── SRR9127455_2_fastqc.html
│       │   ├── SRR9127457_1_fastqc.html
│       │   ├── SRR9127457_2_fastqc.html
│       │   ├── SRR9127584_1_fastqc.html
│       │   └── SRR9127584_2_fastqc.html
│       ├── fastqc_cleaned/           # Reportes FastQC de archivos limpios (14 muestras)
│       │   ├── SRR9126694_1_fastqc.html
│       │   ├── SRR9126694_2_fastqc.html
│       │   ├── SRR9126930_1_fastqc.html
│       │   ├── SRR9126930_2_fastqc.html
│       │   ├── SRR9126963_1_fastqc.html
│       │   ├── SRR9126963_2_fastqc.html
│       │   ├── SRR9127454_1_fastqc.html
│       │   ├── SRR9127454_2_fastqc.html
│       │   ├── SRR9127455_1_fastqc.html
│       │   ├── SRR9127455_2_fastqc.html
│       │   ├── SRR9127457_1_fastqc.html
│       │   ├── SRR9127457_2_fastqc.html
│       │   ├── SRR9127584_1_fastqc.html
│       │   └── SRR9127584_2_fastqc.html
│       └── reporte/
│           ├── Reporte_Nextflow.md  # Este reporte
│           └── Reporte_Nextflow.qmd # Versión Quarto del reporte
├── deseq/                            # Sripts de alineamiento (HISAT2, STAR, Salmon) y datos generados 
│   └── [Contiene: datos crudos, datos limpios, reportes FastQC, alineamientos, conteos de features y análisis de expresión]
```

## Metodología y justificación 

### Pipeline workflow

Primero crearemos las carpeta donde se guardaran todos los datos generados y el script para el pipeline de nextflow.

El flujo de trabajo se codificó para automatizar tres fases secuenciales que transforman el análisis de calidad en un proceso dinámico y reproducible. En la primera fase, se utiliza la función fromFilePairs para agrupar los archivos FASTQ crudos por muestra y generar un diagnóstico inicial con FastQC, identificando posibles sesgos o contaminaciones desde el inicio. Posteriormente, los datos se pasan hacia la segunda fase de limpieza profunda con fastp, donde se ejecutan procesos de detección de adaptadores, recorte de colas de poliguanina y un filtrado de calidad por ventanas deslizantes (--cut_window_size 4, --cut_mean_quality 20) para garantizar la integridad de las lecturas. Finalmente, el pipeline redirige los archivos procesados hacia un segundo control con FastQC, validando la efectividad del trimming y asegurando que los datos limpios tengan la calidad óptima para las etapas posteriores de alineamiento

```bash
# Nos movemos a la carpeta donde crearemos de carpeta para nextflow
cd /export/storage/users/andreavg/transcriptomica/
mkdir nextflow_tarea
# Nos movemos a la carpeta y creamos la carpeta de src 
cd nextflow_tarea
mkdir src
cd src 
# Creamos el script de nextflow
nano limpieza.sh
```
El script es el siguiente:

```groovy
#!/usr/bin/env nextflow

// Definimos que usemos la versión 2 del DSL de Nextflow 
nextflow.enable.dsl = 2

// Definimos los parametros configurables para el pipeline
params.fastq_dir = "/export/storage/users/andreavg/transcriptomica/deseq/data/SRRs/*_{1,2}.fastq"
params.output_dir = "/export/storage/users/andreavg/transcriptomica/nextflow_tarea/results"
// Definimos los subdirectorios para el pipeline
params.fastqc_raw_dir = "${params.output_dir}/fastqc_raw"
params.cleaned_dir = "${params.output_dir}/cleaned"
params.fastqc_cleaned_dir = "${params.output_dir}/fastqc_cleaned" 

// FASE 1: FastQC para los FASTQ originales
// Definimos el proceso FastQC, que se encargará de ejecutar el análisis de calidad de los archivos FASTQ originales. 
process fastqc {
    // Etiquetamos el proceso con el ID de la muestra para facilitar el seguimiento
    tag "$sample_id"

    // Publicamos los resultados de FastQC en el directorio especificado, copiando los archivos
    // NOTA: Esto nos va a garantizar que los resultados se mantengan organizados y accesibles
    publishDir params.fastqc_raw_dir, mode: 'copy'

    // Definimos la entrada como un tuple que contiene el ID de la muestra y la ruta a los archivos FASTQ
    // NOTA: se espera recibir una tupla con un valor de ID de muestra y una ruta a los archivos FASTQ correspondientes.
    input:
    tuple val(sample_id), path(reads)

    // Definimos la salida como los archivos generados por FastQC
    output:
    path "*_fastqc.zip"
    path "*_fastqc.html"

    // El bloque de script ejecuta el comando FastQC para los archivos FASTQ de entrada,
    // generando los resultados en el directorio especificado
    script:
    """
    /export/apps/bioconda/bin/fastqc ${reads}
    """
}

// FASE 2: Limpieza de los FASTQ con FASTP

process fastp {
    // Etiquetamos el proceso con el ID de la muestra para facilitar el seguimiento
    tag "$sample_id"
    // Publicamos los resultados de FastQC en el directorio especificado, copiando los archivos
    publishDir params.cleaned_dir, mode: 'copy'
    
    // Definimos la entrada como un tuple que contiene el ID de la muestra y la ruta a los archivos FASTQ
    input: 
    tuple val(sample_id), path(reads)

    // Definimos los archivos de salida, que serán los archivos fastq generados por el fastp y el reporte html 
    // NOTA: emit se utiliza para etiquetar cada salida con un nombre específico, lo que facilita su referencia en etapas posteriores del pipeline.
    output: 
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

    // Bloque del script que correra fastp 
    script:
    // Extraemos los archivos de lectura 1 y lectura 2 de la tupla de entrada para usarlos en el comando fastp
    def (read1, read2) = reads
    // Flags:
    // -i: archivo de lectura 1
    // -I: archivo de lectura 2
    // -o: archivo de salida para lectura 1
    // -O: archivo de salida para lectura 2
    // --detect_adapter_for_pe: detecta adaptadores para datos de secuenciación pareados
    // --trim_poly_g: recorta las colas de poliguanina
    // --cut_front: recorta bases de baja calidad desde el inicio de las lecturas
    // --cut_tail: recorta bases de baja calidad desde el final de las lecturas
    // --cut_window_size: tamaño de la ventana para el recorte basado en la calidad, elegimos 4 para evaluar la calidad en ventanas de 4 bases para un recorte más preciso
    // --cut_mean_quality: calidad media requerida para el recorte, elegimos 20 para asegurar que las bases recortadas tengan una calidad media de al menos 20 para mantener la integridad de las lecturas
    // --length_required: longitud mínima requerida para las lecturas después del recorte
    // -h: archivo de salida para el reporte HTML de fastp
    """
    /export/apps/bioconda/envs/fastp/bin/fastp \\
      -i ${read1} \\
      -I ${read2} \\
      -o ${sample_id}_1.fastq.gz \\
      -O ${sample_id}_2.fastq.gz \\
      --detect_adapter_for_pe \\
      --trim_poly_g \\
      --cut_front \\
      --cut_tail \\
      --cut_window_size 4 \\
      --cut_mean_quality 20 \\
      --length_required 50 \\
      -h ${sample_id}_fastp_report.html
    """
}

// Fase 3: FastQC para los FASTQ limpios

process fastqc_cleaned {
    // Etiquetamos el proceso con el ID de la muestra para facilitar el seguimiento y 
    // publicamos los resultados de FastQC en el directorio especificado, copiando los archivos
    tag "$sample_id"
    publishDir params.fastqc_cleaned_dir, mode: 'copy'

    // Definimos la entrada como un tuple que contiene el ID de la muestra y las rutas a los archivos FASTQ limpios generados por el proceso fastp
    input:
    tuple val(sample_id), path(clean1), path(clean2)

    // Definimos la salida como los archivos generados por FastQC para los archivos FASTQ limpios
    output:
    path "*_fastqc.zip"
    path "*_fastqc.html"

    // El bloque de script ejecuta el comando FastQC para los archivos FASTQ limpios de entrada
    // Creamos el directorio de salida para los resultados de FastQC limpios si no existe y luego ejecutamos FastQC
    script:
    """
    /export/apps/bioconda/bin/fastqc ${clean1} ${clean2}
    """
}

// Definimos el workflow principal que ejecuta la ejecución de los procesos definidos anteriormente
workflow {
    // Creamos un canal que agrupa los archivos FASTQ en pares de lectura 1 y lectura 2 utilizando la función fromFilePairs de Nextflow.
    // NOTA: size 2 indica que cada par debe contener exactamente dos archivos. 
    // Broadcast se utiliza para compartir el canal entre múltiples procesos, permitiendo que tanto el proceso FastQC para los archivos originales como el proceso FastP para la limpieza accedan a los mismos pares de lectura.
    read_pairs = channel.fromFilePairs(params.fastq_dir, size: 2)
    // Ejecutamos el proceso FastQC para los archivos FASTQ originales, pasando los pares de lectura al proceso FastQC
    fastqc(read_pairs)
    // Ejecutamos el proceso FastP para limpiar los archivos FASTQ y este genera una tupla con el ID de la muestra y las rutas a los archivos FASTQ limpios, 
    cleaned_reads = fastp(read_pairs)
    // Ejecutamos el proceso FastQC para los archivos FASTQ limpios, accedemos a los archivos limpios generados por el proceso FastP a través de la tupla emitida 
    fastqc_cleaned(cleaned_reads)
}

```

### MultiQC

Con el objetivo de verificar la eficacia del pipeline, se ejecutó MultiQC para agregar los resultados de los archivos limpios. Este proceso permite inspeccionar de forma masiva cómo terminaron las secuencias tras el filtrado

```bash
# Nos movemos a la carpeta de los archivos fastqc limpios 
cd /export/storage/users/andreavg/transcriptomica/nextflow_tarea/results/fastqc_cleaned
# Creamos la carpeta para guardar los archivos generados por multiqc
mkdir multiqc
# Corremos multiqc 
multiqc ../* -o ./ 
```

### Discusión 

El análisis final mediante MultiQC muestra que el pipeline en Nextflow sí automatiza correctamente el flujo FastQC -> fastp -> FastQC, pero también permite ver con claridad qué métricas mejoraron y cuáles permanecieron prácticamente igual.

En la figura **Heatmap_datoscrudos**, los archivos sin limpiar presentan el patrón esperado de advertencias de calidad: aparecen módulos en rojo y amarillo, principalmente en **Per Base Sequence Content**, **Sequence Duplication Levels** y algunos casos en **Overrepresented Sequences**. Esto confirma que el punto de partida tenía sesgos de composición y señales de redundancia que requerían limpieza.

![Heatmap_datoscrudos](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/data/fastqc/multiqc_report_data/images/fastqc-status-check-heatmap_raw.png) 

En la figura **Heatmap_datoslimpios_manualmente**, se observa una mejora global tras el curado manual. En particular, **Adapter Content** aparece en verde y disminuyen advertencias asociadas a contaminantes técnicos. De igual manera la parte de **Per Base Sequence Content** aparece en verde, indicando que el recorte manual de los primeros nucleótidos eliminó el sesgo de composición inicial logrando una distribución de bases homogénea a lo largo de toda la lectura. 

![Heatmap_datoslimpios_manualmente](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/data/fastqc/multiqc_report_data/images/fastqc-status-check-heatmap_procesado.png)

En la figura **Heatmap_datoslimpios_nextflow**, el patrón es muy parecido: **Adapter Content** queda en verde para todas las muestras (señal de trimming efectivo de adaptadores), pero **Per Base Sequence Content** sigue en rojo. Es decir, el pipeline de Nextflow corrigió bien contaminación por adaptadores, aunque no cambió de forma fuerte el sesgo de composición por base.

![Heatmap_datoslimpios_nextflow](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/nextflow_tarea/results/fastqc_cleaned/multiqc/multiqc_data/images/fastqc-status-check-heatmap_nextflow.png)

La comparación directa entre **Per_base_sequence_SRR126694_1_manualmente** y **fastqc_per_base_sequence-SRR126694_1_nextflow** refuerza que el perfil de las muestras cambia muy poco entre ambos procesos. En ambos gráficos, la variabilidad inicial de bases en los primeros ciclos se mantiene y luego las curvas se estabilizan. Esto sugiere que, para esta muestra, el trimming aplicado no redujo de manera importante la señal de **Per Base Sequence Content** (no "trimmeó" suficientemente ese componente).

![Per_base_sequence_SRR126694_1_manualmente](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/data/fastqc/multiqc_report_data/images/fastqc_per_base_sequence-SRR126694_1.png)

![fastqc_per_base_sequence-SRR126694_1_nextflow](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/nextflow_tarea/results/fastqc_cleaned/multiqc/multiqc_data/images/fastqc_per_base_sequence-SRR126694_1_nextflow.png)

En conjunto, los resultados indican que el flujo en Nextflow fue exitoso para limpieza técnica, sobre todo adaptadores, pero conservó un sesgo de composición que probablemente está asociado a la biología y al protocolo, además de que parte de ese sesgo podría requerir un recorte más agresivo de bases iniciales si el objetivo fuera optimizar estrictamente ese módulo de FastQC. Aun así, considerando la buena calidad general por base y la eliminación de adaptadores, las lecturas siguen siendo utilizables para alineamiento y cuantificación.


### Conclusión 

El pipeline implementado en Nextflow permitió automatizar de forma reproducible el control de calidad y la limpieza de lecturas RNA-seq, cumpliendo el objetivo principal del trabajo. En términos técnicos, la estrategia fue efectiva para remover secuencias de adaptadores y mantener una calidad general por base adecuada para análisis posteriores. No obstante, el módulo **Per Base Sequence Content** permaneció con alerta en los datos procesados por Nextflow, y la comparación de perfiles por muestra indica que el sesgo inicial de composición no se corrigió de manera sustancial en todas las lecturas. Por ello, aunque los datos obtenidos son válidos para alineamiento y cuantificación, una mejora futura del flujo consistiría en evaluar un recorte inicial más específico (por ejemplo, de los primeros nucleótidos) para disminuir ese sesgo sin comprometer la longitud ni la calidad global de las secuencias.










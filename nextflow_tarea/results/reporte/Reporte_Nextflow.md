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

Primero crearemos las carpeta donde se guardaran todos los datos generados y el script para el ppipeline de nextflow

```bash
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



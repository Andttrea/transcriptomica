#!/usr/bin/env nextflow

// Definimos que usemos la versión 2 del DSL de Nextflow 
nextflow.enable.dsl = 2

// Definimos los parametros configurables para el pipeline
params.fastq_dir = '/export/storage/users/andreavg/transcriptomica/deseq/data/SRRs/*_{1,2}.fastq'
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
    path "*_{1,2}_fastqc.{zip,html}"

    // El bloque de script ejecuta el comando FastQC para los archivos FASTQ de entrada,
    // generando los resultados en el directorio especificado
    script:
    """
    mkdir -p ${params.fastqc_raw_dir}
    fastqc -o ${params.fastqc_raw_dir} ${reads}
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
    output: 
    tuple val(sample_id), path("*_{1,2}.fastq.gz")
    path "*_fastp_report.html"

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
    fastp \\
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
    path "*_{1,2}_cleaned_fastqc.{zip,html}"

    // El bloque de script ejecuta el comando FastQC para los archivos FASTQ limpios de entrada
    // Creamos el directorio de salida para los resultados de FastQC limpios si no existe y luego ejecutamos FastQC
    script:
    """
    mkdir -p ${params.fastqc_cleaned_dir}
    fastqc -o ${params.fastqc_cleaned_dir} ${clean1} ${clean2}
    """
}

// Definimos el workflow principal que ejecuta la ejecución de los procesos definidos anteriormente
workflow {
    // Creamos un canal que agrupa los archivos FASTQ en pares de lectura 1 y lectura 2 utilizando la función fromFilePairs de Nextflow.
    // NOTA: size 2 indica que cada par debe contener exactamente dos archivos. 
    read_pairs = Channel.fromFilePairs(params.fastq_dir, size: 2)
    // Ejecutamos el proceso FastQC para los archivos FASTQ originales, pasando los pares de lectura al proceso FastQC
    fastqc(read_pairs)
    // Ejecutamos el proceso FastP para limpiar los archivos FASTQ y este genera una tupla con el ID de la muestra y las rutas a los archivos FASTQ limpios, que luego se pasan al proceso FastQC para los archivos limpios
    fastp_out = fastp(read_pairs)
    // Ejecutamos el proceso FastQC para los archivos FASTQ limpios
    fastqc_cleaned(fastp_out)
}



#!/bin/bash
 
 # Asignamos la variable para el directorio de Hisat2
 SAM_Dir=("/export/storage/users/andreavg/transcriptomica/results/hisat2/single_end"
 "/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end"
 "/export/storage/users/andreavg/transcriptomica/results/star/single_end"
 "/export/storage/users/andreavg/transcriptomica/results/star/paired_end")
 
 # Recorremos cada directorio de la variable SAM_Dir y convertiremos los archivos SAM a BAM
 for DIR in "${SAM_Dir[@]}"; do
	 
     echo "Procesando directorio: $DIR"
	 
     # Recorremos cada archivo SAM en el directorio actual
	 for SAM in "$DIR"/*.sam ; do
    # Declaramos los archivos de salida BAM y el nombre base para cada archivo SAM
        base="${SAM%.sam}"
		BAM="${base}.bam"
	 
		echo "Convirtiendo $(basename $SAM)"
	 
	    #Conversión de SAM a BAM y utilizamos samtools sort para ordenar el BAM resultante. La ordenación es 
        #importante para la indexación posterior y para mejorar el rendimiento en análisis posteriores.
        # Flags:
        # -bSq 10: Convertir a BAM y filtrar por calidad de mapeo (MAPQ >= 10)
        # --threads 5: Usar 5 hilos para la conversión
        # -T : Usar un archivo temporal para la ordenación
        # -o : Especificar el archivo de salida BAM
        samtools view -bSq 10 "$SAM" | samtools sort -- threads 5 -T "${SAM}_temp" -o "$BAM"
	 
	 # Indexamos 
		samtools index "$BAM"
	 
		echo "Se convirtió con éxito $(basename $BAM)"
	 
	 done 
 done 
 
 echo "Conversión completada para todos los directorios"
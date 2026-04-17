# Reporte de pseudoalineamientos 

## Introducción 

Este reporte documenta el análisis bioinformático, incluyendo la 

Todos los datos y resultados se encuentran ubicados en el servidor chaac de la licenciatura, en el directorio `/export/storage/users/andreavg/transcriptomica`.

## Estructura del repostorio 



## Conversión de SAM a BAM 

Debido a que los archivos BAM generados por los alineadores Hisat2 y STAR pesan demasiado, lo que debemos de hacer es convertirlos a SAM. Para realizar esta tarea, utilizaremos el siguiente script. Que fue generado en la carpeta src

```bash
 # Creamos el script para hacer la conversión de SAM a BAM
 nano sam_to_bam.sh 
```

El script es el siguiente: 

```bash
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
        samtools view -bSq 10 "$SAM" | samtools sort -@ 5 -T "${base}_temp" -o "$BAM"
	 
	 # Indexamos 
		samtools index "$BAM"
	 
		echo "Se convirtió con éxito $(basename $BAM)"
	 
	 done 
 done 
 
 echo "Conversión completada para todos los directorios"
```

```bash
# Le damos permisos al script
chomd +x sam_to_bam.sh
# Ejecutamos el script
./sam_to_bam.sh
```


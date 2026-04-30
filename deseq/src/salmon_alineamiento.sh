#!/bin/bash

# Asignamos las variables para nuestro indice de referencia
index="/export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.salmon/"
# Asignamos nuestras variables para nuestros archivos de tiempo finlaes
Salmon_time_file_SE="/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/Salmon_time_SE.txt"
Salmon_time_file_PE="/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/Salmon_time_PE.txt"

# Iniciamos nuestro ciclo
for R1 in /export/storage/users/andreavg/transcriptomica/data/clean/*_1_clean.fastq.gz; do

# Asignamos el nombre de la muestra a partir del nombre del archivo

	base=$(basename $R1 _1_clean.fastq.gz)

# Definimos donde se encuentran los archivos R2

	R2="/export/storage/users/andreavg/transcriptomica/data/clean/${base}_2_clean.fastq.gz"

	echo "Procesando muestra: $base"
# Ejecutamos el comando de alineamiento con salmon para single-end
	echo "Single-end para $base"
#Flags:
# -i: Nos indica en que lugar se encuentra el index
# -l: Se utiliza este parámetro para definir las propiedades de la biblioteca de secuenciación, en este caso utilizamos `A` para que Salmon analice una parte de las lecturas para decidir por sí mismo cuál es la configuración correcta
# -p: Número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# -r: Nos indica que las lecturas son unpaired y especifica 
# -o: Indicamos donde guardaremos los resultados
# time se utiliza para medir el tiempo de ejecución del comando salmon y se redirige la salida de error estándar, que es donde se imprime el tiempo, al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando salmon y no para otros comandos que puedan estar en el script.
	{ time salmon quant -i "$index" \
		-l A -p 5 \
		-r "$R1" \
		-o "/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/${base}_SE_" ; } 2>> "$Salmon_time_file_SE" \

# Ejecutamos hisat2 con paired-end
	#echo "Paired-end para $base"
#Flags:
# -i: Nos indica en que lugar se encuentra el index
# -l: Se utiliza este parámetro para definir las propiedades de la biblioteca de secuenciación, en este caso utilizamos `A` para que Salmon analice una parte de las lecturas para decidir por sí mismo cuál es la configuración correcta
# -p: Número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# -1: indica el archivo que contiene las lecturas hacia adelante 
# -2: Indica el archivo de las lecturas reversas
# -o: Indicamos donde guardaremos los resultados
# time se utiliza para medir el tiempo de ejecución del comando salmon y se redirige la salida de error estándar, que es donde se imprime el tiempo, al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando salmon y no para otros comandos que puedan estar en el script.
	{ time salmon quant -i "$index" \
		-l A -p 5 \
		-1 "$R1" \
		-2 "$R2" \
		-o "/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/${base}_PE_" ; } 2>> "$Salmon_time_file_PE"
			

done

echo "Alineamiento completado para todas las muestras"
# Cargamos las librerias necesarias
library(GenomicFeatures)
library(tximport)
library(txdbmaker)

# Para convertir las abundancias de transcritos a genes, primero se generará una tabla de mapeo (tx2gene) 
# que relaciona cada transcrito con su gen correspondiente usando la anotación de ratón Gencode vM36.

# Primero creamos la base de datos de TxDb a partir del archivo GFF, el cual tambien utilizamos para los featurecounts
txdb <- makeTxDbFromGFF("/export/storage/users/andreavg/transcriptomica/data/gencode/gencode.vM36.chr_patch_hapl_scaff.annotation.gff3.gz")

# Extraemos los IDs de transcritos que existen en txdb
# NOTA: keytype = "TXNAME" nos indica que solo queremos los nombres de los transcritos
tx_ids <- keys(txdb, keytype = "TXNAME")

# Creamos un data frame con los IDs de transcritos y sus correspondientes genes
# NOTA: keys = tx_ids nos indica que solo queremos los IDs de transcritos que existen en txdb
# columns = "GENEID" nos indica que solo queremos los IDs de genes correspondientes a cada transcrito
tx2gene <- select(txdb, keys = tx_ids, keytype = "TXNAME", columns = "GENEID")

# Quitamos la versión de los IDs de transcritos en tx2gene
# Esto es necesario porque los quant.sf de Salmon no incluyen la versión
# NOTA: sub("\\..*", "", tx2gene$TXNAME) nos indica que se reemplace cualquier cosa después del punto (incluyendo el punto) por una cadena vacía
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)

# Vamos a crear una varible con la ruta de los archivos de abundancia de transcritos generados por Salmon
# Creamos una lista de las muestras que queremos procesar
muestras <- c("SRR9126694", "SRR9126930", "SRR9126963", "SRR9127454", "SRR9127455", "SRR9127457", "SRR9127584")

# Construimos las rutas de los archivos single-end
files_SE <- file.path("/export/storage/users/andreavg/transcriptomica/results/salmon/single_end", paste0(muestras, "_SE_"), "quant.sf")
# Asignamos los nombres de las muestras a los archivos
names(files_SE) <- muestras


# Construimos las rutas de los archivos paired-end
files_PE <- file.path("/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end", paste0(muestras, "_PE_"), "quant.sf")
# Asignamos los nombres de las muestras a los archivos
names(files_PE) <- muestras

# Ejecutamos tximport para convertir las abundancias de transcritos a genes utilizando el mapeo tx2gene

# Para los archivos single-end
# NOTA: ignoreTxVersion = TRUE nos indica que se ignore la versión de los transcritos al hacer el mapeo, 
# lo cual es útil si los nombres de los transcritos en los archivos de Salmon no incluyen la versión.
txi_SE <- tximport(files_SE, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Para los archivos paired-end
# NOTA: ignoreTxVersion = TRUE nos indica que se ignore la versión de los transcritos al hacer el mapeo, 
# lo cual es útil si los nombres de los transcritos en los archivos de Salmon no incluyen la versión.
txi_PE <- tximport(files_PE, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Guardamos los resultados en archivos
saveRDS(txi_SE, file = "/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/txi/txi_SE.rds")
saveRDS(txi_PE, file = "/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/txi/txi_PE.rds")
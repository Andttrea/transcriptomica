# Exportamos las tablas guardadas anteriormente
txi_SE <- readRDS("/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/txi/txi_SE.rds")
txi_PE <- readRDS("/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/txi/txi_PE.rds")

# Hacemos una función para convertir las tablas a formato TSV solamente de counts y abundance
rds_to_tsv <- function(txi, base_rute) {
    # Extraemos la base del nombre del archivo sin la extensión .rds
    base <- sub("\\.rds$", "", base_rute)
    # Escribimos las tablas de counts y abundance en formato TSV
    write.table(txi$counts, file = paste0(base, "_counts.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA )
    write.table(txi$abundance, file = paste0(base, "_abundance.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}

# Convertimos las tablas a formato TSV
rds_to_tsv(txi_SE,"/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/txi/txi_SE.rds")
rds_to_tsv(txi_PE,"/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/txi/txi_PE.rds")

library(dplyr)
# Cargamos el objeto txi generado por Salmon para los datos de paired-end

txi_pe<- readRDS("/export/storage/users/andreavg/transcriptomica/deseq/results/salmon/paired_end/txi/txi_PE.rds")

# Limpiamos el objeto ppara que en el geneid se quite el decimal, por ejemplo ENSG000001234.5 se convierta en ENSG000001234
rownames(txi_pe$counts) <- sub("\\..*", "", rownames(txi_pe$counts))
rownames(txi_pe$abundance) <- sub("\\..*", "", rownames(txi_pe$abundance))
rownames(txi_pe$length) <- sub("\\..*", "", rownames(txi_pe$length))

# Cargamos los metadatos tal como lo hiciste
metadata_table <- read.csv("/export/storage/users/andreavg/transcriptomica/deseq/data/metadata/GSE132040_MACA_Bulk_metadata.csv", check.names = FALSE)

# Preparamos el mapeo SRR -> Condición
mapping_data <- metadata_table %>%
  select(srr_id = "raw file", edad = "characteristics: age", sexo = "characteristics: sex") %>%
  mutate(condicion = paste0("m", sexo, "_", edad))

  # 2. Obtener los nombres actuales (SRRs)
  srr_names <- colnames(txi_pe$counts)
  
  # 3. Mapear a condiciones usando tu lógica de match
  condition_names <- mapping_data$condicion[match(srr_names, mapping_data$srr_id)]
  
  # 4. Crear nombres únicos para las réplicas
  new_names <- make.unique(condition_names, sep = "_rep")
  
  # 5. Reemplazar nombres en todas las matrices del objeto
  colnames(txi_pe$counts)    <- new_names
  colnames(txi_pe$abundance) <- new_names
  colnames(txi_pe$length)    <- new_names
  

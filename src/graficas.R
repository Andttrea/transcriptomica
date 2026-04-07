#!/usr/bin/env Rscript

# Suprimir el stdin sobre el cargado de datos
suppressPackageStartupMessages({
	library(ggplot2)
	library(dplyr)
})

# Nombre de los SRR 
srr <- c(
	"SRR9126694", "SRR9126930", "SRR9126963", "SRR9127454",
	"SRR9127455", "SRR9127457", "SRR9127584"
)

# STAR: tasa total mapeada = uniquely + multiple + too many loci
# Tomando todos los datos de STAR
# Formando un data.frame
star_single <- data.frame(
	SRR = srr,
	uniquely = c(73.60, 75.17, 73.67, 80.67, 74.63, 76.63, 78.27),
	multiple = c(14.77, 15.58, 17.20, 13.34, 16.04, 15.46, 14.32),
	too_many = c(0.83, 0.31, 0.58, 0.81, 0.40, 0.34, 1.77)
) |>
	mutate(
		Alineador = "STAR",
		Modo = "Single-end",
		Mapeo = uniquely + multiple + too_many
	) |>
	select(SRR, Alineador, Modo, Mapeo)

star_paired <- data.frame(
	SRR = srr,
	uniquely = c(67.40, 69.33, 65.95, 77.82, 70.45, 71.88, 75.03),
	multiple = c(5.89, 6.42, 7.93, 5.72, 6.55, 6.52, 7.78),
	too_many = c(0.61, 0.23, 0.40, 0.61, 0.27, 0.25, 1.38)
) |>
	mutate(
		Alineador = "STAR",
		Modo = "Paired-end",
		Mapeo = uniquely + multiple + too_many
	) |>
	select(SRR, Alineador, Modo, Mapeo)

# HISAT2: usar la tasa global de alineamiento reportada
hisat_single <- data.frame(
	SRR = srr,
	Alineador = "HISAT2",
	Modo = "Single-end",
	Mapeo = c(79.10, 81.91, 83.27, 88.13, 81.84, 84.09, 87.80)
)

hisat_paired <- data.frame(
	SRR = srr,
	Alineador = "HISAT2",
	Modo = "Paired-end",
	Mapeo = c(63.92, 66.52, 64.77, 75.25, 67.38, 69.55, 74.92)
)

# Organizar todos los datos
datos <- bind_rows(star_single, star_paired, hisat_single, hisat_paired) |>
	mutate(
		Grupo = factor(
			paste(Alineador, Modo, sep = " - "),
			levels = c(
				"STAR - Single-end", "STAR - Paired-end",
				"HISAT2 - Single-end", "HISAT2 - Paired-end"
			)
		)
	)

# Obtener el promedio de los datos
promedios <- datos |>
	group_by(Grupo) |>
	summarise(Media = mean(Mapeo), .groups = "drop")

# Formar el gráfico
grafica <- ggplot() +
	geom_col(
		data = promedios,
		aes(x = Grupo, y = Media, fill = Grupo),
		width = 0.70,
		alpha = 0.80,
		show.legend = FALSE
	) +
	geom_point(
		data = datos,
		aes(x = Grupo, y = Mapeo),
		position = position_jitter(width = 0.12, height = 0),
		size = 2.6,
		alpha = 0.95,
		color = "black"
	) +
	stat_summary(
		data = datos,
		aes(x = Grupo, y = Mapeo),
		fun.data = mean_se,
		geom = "errorbar",
		width = 0.15,
		linewidth = 0.5,
		color = "black"
	) +
	scale_y_continuous(
		limits = c(0, 100),
		expand = expansion(mult = c(0, 0.03))
	) +
	labs(
		title = "Tasa de mapeo por alineador y modo",
		subtitle = "Barras: promedio por grupo | Puntos: muestras SRR individuales",
		x = "Grupo",
		y = "Lecturas mapeadas (%)"
	) +
	theme_minimal(base_size = 12) +
	theme(
		axis.text.x = element_text(angle = 20, hjust = 1),
		panel.grid.major.x = element_blank(),
		plot.title = element_text(face = "bold")
	)

print(grafica)

ggsave(
	filename = "results/tasa_mapeo_barras_puntos.png",
	plot = grafica,
	width = 10,
	height = 6,
	dpi = 300,
    bg = "white"
)
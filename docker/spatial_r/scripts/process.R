library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(scales)

parser <- ArgumentParser(description = "Process merged RDS object by filtering and normalization")

add_argument(
	parser,
	"--team-id",
	required=TRUE,
	help="Team ID"
)
add_argument(
	parser,
	"--paths-input",
	nargs="+",
	required=TRUE,
	help="List of RDS objects to merge"
)
add_argument(
	parser,
	"--output",
	required=TRUE,
	help="Output file name for the merged RDS object"
)

args <- parser$parse_args()


#############################
## LIMIT OF QUANTIFICATION ##
#############################
target_geomxdata <- readRDS(args$input)
pkc <- annotation(target_geomxdata)
module <- gsub(".pkc", "", pkc)

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

LOQ <- data.frame(row.names = colnames(target_geomxdata))
vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
if(all(vars[1:2] %in% colnames(pData(target_geomxdata)))) {
	LOQ[, module] <-
		pmax(minLOQ,
			pData(target_geomxdata)[, vars[1]] * 
				pData(target_geomxdata)[, vars[2]] ^ cutoff)
}
pData(target_geomxdata)$LOQ <- LOQ


##################################################
## FILTER SEGMENTS AND/OR GENES WITH LOW SIGNAL ##
##################################################
LOQ_mat <- c()
ind <- fData(target_geomxdata)$Module == module
mat_i <- t(esApply(target_geomxdata[ind, ], MARGIN = 1,
                   FUN = function(x) {
                       x > LOQ[, module]
                   }))
LOQ_mat <- rbind(LOQ_mat, mat_i)
LOQ_mat <- LOQ_mat[fData(target_geomxdata)$TargetName, ]

# Save detection rate information to pheno data
pData(target_geomxdata)$GenesDetected <- colSums(LOQ_mat, na.rm = TRUE)
pData(target_geomxdata)$GeneDetectionRate <- pData(target_geomxdata)$GenesDetected / nrow(target_geomxdata)

pData(target_geomxdata)$DetectionThreshold <- cut(
	pData(target_geomxdata)$GeneDetectionRate,
	breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
	labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%")
)

segment_gene_detection_plot <- ggplot(pData(target_geomxdata), aes(x = DetectionThreshold)) +
	geom_bar(aes(fill = region)) + 
	geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
	theme_bw() +
	scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
	labs(x = "Gene Detection Rate",
		y = "Segments, #",
		fill = "Segment Type")

segment_gene_detection_plot_output = paste0(args$team_id + ".segment_gene_detection_plot.png")
ggsave(segment_gene_detection_plot_output, plot = segment_gene_detection_plot, width = 6, height = 4, dpi = 300)

# Remove segments with less than 10% of the genes detected
target_geomxdata <- target_geomxdata[, pData(target_geomxdata)$GeneDetectionRate >= .1]


saveRDS(target_geomxdata, file = args$output)

library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot) 

parser <- ArgumentParser(description = "Process merged RDS object by filtering and normalization")

parser$add_argument(
	"--cohort-id",
	required=TRUE,
	help="Cohort ID"
)
parser$add_argument(
	"--input",
	required=TRUE,
	help="The GeoMx data to filter and normalize"
)
parser$add_argument(
	"--celltype-markers",
	required=TRUE,
	help="Cell type marker list"
)
parser$add_argument(
	"--min-segment",
	required=TRUE,
	help="Minimum % of segments that detect the genes [0.1]"
)
parser$add_argument(
	"--output",
	required=TRUE,
	help="Output file name for the processed RDS object"
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

segment_gene_detection_plot_output <- paste0(args$cohort_id + ".segment_gene_detection_plot.png")
ggsave(segment_gene_detection_plot_output, plot = segment_gene_detection_plot, width = 6, height = 4, dpi = 300)

# Remove segments with less than 10% of the genes detected
target_geomxdata <- target_geomxdata[, pData(target_geomxdata)$GeneDetectionRate >= .1]

LOQ_mat <- LOQ_mat[, colnames(target_geomxdata)]
fData(target_geomxdata)$DetectedSegments <- rowSums(LOQ_mat, na.rm = TRUE)
fData(target_geomxdata)$DetectionRate <-
    fData(target_geomxdata)$DetectedSegments / nrow(pData(target_geomxdata))

cell_type_markers <- read.csv(args$celltype_markers, header = TRUE, stringsAsFactors = FALSE)
gene_list <- as.list(cell_type_markers$marker)
gene_list_df <- data.frame(
	Gene = gene_list,
	Number = fData(target_geomxdata)[gene_list, "DetectedSegments"],
	DetectionRate = percent(fData(target_geomxdata)[gene_list, "DetectionRate"]))
gene_list_df_output <- paste0(args$cohort_id + ".gene_detection_rate.csv")
write.csv(gene_list_df, gene_list_df_output, row.names = FALSE)

neg_probe_fData <- subset(fData(target_geomxdata), CodeClass == "Negative")
neg_probes <- unique(neg_probe_fData$TargetName)
target_geomxdata <- target_geomxdata[fData(target_geomxdata)$DetectionRate >= args$min_segment |
					fData(target_geomxdata)$TargetName %in% neg_probes, ]

# Retain only detected genes of interest
# gene_list <- gene_list[gene_list %in% rownames(target_geomxdata)]


###################
## NORMALIZATION ##
###################
# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "region"
stat_data <- data.frame(row.names = colnames(exprs(target_geomxdata)),
						Segment = colnames(exprs(target_geomxdata)),
						Annotation = pData(target_geomxdata)[, ann_of_interest],
						Q3 = unlist(apply(exprs(target_geomxdata), 2,
							quantile, 0.75, na.rm = TRUE)),
						NegProbe = exprs(target_geomxdata)[neg_probes, ])
stat_data_m <- melt(stat_data, measure.vars = c("Q3", "NegProbe"),
					variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(stat_data_m,
			aes(x = Value, fill = Statistic)) +
			geom_histogram(bins = 40) + theme_bw() +
			scale_x_continuous(trans = "log2") +
			facet_wrap(~Annotation, nrow = 1) + 
			scale_fill_brewer(palette = 3, type = "qual") +
			labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(stat_data,
			aes(x = NegProbe, y = Q3, color = Annotation)) +
			geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
			geom_point() + guides(color = "none") + theme_bw() +
			scale_x_continuous(trans = "log2") + 
			scale_y_continuous(trans = "log2") +
			theme(aspect.ratio = 1) +
			labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(stat_data,
			aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
			geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
			geom_point() + theme_bw() +
			scale_x_continuous(trans = "log2") + 
			scale_y_continuous(trans = "log2") +
			theme(aspect.ratio = 1) +
			labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
			rel_widths = c(0.43,0.57))
combined_plt <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
combined_plot_output <- paste0(args$cohort_id + ".q3_negprobe_plot.png")
ggsave(combined_plot_output, plot = combined_plt, width = 8, height = 6, dpi = 300)

# Q3 norm (75th percentile) for WTA/CTA with or without custom spike-ins
target_geomxdata <- normalize(
	target_geomxdata,
	norm_method = "quant",
	desiredQuantile = .75,
	toElt = "q_norm"
)

# Background normalization for WTA/CTA without custom spike-in
target_geomxdata <- normalize(
	target_geomxdata,
	norm_method = "neg",
	romElt = "exprs",
	toElt = "neg_norm"
)

saveRDS(target_geomxdata, file = args$output)

# Visualize the first 10 segments with each normalization method
normalization_plot_output <- paste0(args$cohort_id + ".normalization_plot.png")
png(normalization_plot_output, width = 800, height = 600, res = 300)

par(mfrow = c(1, 3))

boxplot(exprs(target_geomxdata)[,1:10],
		col = "#9EDAE5", main = "Raw Counts",
		log = "y", names = 1:10, xlab = "Segment",
		ylab = "Counts, Raw")

boxplot(assayDataElement(target_geomxdata[,1:10], elt = "q_norm"),
		col = "#2CA02C", main = "Q3 Norm Counts",
		log = "y", names = 1:10, xlab = "Segment",
		ylab = "Counts, Q3 Normalized")

boxplot(assayDataElement(target_geomxdata[,1:10], elt = "neg_norm"),
		col = "#FF7F0E", main = "Neg Norm Counts",
		log = "y", names = 1:10, xlab = "Segment",
		ylab = "Counts, Neg. Normalized")

dev.off()

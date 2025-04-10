library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)

parser <- ArgumentParser(description = "Perform Segment and Probe QC on GeoMx data")

parser$add_argument(
	"--sample-id",
	required=TRUE,
	help="Sample ID"
)
parser$add_argument(
	"--input",
	required=TRUE,
	help="The GeoMx data to perform QC on"
)
parser$add_argument(
	"--min-reads",
	type = "integer",
	required=TRUE,
	help="Minimum number of reads [1000]"
)
parser$add_argument(
	"--percent-trimmed",
	type = "integer",
	required=TRUE,
	help="Minimum % of reads trimmed [80]"
)
parser$add_argument(
	"--percent-stitched",
	type = "integer",
	required=TRUE,
	help="Minimum % of reads stitched [80]"
)
parser$add_argument(
	"--percent-aligned",
	type = "integer",
	required=TRUE,
	help="Minimum % of reads aligned [80]"
)
parser$add_argument(
	"--percent-saturation",
	type = "integer",
	required=TRUE,
	help="Minimum sequencing saturation [50]"
)
parser$add_argument(
	"--min-neg-count",
	type = "integer",
	required=TRUE,
	help="Minimum negative control counts [1]"
)
parser$add_argument(
	"--max-ntc-count",
	type = "integer",
	required=TRUE,
	help="Maximum counts observed in NTC well [1000]"
)
parser$add_argument(
	"--min-nuclei",
	type = "integer",
	required=TRUE,
	help="Minimum # of nuclei estimated [100]"
)
parser$add_argument(
	"--min-area",
	type = "integer",
	required=TRUE,
	help="Minimum segment area [5000]"
)
parser$add_argument(
	"--output",
	required=TRUE,
	help="Output file name for the QC'ed RDS object"
)

args <- parser$parse_args()


# Shift counts to one for downstream transformations
geomxdata <- readRDS(args$input)
geomxdata <- shiftCountsOne(geomxdata, useDALogic = TRUE)

################
## SEGMENT QC ##
################
# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html#411_Select_Segment_QC
qc_params <- list(
	minSegmentReads = args$min_reads,				# Minimum number of reads (1000)
	percentTrimmed = args$percent_trimmed,			# Minimum % of reads trimmed (80%)
	percentStitched = args$percent_stitched,		# Minimum % of reads stitched (80%)
	percentAligned = args$percent_aligned,			# Minimum % of reads aligned (80%)
	percentSaturation = args$percent_saturation,	# Minimum sequencing saturation (50%)
	minNegativeCount = args$min_neg_count,			# Minimum negative control counts (10) - chose 1
	maxNTCCount = args$max_ntc_count,				# Maximum counts observed in NTC well (1000)
	minNuclei = args$min_nuclei,					# Minimum # of nuclei estimated (100)
	minArea = args$min_area							# Minimum segment area (5000)
)

geomxdata <- setSegmentQCFlags(geomxdata, qcCutoffs = qc_params)

# Collate QC Results
qc_results <- protocolData(geomxdata)[["QCFlags"]]
flag_columns <- colnames(qc_results)
qc_summary_df <- data.frame(
	Pass = colSums(!qc_results[, flag_columns]),
	Warning = colSums(qc_results[, flag_columns])
)
qc_results$QCStatus <- apply(qc_results, 1L, function(x) {
	ifelse(sum(x) == 0L, "PASS", "WARNING")
})
qc_summary_df["TOTAL FLAGS", ] <-
	c(sum(qc_results[, "QCStatus"] == "PASS"),
	sum(qc_results[, "QCStatus"] == "WARNING"))

qc_summary_output = paste0(args$sample_id, ".segment_qc_summary.csv")
write.csv(qc_summary_df, qc_summary_output, row.names = FALSE, quote = FALSE)

# Remove flagged segments
geomxdata <- geomxdata[, qc_results$QCStatus == "PASS"]


##############
## PROBE QC ##
##############
# Generally keep the qcCutoffs parameters unchanged
# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html#421_Set_Probe_QC_Flags
geomxdata <- setBioProbeQCFlags(
	geomxdata,
	qcCutoffs = list(
		minProbeRatio = 0.1,
		percentFailGrubbs = 20
	),
	removeLocalOutliers = TRUE
)

probe_qc_results <- fData(geomxdata)[["QCFlags"]]

probe_qc_df <- data.frame(
	Passed = sum(rowSums(probe_qc_results[, -1]) == 0),
	Global = sum(probe_qc_results$GlobalGrubbsOutlier),
	Local = sum(rowSums(probe_qc_results[, -2:-1]) > 0
		& !probe_qc_results$GlobalGrubbsOutlier)
)

probe_qc_output = paste0(args$sample_id, ".probe_qc_summary.csv")
write.csv(probe_qc_df, probe_qc_output, row.names = FALSE, quote = FALSE)

# Exclude outlier probes
probe_qc_passed <-  subset(
	geomxdata,
	fData(geomxdata)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
		fData(geomxdata)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE
)

geomxdata <- probe_qc_passed


###########################
## GENE-LEVEL COUNT DATA ##
###########################
target_geomxdata <- aggregateCounts(geomxdata)
gene_count <- data.frame(exprs(target_geomxdata))
gene_count_output = paste0(args$sample_id, ".gene_count.csv")
write.csv(gene_count, gene_count_output, row.names = FALSE, quote = FALSE)

saveRDS(target_geomxdata, file = args$output)

library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(openxlsx)
library(dplyr)
library(ggforce)
library(networkD3)

parser <- ArgumentParser(description = "Convert DCC files to a NanoStringGeoMxSet object (.rds)")

parser$add_argument(
	"--team-id",
	required=TRUE,
	help="Team ID"
)
parser$add_argument(
	"--dataset-id",
	required=TRUE,
	help="Dataset ID"
)
parser$add_argument(
	"--sample-id",
	required=TRUE,
	help="Sample ID"
)
parser$add_argument(
	"--batch",
	required=TRUE,
	help="Batch"
)
parser$add_argument(
	"--dcc-dir",
	required=TRUE,
	help="Path to DCC files directory"
)
parser$add_argument(
	"--pkc-file",
	required=TRUE,
	help="Path to PKC file"
)
parser$add_argument(
	"--annotation-file",
	required=TRUE,
	help="Path to annotation file"
)
parser$add_argument(
	"--output",
	required=TRUE,
	help="Output file name for the NanoStringGeoMxSet object"
)

args <- parser$parse_args()


##########################################
## GENERATE NANOSTRING GEOMX SET OBJECT ##
##########################################
dcc_files <- list.files(path = args$dcc_dir, full.names = TRUE)

# The Sample ID in the annotation file and Sample ID in the FASTQ file name are different (dashes vs. underscores)
original_annotation_file_df <- read.xlsx(args$annotation_file)
if ("Sample_ID" %in% colnames(original_annotation_file_df)) {
	original_annotation_file_df$Fastq_Sample_ID <- gsub("-", "_", original_annotation_file_df$Sample_ID)
} else {
	stop("Column 'Sample_ID' not found in the dataset.")
}

mod_annotation_file <- "modified_annotation_file.xlsx"
write.xlsx(original_annotation_file_df, mod_annotation_file)
mod_annotation_sheet_name <- getSheetNames(mod_annotation_file)

geomxdata <- readNanoStringGeoMxSet(dccFiles = dcc_files,
									pkcFiles = args$pkc_file,
									phenoDataFile = mod_annotation_file,
									phenoDataSheet = mod_annotation_sheet_name,
									phenoDataDccColName = "Fastq_Sample_ID",
									protocolDataColNames = c("aoi", "roi"),
									experimentDataColNames = c("panel"))

# Add metadata
pData(geomxdata)$team <- args$team_id
pData(geomxdata)$dataset <- args$dataset_id
pData(geomxdata)$sample <- args$sample_id
pData(geomxdata)$batch <- args$batch
pData(geomxdata)$batch_id <- paste0(args$team_id, "_", args$dataset_id, "_", args$batch)


#####################
## SAMPLE OVERVIEW ##
#####################
sankey_cols <- c("source", "target", "value")

link1 <- count(pData(geomxdata), `slide name`, class)
link2 <- count(pData(geomxdata), class, region)
link3 <- count(pData(geomxdata), region, segment)

colnames(link1) <- sankey_cols
colnames(link2) <- sankey_cols
colnames(link3) <- sankey_cols

links <- rbind(link1,link2,link3)
nodes <- unique(data.frame(name=c(links$source, links$target)))

# sankeyNetwork is 0 based, not 1 based
links$source <- as.integer(match(links$source,nodes$name)-1)
links$target <- as.integer(match(links$target,nodes$name)-1)

sankey <- sankeyNetwork(
	Links = links,
	Nodes = nodes,
	Source = "source",
	Target = "target",
	Value = "value",
	NodeID = "name",
	units = "TWh",
	fontSize = 12,
	nodeWidth = 30
)

output_sankey_filename <- paste0(args$sample_id, ".sankey_diagram.html")
saveWidget(sankey, output_sankey_filename, selfcontained = FALSE)

saveRDS(geomxdata, file = args$output)

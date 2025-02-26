library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(dplyr)
library(ggforce)
library(networkD3)

parser <- ArgumentParser(description = "Convert DCC files to a NanoStringGeoMxSet object (.rds)")

add_argument(
	parser,
	"--team-id",
	required=TRUE,
	help="Team ID"
)
add_argument(
	parser,
	"--dataset-id",
	required=TRUE,
	help="Dataset ID"
)
add_argument(
	parser,
	"--sample-id",
	required=TRUE,
	help="Sample ID"
)
add_argument(
	parser,
	"--batch",
	required=TRUE,
	help="Batch from which the sample/dataset originated"
)
add_argument(
	parser,
	"--dcc-dir",
	required=TRUE,
	help="Path to DCC files directory"
)
add_argument(
	parser,
	"--pkc-file",
	required=TRUE,
	help="Path to PKC file"
)
add_argument(
	parser,
	"--annotation-file",
	required=TRUE,
	help="Path to annotation file"
)
add_argument(
	parser,
	"--output",
	required=TRUE,
	help="Output file name for the NanoStringGeoMxSet object"
)

args <- parser$parse_args()


##########################################
## GENERATE NANOSTRING GEOMX SET OBJECT ##
##########################################
dcc_files <- list.files(path = args$dcc_dir, full.names = TRUE)

geomxdata <- readNanoStringGeoMxSet(dccFiles = dcc_files,
									pkcFiles = args$pkc_file,
									phenoDataFile = args$annotation_file,
									phenoDataSheet = "Template", # TODO
									phenoDataDccColName = "Sample_ID",
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

output_sankey_filename <- paste0(args$team_id, "_sankey_diagram.html")
saveWidget(sankey, output_sankey_filename, selfcontained = FALSE)

saveRDS(geomxdata, file = args$output)

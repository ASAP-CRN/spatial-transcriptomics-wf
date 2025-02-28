library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

parser <- ArgumentParser(description = "Merge GeoMx data RDS objects into one")

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


geomx_data_list <- lapply(args$paths_input, readRDS)
combined_geomx_data <- do.call(combine, geomx_data_list)

saveRDS(combined_geomx_data, file = args$output)

library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(purrr)

parser <- ArgumentParser(description = "Merge GeoMx data RDS objects into one")

parser$add_argument(
	"--paths-input",
	nargs="+",
	required=TRUE,
	help="List of RDS objects to merge"
)
parser$add_argument(
	"--output",
	required=TRUE,
	help="Output file name for the merged RDS object"
)

args <- parser$parse_args()


geomx_data_list <- lapply(args$paths_input, readRDS)
combined_geomx_data  <- map_dfr(geomx_data_list, readRDS)

combined_geomx_data <- lapply(geomx_data_list, readRDS)
raster_stack_combined_geomx_data <- stack(combined_geomx_data)

saveRDS(combined_geomx_data, file = args$output)

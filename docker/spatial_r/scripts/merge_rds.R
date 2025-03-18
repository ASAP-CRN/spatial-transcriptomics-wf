library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

parser <- ArgumentParser(description = "Merge normalized GeoMx data RDS objects into one")

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


geomx_data_list <- future.apply::future_lapply(args$paths_input, readRDS)

message('Merging samples...')
combined_geomx_data <- merge(x=geomx_data_list[[1]], y=geomx_data_list[-1])
message('Done.')

saveRDS(combined_geomx_data, file = args$output)

library(argparse)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

parser <- ArgumentParser(description = "Convert processed RDS object to AnnData object")

parser$add_argument(
	"--input",
	required=TRUE,
	help="The RDS object to convert"
)
parser$add_argument(
	"--output",
	required=TRUE,
	help="Output file name for the AnnData object"
)

args <- parser$parse_args()


################
## CONVERSION ##
################
geomxdata <- readRDS(args$input)

SaveH5Seurat(geomxdata, filename = "processed.h5Seurat")
Convert("processed.h5Seurat", dest = "h5ad")

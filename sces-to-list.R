
library(SingleCellExperiment)
library(DropletUtils)
library(argparse)

parser <- ArgumentParser(description = "Quantify expression to SingleCellExperiment")

parser$add_argument('--input_dir', type='character',
                    help="Path to rds files")
parser$add_argument('--output', type = 'character', metavar = 'FILE',
                    help="Output SCE")
args <- parser$parse_args()

files <- dir(args$input_dir, 
            full.names=TRUE,
            recursive=TRUE)

sces <- lapply(files, readRDS)

# sce <- do.call('cbind', sces)

saveRDS(sces, args$output)
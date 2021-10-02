
library(SingleCellExperiment)
library(DropletUtils)
library(argparse)

parser <- ArgumentParser(description = "Quantify expression to SingleCellExperiment")

parser$add_argument('--input_dir', type='character',
                    help="Path to Cellranger output")
parser$add_argument('--id', type='character', help='Sample ID')
parser$add_argument('--output', type = 'character', metavar = 'FILE',
                    help="Output path for HLA genotypes.")
args <- parser$parse_args()

sce <- read10xCounts(args$input_dir)
sce$id <- args$id

saveRDS(sce, args$output)
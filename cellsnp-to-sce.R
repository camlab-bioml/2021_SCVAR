suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(vcfR)
  library(Matrix)
  library(glue)
  library(tidyverse)
})

ad <- readMM(snakemake@input[['ad_mat']])
dp <- readMM(snakemake@input[['dp_mat']])

cell_ids <- read_tsv(snakemake@input[['sample']], col_names = FALSE)
snp_ids <- read.vcfR(snakemake@input[['snp_id']])
cell_ids <- cell_ids$X1

baf <- ad / dp
baf[is.na(baf)] <- 0

sce <- SingleCellExperiment(assays = list(ad = ad,
                                          dp = dp,
                                          baf = baf),
                            colData = data.frame(barcode = cell_ids,
                                                 id = snakemake@wildcards[['id']]),
                                                 rowData=as.data.frame(snp_ids@fix))                            

saveRDS(sce, snakemake@output[['sce']])
# adapt from Noel Yue extract_TCGA_gendoo19.R script
# ncyue2@github:db-hgsoc-cisp/tcga-ov/microarray/extract_TCGA_gendoo19.R

# extract the TCGA data from gendoo et al., 2019 -- MetaGxOvarian package
# the data are provided as [HT_HG-U133A] Affymetrix HT Human Genome U133A Array
# https://bioconductor.org/packages/release/data/experiment/manuals/MetaGxOvarian/man/MetaGxOvarian.pdf

#BiocManager::install("MetaGxOvarian")

library(MetaGxOvarian)
library(io)

esets <- MetaGxOvarian::loadOvarianEsets()[[1]];

# micro-array data
data <- esets$TCGAOVARIAN;
exprs_data <-  exprs(data);
dim(exprs_data); # [1] 21260   536

qwrite(exprs_data, "tcga_ov_microarray_exprs.rds");

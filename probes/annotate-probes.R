# Create probe information table for filtering

# microrray: Affymetrix U133 (fRMA normalized), OVCA
# RNA-seq: TCGA-OVCA from PanCanAtlas
# TODO Expand to other cancer types

library(io)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# for production use
#library(array2rnaseq)

# for development only
library(usethis)
library(devtools)
load_all("../../array2rnaseq")


db <- org.Hs.eg.db
marr <- as.matrix(qread("../expr/micro_exprs_matched.rds"));
rseq <- as.matrix(qread("../expr/rna_seq_exprs_t_matched.rds"));

dim(marr) # 12664   294
dim(rseq) # 12664   294

# sanity check to ensure that rows and columns are matched
stopifnot(rownames(marr) == rownames(rseq))
stopifnot(colnames(marr) == colnames(rseq))

# map gene names
entrezs <- rownames(marr)
genes <- mapIds(db, keys = entrezs, keytype = "ENTREZID", column = "SYMBOL")
names(genes) <- genes
# check that there are no duplicated gene symbols
stopifnot(sum(duplicated(genes)) == 0)

rownames(marr) <- genes
rownames(rseq) <- genes

probes.info <- data.frame(
  entrez = entrezs,
  gene = genes
);
rownames(probes.info) <- NULL;

# visualize mapping

smoothScatter(marr, rseq)
plot(marr, rseq, pch = ".")

# mean expressions for each probes
rseq.m <- rowMeans(rseq)
marr.m <- rowMeans(marr)
plot(marr.m, rseq.m, pch = ".")

probes.info$marr_mean <- marr.m;
probes.info$rseq_mean <- rseq.m;


# constant genes with only one unique value
# Delete genes expressed 0 among all samples
n.uniqs <- apply(rseq, 1, function(x) length(unique(x)))
probes.info$constant <- n.uniqs == 1;

# identify genes which have high NA: pipeline deficiency?
probes.info$p_na <- apply(rseq, 1, function(z) mean(is.na(z)))
table(probes.info$p_na)

# calculate mutual information for each probes
# but we can only calculate for probes with no NA and not constant
idx <- !probes.info$constant & probes.info$p_na == 0;
probes.info$mi <- NA;
probes.info$mi[idx] <- mutual_info(marr[idx, ], rseq[idx, ]);
hist(probes.info$mi, breaks = 100)

# identify probes that are too noisy
mi.cut <- quantile(probes.info$mi, 0.15, na.rm=TRUE);
probes.info$noisy <- probes.info$mi < mi.cut;

probes.info$keep <- with(probes.info, !constant & p_na == 0 & !noisy);

head(probes.info)
table(probes.info$keep)

# write out probe annotation for all probes, indicating filter status
qwrite(probes.info, "probes_info.rds");


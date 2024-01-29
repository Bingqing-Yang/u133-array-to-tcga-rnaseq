# Learn mapping from microarray to RNA-seq

library(io)
library(ggplot2)
library(dplyr)
library(scam)
library(foreach)
library(doParallel)


# for production use
# library(array2rnaseq)

# for development only
library(devtools)


getwd()
# read in data
load_all("../../array2rnaseq")
probes <- qread("../probes/probes_info.rds");
rseq <- as.matrix(qread("../expr/rna_seq_exprs_t_matched.rds"));
marr <- as.matrix(qread("../expr/micro_exprs_matched.rds"));

# ensure that probes annotation are in correct order
stopifnot(probes$entrez == rownames(rseq))
stopifnot(probes$entrez == rownames(marr))

# annotate gene names
rownames(rseq) <- probes$gene;
rownames(marr) <- probes$gene;

probes.f <- probes[probes$keep, ];
marr.f <- marr[probes$keep, ];
rseq.f <- rseq[probes$keep, ];
J <- nrow(marr.f);
models <- select_models(marr.f, rseq.f);
table(models)

# select gene fitted on scam 
x_scam <- marr.f[which(models =='scam'), ][1:10, ]
y_scam <- rseq.f[which(models =='scam'), ][1:10, ]
probes_scam <- probes[which(models =='scam'), ][1:10, ]
models_scam <- select_models(x_scam, y_scam);
table(models_scam) # lm:0  scam: 10298
# saveRDS(x_scam, "./data/x_500.rds")
# saveRDS(y_scam, "./data/y_500.rds")


# fit models
fit_scam_r2_loess <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res2", model_type = "loess", log = TRUE)
fit_scam_r2_loess_noexp <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res2", model_type = "loess", log = FALSE)
fit_scam_r_loess <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "loess", log = TRUE)
fit_scam_r_loess_noexp <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "loess", log = FALSE)

# ---- save model ----

saveRDS(fit_scam_r2_loess, "./models/fit_scam_r2_loess_1000_lowss_aic.rds")
saveRDS(fit_scam_r2_loess_noexp, "./models/fit_scam_r2_loess_noexp_1000_lowss_aic.rds")
saveRDS(fit_scam_r_loess, "./models/fit_scam_r_loess_1000_lowss_aic.rds")
saveRDS(fit_scam_r_loess_noexp, "./models/fit_scam_r_loess_noexp_1000_lowss_aic.rds")


fit_scam_r2_loess <- readRDS("./models/fit_scam_r2_loess.rds")
fit_scam_r2_loess_noexp <- readRDS("./models/fit_scam_r2_loess_noexp.rds")
fit_scam_r_loess <- readRDS("./models/fit_scam_r_loess.rds")
fit_scam_r_loess_noexp <- readRDS("./models/fit_scam_r_loess_noexp.rds")


# ----------------

# # predict interval for scam
pred_scam_r2_loess <- predict.array2rnaseq(fit_scam_r2_loess$maps, new_data = x_scam, res_type = "res2", model_type = "loess", log = TRUE, level = 0.95)
pred_scam_r2_loess_noexp <- predict.array2rnaseq(fit_scam_r2_loess_noexp$maps, new_data = x_scam, res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
pred_scam_r_loess <- predict.array2rnaseq(fit_scam_r_loess$maps, new_data = x_scam, res_type = "res_abs", model_type = "loess", log = TRUE, level = 0.95)
pred_scam_r_loess_noexp <- predict.array2rnaseq(fit_scam_r_loess_noexp$maps, new_data = x_scam, res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)


# calculate rates that points fall in the interval
ratio_r2_loess <- rate_interval(object = fit_scam_r2_loess$maps, X = x_scam, Y = y_scam,
                              res_type = "res2", model_type = "loess", level = 0.95)
mean(ratio_r2_loess)  # 0.285102   0.7146769

ratio_r2_loess_noexp <- rate_interval(object = fit_scam_r2_loess_noexp$maps, X = x_scam, Y = y_scam,
                              res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
mean(ratio_r2_loess_noexp)  # 0.5235544   0.9515816

ratio_r_loess <- rate_interval(object = fit_scam_r_loess$maps, X = x_scam, Y = y_scam,
                                       res_type = "res_abs", model_type = "loess", log = TRUE, level = 0.95)
mean(ratio_r_loess)  # 0.285102   0.7146769

ratio_r_loess_noexp <- rate_interval(object = fit_scam_r_loess_noexp$maps, X = x_scam, Y = y_scam,
                                       res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)
mean(ratio_r_loess_noexp)  # 0.4205952   0.8885374



# models selection
x <- x_scam
y <- y_scam

aic_r2_loess <- aic_r2_loess_noexp <- aic_r_loess <- aic_r_loess_noexp <-  rep(NULL, nrow(x))

for (i in 1:nrow(x)) {
  
  aic_r2_loess[i] <- aic_output(fit_scam_r2_loess$maps[[i]]$res_model, model_type = "loess")
  aic_r2_loess_noexp[i] <- aic_output(fit_scam_r2_loess_noexp$maps[[i]]$res_model, model_type = "loess")
  aic_r_loess[i] <- aic_output(fit_scam_r_loess$maps[[i]]$res_model, model_type = "loess")
  aic_r_loess_noexp[i] <- aic_output(fit_scam_r_loess_noexp$maps[[i]]$res_model, model_type = "loess")
}

# mean: r2_loess_noexp is best
aic_lst <- c("aic_r2_loess",  "aic_r2_loess_noexp",  "aic_r_loess", "aic_r_loess_noexp" )
aic_data <- data.frame(sapply(aic_lst, get))
colMeans(aic_data)
barplot(colMeans(aic_data), main = "Model selection on AICc", ylab = "AICc", las=2)
# aic_r2_loess aic_r2_loess_noexp        aic_r_loess  aic_r_loess_noexp 
# 2.616763          -4.720887           1.230468          -3.277139


# count: r2_loess_noexp is best
table(apply(aic_data, 1, which.min)) 
# 2    4 
# 193  7 





# prepare data
aic_all <- c(aic_r2_loess, aic_r2_loess_noexp, aic_r_loess, aic_r_loess_noexp)
model_name <- rep(aic_lst, rep(200, 4))
data_test_aic <- data.frame(AIC = aic_all, Model_name = model_name)

boxplot_res_model <- ggplot(data = data_test_aic, mapping = aes(
  x = Model_name, y = AIC)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 16))
print(boxplot_res_model)
ggsave("./plots/boxplot_res_model.png",width = 10, height = 6, dpi = 300)




# non-parametric alternative to one-way ANOVA test
# https://bookdown.org/thomas_pernet/Tuto/non-parametric-tests.html
# p-value < alpha, there are significant differences between the residual models
# but do not know which pairs of models are different
kruskal.test(AIC ~ Model_name, data = data_test_aic)  

# calculate pairwise comparisons between residual models
# all pairwise comparisons are different
pairwise.wilcox.test(data_test_aic$AIC, data_test_aic$Model_name,
                     p.adjust.method = "BH")





# # generate scatter plot of diff model
# pdf(file = "fit_scam_r2_loess200.pdf")
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r2_loess, res_type = "res2", model_type = "loess", log = TRUE, level = 0.95)
#   readline("Press to continue ...")
#   print(pic)
# }
# dev.off()
# 
# pdf(file = "fit_scam_r2_loess_noexp200.pdf")
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r2_loess_noexp, res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
#   print(pic)
# }
# dev.off()
# 
# pdf(file = "fit_scam_r_loess200.pdf")
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r_loess, res_type = "res_abs", model_type = "loess", log = TRUE, level = 0.95)
#   print(pic)
# }
# dev.off()
# 
# 
# pdf(file = "fit_scam_r_loess_noexp200.pdf")
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r_loess_noexp, res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)
#   print(pic)
# }
# dev.off()



# # generate scatter plot for special gene
# # r2_loess_noexp VS r2_loess
# pic1 <- c(2)
# pic2 <- c(28,38)
# 
# pdf(file = "r2_loess_noexpVSr2_loess.pdf")
# for (i in c(pic1, pic2)) {
# 
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r2_loess_noexp, res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
#   print(pic)
# 
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r2_loess, res_type = "res2", model_type = "loess", log = TRUE, level = 0.95)
#   print(pic)
# 
# }
# dev.off()
# 
# # r2_loess_noexp VS r_loess_noexp
# pic1 <- c(36, 49, 78, 94, 134)
# pdf(file = "r2_loess_noexpVSr_loess_noexp.pdf")
# for (i in pic1) {
#   
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r2_loess_noexp, res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
#   print(pic)
#   
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r_loess_noexp, res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)
#   print(pic)
#   
# }
# dev.off()










# #############
# a <- c(909,823,406,370,631,878,135,661,659,356,187,405,763,524,819,706,388,764,917,344,225)
# b <- c(356,405,706,763,894,909)
# gene_lst <- union(a, b)
# 








# # predictions
# plot(marr.f, preds$mean, pch=".")
# 
# # calibration plot
# plot(preds$mean, rseq.f, pch=".", col=models)
# abline(a=0, b=1, col="grey")
# cor(c(preds$mean), c(rseq.f), use="complete.obs")
# 
# # examine specific genes
# gene <- 64;
# plot(marr.f[gene, ], preds$mean[gene, ], pch=".")
# points(marr.f[gene, ], preds$lower[gene, ], pch=".", col="blue")
# points(marr.f[gene, ], preds$upper[gene, ], pch=".", col="blue")
# points(marr.f[gene, ], rseq.f[gene, ], pch=".", col="red")
# 
# # modeify J ——> 100
# rs <- unlist(lapply(1:100,
#   function(j) cor(rseq.f[j, ], preds$mean[j, ])
# ));
# summary(rs)
# 
# fves <- unlist(lapply(1:100,
#   function(j) fve(rseq.f[j, ], preds$mean[j, ])
# ));
# summary(fves)
# 
# 
# saveRDS(maps, "../map/mapping.rds")
# write.table(preds, "../map/prediction_interval.txt", sep = "\t")
# write.table(models, "../map/model_type.txt", sep = "\t")








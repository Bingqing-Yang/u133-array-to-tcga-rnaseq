# Learn mapping from microarray to RNA-seq

library(io)
library(ggplot2)
library(dplyr)
library(scam)

# install.packages(c("io", "ggplot2", "scam","devtools"), dependencies = TRUE)

# for production use
# library(array2rnaseq)

# for development only
library(devtools)
load_all("../../array2rnaseq")

# read in data
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
#load_all("../../array2rnaseq")
x_scam <- marr.f[which(models =='scam'), ][1:1000, ]
y_scam <- rseq.f[which(models =='scam'), ][1:1000, ]
probes_scam <- probes[which(models =='scam'), ][1:1000, ]
models_scam <- select_models(x_scam, y_scam);
table(models_scam) # lm:0  scam: 10298
# saveRDS(x_scam, "./data/x_1000.rds")
# saveRDS(y_scam, "./data/y_1000.rds")


# # save scam model that the residual model based on lm 
fit_scam_r2_lm <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res2", model_type = "lm")
fit_scam_r_lm <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "lm")
fit_scam_r2_loess <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res2", model_type = "loess")
fit_scam_r_loess <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "loess")
# saveRDS(fit_scam_r2_lm, "./models/fit_scam_r2_lm.rds")
# saveRDS(fit_scam_r2_loess, "./models/fit_scam_r2_loess.rds")
# saveRDS(fit_scam_r_lm, "./models/fit_scam_r_lm.rds")
# saveRDS(fit_scam_r_loess, "./models/fit_scam_r_loess.rds")



fit_scam_r2_lm <- readRDS("./models/fit_scam_r2_lm.rds")
fit_scam_r2_loess <- readRDS("./models/fit_scam_r2_loess.rds")
fit_scam_r_lm <- readRDS("./models/fit_scam_r_lm.rds")
fit_scam_r_loess <- readRDS("./models/fit_scam_r_loess.rds")


# predict interval for scam
pred_scam_r2_lm <- predict.array2rnaseq(fit_scam_r2_lm$maps, new_data = x_scam, res_type = "res2", model_type = "lm")
pred_scam_r2_loess <- predict.array2rnaseq(fit_scam_r2_loess$maps, new_data = x_scam, res_type = "res2", model_type = "loess")
pred_scam_r_lm <- predict.array2rnaseq(fit_scam_r_lm$maps, new_data = x_scam, res_type = "res_abs", model_type = "lm")
pred_scam_r_loess <- predict.array2rnaseq(fit_scam_r_loess$maps, new_data = x_scam, res_type = "res_abs", model_type = "loess")

# saveRDS(pred_scam_r2_lm, "./pred/pred_scam_r2_lm.rds")
# saveRDS(pred_scam_r2_loess, "./pred/pred_scam_r2_loess.rds")
# saveRDS(pred_scam_r_lm, "./pred/pred_scam_r_lm.rds")
# saveRDS(pred_scam_r_loess, "./pred/pred_scam_r_loess.rds")


pred_scam_r2_lm <- readRDS("./pred/pred_scam_r2_lm.rds")
pred_scam_r2_loess <- readRDS("./pred/pred_scam_r2_loess.rds")
pred_scam_r_lm <- readRDS("./pred/pred_scam_r_lm.rds")
pred_scam_r_loess <- readRDS("./pred/pred_scam_r_loess.rds")

#### summary:
# 1. 不需要关注基因是否存在异质性，因为同方差其实就是异质性的一个特例；
# 2. 可以focus on 非线性的基因，肉眼看是否存在异质性，并且将其标记，记录；

x <- x_scam
y <- y_scam

R2_r2_lm <- NULL
R2_r2_loess <- NULL
R2_r_lm <- NULL
R2_r_loess <- NULL

for (i in 1:nrow(x)) {
  
  res2_model_lm <- fit_scam_r2_lm$maps[[i]]$res_model
  R2 <- summary(res2_model_lm)$r.squared
  R2_r2_lm <- c(R2_r2_lm, R2)


  res2_model_loess <- fit_scam_r2_loess$maps[[i]]$res_model
  res2_hat <- predict(res2_model_loess, x_scam[i, ])
  res2_true <- fit_scam_r2_loess$maps[[i]]$res_model$y
  R2 <- fev(res2_true ,res2_hat)
  R2_r2_loess <- c(R2_r2_loess, R2)
  
  
  res_model_lm <- fit_scam_r_lm$maps[[i]]$res_model
  R2 <- summary(res_model_lm)$r.squared
  R2_r_lm <- c(R2_r_lm, R2)
  
  
  
  res_model_loess <- fit_scam_r_loess$maps[[i]]$res_model
  res_hat <- predict(res_model_loess, x_scam[i, ])
  res_true <- fit_scam_r_loess$maps[[i]]$res_model$y
  R2 <- fev(res_true ,res_hat)
  R2_r_loess <- c(R2_r_loess, R2)
  # # 生成残差散点图
  # log_res2 <- log(model$maps[[i]]$model$residuals^2 + 1e-5)
  # # plot(x[i, ], log_res2, main = paste0("R^2: ", sprintf("%.3e", R2), "  p_value: ", sprintf("%.3e", p_value)), cex.main = 0.8)
  # 
  # # 生成残差的拟合线
  # idx <- order(x[i, ])
  # d <- data.frame(x = x[i, ][idx])
  # # lines(d$x, predict(res_model_lm, newdata = d), col = "red", lwd = 2)
  # 

  
  
  # readline("Enter to continue: ")
}

# R2 <- data.frame(
#   R2_r2_lm = R2_r2_lm,
#   R2_r2_loess = R2_r2_loess,
#   R2_r_lm = R2_r_lm,
#   R2_r_loess = R2_r_loess,
#   row.names = rownames(x)
# )

# write.csv(R2, file = "./pred/R2_summary.csv", row.names = 1)


# 1. log(res^2) VS log(|res|) which is better? ---------------------
mean(R2_r2_lm)  # 0.0114345
mean(R2_r2_loess) # 0.07986962
# mean(R2_r_lm)   # 0.01095803
# mean(R2_r_loess) # 0.07927829
## summary: 0.0114345 > 0.01096105, so log(res^2) as dependent variable is better;

idx_R0.1 <- which(R2_r2_lm > 0.1)
R2_r2_lm0.1 <- R2_r2_lm[idx_R0.1]
idx <- order(R2_r2_lm0.1, decreasing = TRUE)
idx_dec_R0.1 <- idx_R0.1[idx]
# 356 909 405 706 894 763

idx_R0.2 <- which(R2_r2_loess > 0.2)
idx <- order(R2_r2_loess[idx_R0.2], decreasing = TRUE)
idx_dec_R0.2 <- idx_R0.2[idx]
# 909 823 406 370 631 878 135 661 659 356 187 405 763 524 819 706 388 764 917 344 225


# 2. look as genes with heteroscedasticity ---------------------
length(which(R2_r2_lm > 0.1)) # 6 
length(which(R2_r2_loess > 0.2)) # 21


which(R2_r2_lm > 0.1)
# 356 405 706 763 894 909
which(R2_r2_loess > 0.2)
# 909 823 406 370 631 878 135 661 659 356 187 405 763 524 819 706 388 764 917 344 225


# look at scatter plot of genes selected by lm, R^2 > 0.1
pdf(file = "scatter_plots_genes_lm_.1.pdf")

# draw scatter plots 
for (i in idx_dec_R0.1) {
  pic <- plot(x[i, ], y[i, ], main = paste0(" Scatter plot of Gene: ", rownames(x)[i]) ,
              xlab = "Log microarray intensity", ylab = "Log RNA-seq data")
  print(pic)
}
dev.off()



# look at scatter plot of genes selected by loess, R^2 > 0.2
pdf(file = "scatter_plots_genes_loess_.2.pdf")

for (i in idx_dec_R0.2) {
  pic <- plot(x[i, ], y[i, ], main = paste0(" Scatter plot of Gene: ", rownames(x)[i]) ,
              xlab = "Log microarray intensity", ylab = "Log RNA-seq data")
  print(pic)
}
dev.off()

# 3. look at residual plot, loess would be suitable ---------------------
# # 356 405 706 763 894 909
# # 909 823 406 370 631 878 135 661 659 356 187 405 763 524 819 706 388 764 917 344 225
# # 823 388
# i = 356
# # residual plot for lm, which is monotonous
# log_res2_lm <- log(fit_scam_r2_lm$maps[[i]]$model$residuals^2 + 1e-5)
# plot(x[i, ], log_res2_lm, main = paste0("R^2: ", sprintf("%.3e", R2_r2_lm[i])), cex.main = 0.8)
# idx <- order(x[i, ])
# lines(x[i, ][idx], fit_scam_r2_lm$maps[[i]]$res_model$fitted.values[idx], col = "red", lwd = 2)
# 
# # residual plot for loess, which is not monotonous
# log_res2_loess <- log(fit_scam_r2_loess$maps[[i]]$model$residuals^2 + 1e-5)
# plot(x[i, ], log_res2_loess, main = paste0("R^2: ", sprintf("%.3e", R2_r2_lm[i])), cex.main = 0.8)
# idx <- order(x[i, ])
# lines(x[i, ][idx], fit_scam_r2_loess$maps[[i]]$res_model$fitted[idx], col = "red", lwd = 2)
# 
# 
# 
# # loess as residual model is much suitable 
# scatter(i, x_scam, y_scam, pred = pred_scam_r2_loess)
# scatter(i, x_scam, y_scam, pred = pred_scam_r2_lm)

 
pdf(file = "residual_plots_vs_plots_PI_lm_.1.pdf")
for (i in idx_dec_R0.1) {
  log_res2_lm <- log(fit_scam_r2_lm$maps[[i]]$model$residuals^2 + 1e-5)
  log_res2_loess <- log(fit_scam_r2_loess$maps[[i]]$model$residuals^2 + 1e-5)
  plot(x[i, ], log_res2_lm, main = paste0("Residual scatter plot of Gene: ", rownames(x)[i]),  
       xlab = "Log Microarray intensity", ylab = "Log square of residual in scam ",cex.main = 0.8)
  idx <- order(x[i, ])
  lines(x[i, ][idx], fit_scam_r2_lm$maps[[i]]$res_model$fitted.values[idx], col = "blue", lwd = 2)
  lines(x[i, ][idx], fit_scam_r2_loess$maps[[i]]$res_model$fitted[idx], col = "red", lwd = 2)
  legend("topright", legend = c(paste0("lm", "  R^2: ", sprintf("%.3e", R2_r2_lm[i])),  
                                paste0("loess", "  R^2: ",sprintf("%.3e", R2_r2_loess[i]))), 
         col = c("blue", "red"), lty = 1, title = "Model type", cex = 0.8)
  
  scatter(i, x_scam, y_scam, pred = pred_scam_r2_lm, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by linear model"))
    
                       
  scatter(i, x_scam, y_scam, pred = pred_scam_r2_loess, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by loess model")) 
    
}
dev.off()






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








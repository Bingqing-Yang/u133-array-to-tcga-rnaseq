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

# # select data fitted on linear model -------------------------------------
# 
# x_lm <- marr.f[which(models =='lm'), ]
# y_lm <- rseq.f[which(models =='lm'), ]
# probes_lm <- probes[which(models =='lm'), ]
# models_lm <- select_models(x_lm, y_lm);
# 
# # save model that the residual model based on lm
# fit_lm <- array2rnaseq(x_lm, y_lm, models = models_lm, residual_model = "lm")
# saveRDS(fit_lm, "./models/fit_lm.rds")
# 
# # # save model that the residual model based on loess
# fit_loess <- array2rnaseq(x_lm, y_lm, models = models_lm, residual_model = "loess")
# saveRDS(fit_loess, "./models/fit_loess.rds")
# 
# # read model
# # fit_loess <- readRDS("./model/fit_loess.rds")
# # fit_lm <- readRDS("./model/fit_lm.rds")
# 
# # predict interval for linear model
# pred_lm <- predict.array2rnaseq(fit_lm$maps, X = x_lm, new_data = x_lm)
# pred_loess <- predict.array2rnaseq(fit_loess$maps, x_lm, new_data = x_lm)
# 
# # 1. check: look at gene with heteroscedasticity ------------------------------------------------------------
# 
# hetero_idx <- fit_lm$hetero_idx
# x <- x_lm[hetero_idx, ]
# y <- y_lm[hetero_idx, ]
# library(lmtest)
# library(car)
# test_p <- data.frame()
# for (i in 1:nrow(x)) {
#   par(mfrow = c(1, 2))
#   model <- lm(y[i, ] ~ x[i, ])
#   
#   # heteroscedasticity test
#   bptest_p <- bptest(model)$p.value
#   white_p <- ncvTest(model)$p
#   test_p <- rbind(test_p, c(bptest_p, white_p))
#   
#   plot(x[i, ], y[i, ], main = paste0("bptest_p: ", sprintf("%.3e", bptest_p), "   white_p: ", sprintf("%.3e", white_p)), cex.main = 0.8)
#   plot(model, which = 1)
#   
#   # bptest manually
#   aux_lm <- lm(model$residuals^2 ~ x[i, ])
#   pchisq(ncol(x) * summary(aux_lm)$r.squared, 1, lower.tail=F)
#   
#   readline("Enter to continue: ")
# }
# 
# colnames(test_p) <- c("bptest", "whitetest")
# mean(test_p$bptest)
# mean(test_p$white)
# 
# # 2. check: which is better: log(res^2) vs log(|res|) --------------------------------------------------------
# 
# hetero_idx <- fit_lm$hetero_idx
# R2_res2 <- NULL
# R2_res <- NULL
# for (i in hetero_idx){
#   res <- lm(y_lm[i, ] ~ x_lm[i, ])$residuals
#   log_res_abs <- log(abs(res) + 1e-5)
#   log_res2 <- log(res^2 + 1e-5)
#   
#   par(mfrow = c(1, 2))
#   
#   # log(res^2)
#   residual_model_log_2 <- lm(log_res2 ~ x_lm[i, ])
#   R2 <- summary(residual_model_log_2)$r.squared
#   R2_res2 <- c(R2_res2, R2)
#   p_value <- anova(residual_model_log_2)$"Pr(>F)"[1]
#   plot(x_lm[i, ], log_res2, main = paste0("R^2: ", round(R2,8), " p: ", round(p_value, 8)))
#   abline(residual_model_log_2, col = "blue")
#   
#   # log(|res|)
#   residual_model_log_abs <- lm(log_res_abs ~ x_lm[i, ])
#   R2 <- summary(residual_model_log_abs)$r.squared
#   p_value <- anova(residual_model_log_abs)$"Pr(>F)"[1]
#   R2_res <- c(R2_res, R2)
#   
#   min_y <- min(range(log_res2), range(log_res_abs))
#   max_y <- max(range(log_res2), range(log_res_abs))
#   y_range <- c(min_y, max_y)
#   plot(x_lm[i, ], log_res_abs,ylim = y_range, main = paste0("R^2: ", round(R2,8), " p: ", round(p_value, 8)))
#   abline(residual_model_log_abs, col = "red")
# 
#   print(paste0("This is: ", i))
#   readline("Enter to contunue: ")
# }
# mean(R2_res) 
# mean(R2_res2) 


## summary:
## 1. R^2 of log(res^2) is higher ;
## 2. problem: there are log(|res|) who appear nonlinear trend; 
## 3. negative slope
# index of genes:5, 12, 14, 129, 132, 193, 204

# 3. check: lm(log(res^2) ~ x) VS loess(log(res^2) ~ x) -------------------------------------
# hetero_idx <- fit_lm$hetero_idx
# R2_lm <- NULL
# R2_loess <- NULL
# for (i in hetero_idx){
# 
#   x <- x_lm
#   y <- y_lm
#   
#   model <- lm(y[i, ] ~ x[i, ])
#   res <- model$residuals
#   log_res2 <- log(res^2 + 1e-5)
#   
#   par(mfrow = c(1, 2))
#   
#   # residual model based on linear model
#   plot(x[i, ], log_res2)
#   res2_model_lm <- fit_lm$maps[[i]]$res_model$res_model  
#   idx <- order(x[i, ])
#   d <- data.frame(x = x[i, ][idx])
#   y_hat_lm <- predict(res2_model_lm, d)
#   lines(d$x, y_hat_lm, col = "blue")
#   
#   R2 <- summary(res2_model_lm)$r.squared
#   R2_lm <- c(R2_lm, R2)
#   
#   
#   # residual model based on loess
#   plot(x[i, ], log_res2)
#   res2_model_loess <- fit_loess$maps[[i]]$res_model$res_model
#   y_hat_loess <- predict(res2_model_loess, d)
#   lines(d$x, y_hat_loess, col = "red")
# 
#   
#   R2_pseudo <- fev(y[i, ], y_hat_loess)
#   R2_loess <- c(R2_loess, R2_pseudo)
#   
#   readline("Press to continue: ")
# 
# }
# # compare R^2 between lm and loess
# mean(R2_loess) # -42.47051
# mean(R2_lm) # 0.0365966

# SUMMARY:
# 1. R^2 of loess is less;


# 4. log(|res|) scatter plot based on lm or loess -----------------------------------------------

# library(gridExtra)
# 
# for (i in 1:dim(x_lm)[1]){
#   pic1 <- scatter(i, x_lm[i, ] , y_lm[i, ], pred = pred_lm)
#   pic2 <- scatter(i, x_lm[i, ], y_lm[i, ], pred = pred_loess)
#   grid.arrange(pic1, pic2, ncol = 2)
#   readline("Enter to contunue: ")
# }
## summary:
## 实际上两者肉眼看没有什么区别，但是数据有差异
## 生成的PI也是直线，不是曲线
## 曲线好还是直线好呢？

# scam ------------------------------------------------------------------------------------------
# select gene fitted on scam 
#load_all("../../array2rnaseq")
x_scam <- marr.f[which(models =='scam'), ][1:1000, ]
y_scam <- rseq.f[which(models =='scam'), ][1:1000, ]
probes_scam <- probes[which(models =='scam'), ][1:1000, ]
models_scam <- select_models(x_scam, y_scam);
table(models_scam) # lm:0  scam: 10298
saveRDS(x_scam, "./data/x_1000.rds")
saveRDS(y_scam, "./data/y_1000.rds")


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

saveRDS(pred_scam_r2_lm, "./pred/pred_scam_r2_lm.rds")
saveRDS(pred_scam_r2_loess, "./pred/pred_scam_r2_loess.rds")
saveRDS(pred_scam_r_lm, "./pred/pred_scam_r_lm.rds")
saveRDS(pred_scam_r_loess, "./pred/pred_scam_r_loess.rds")


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

R2_summary <- data.frame(
  R2_r2_lm = R2_r2_lm,
  R2_r2_loess = R2_r2_loess,
  R2_r_lm = R2_r_lm,
  R2_r_loess = R2_r_loess,
  row.names = rownames(x)
)
# 
write.csv(R2_summary, file = "./pred/R2_summary.csv", row.names = TRUE)

# 

res2_model_lm <- fit_scam_r2_lm$maps[[i]]$res_model











# 1. log(res^2) VS log(|res|) which is better? ---------------------
mean(R2_r2_lm)  # 0.0114345
mean(R2_r2_loess) # 0.07986962
# mean(R2_r_lm)   # 0.01095803
# mean(R2_r_loess) # 0.07927829
R2 <- read.csv("./pred/R2_summary.csv", row.names = 1)
## summary: 0.0114345 > 0.01096105, so log(res^2) as dependent variable is better;

idx_R0.1 <- which(R2_r2_loess > 0.2)
R2_r2_lm0.1 <- R2_r2_loess[idx_R0.1]
idx <- order(R2_r2_lm0.1, decreasing = TRUE)
idx_dec_R0.1 <- idx_R0.1[idx]
# 356 909 405 706 894 763
# 909 823 406 370 631 878 135 661 659 356 187 405 763 524 819 706 388 764 917 344 225


# 2. look as genes with heteroscedasticity ---------------------
length(which(R2_r2_lm > 0.1)) # 6 
length(which(R2_r2_loess > 0.2)) # 21


which(R2_r2_lm > 0.1)
# 356 405 706 763 894 909
which(R2_r2_loess > 0.2)
# 909 823 406 370 631 878 135 661 659 356 187 405 763 524 819 706 388 764 917 344 225


i = 370
R2_r2_loess[i]
# raw data scatter plot
plot(x[i, ], y[i, ], main = paste0("R^2: ", sprintf("%.3e", R2_r2_loess[i]), "  Gene index: ", i))
idx <- order(x[i, ])
lines(x[i, ][idx], fit_scam_r2_loess$maps[[i]]$model$fitted[idx], col = "red", lwd = 2)

library(gridExtra)


# Open pdf file 
pdf(file = "scatter_plots_genes_loess_0.2.pdf")

# draw plots 

for (i in idx_dec_R0.1) {
  pic <- plot(x[i, ], y[i, ], main = paste0(" Scatter plot of Gene: ", rownames(x)[i]) ,
              xlab = "Log microarray intensity", ylab = "Log RNA-seq data")
  print(pic)
}
dev.off()

# 3. look at residual plot, loess would be suitable ---------------------
# 356 405 706 763 894 909
# 909 823 406 370 631 878 135 661 659 356 187 405 763 524 819 706 388 764 917 344 225
# 823 388
i = 356
# residual plot for lm, which is monotonous
log_res2_lm <- log(fit_scam_r2_lm$maps[[i]]$model$residuals^2 + 1e-5)
plot(x[i, ], log_res2_lm, main = paste0("R^2: ", sprintf("%.3e", R2_r2_lm[i])), cex.main = 0.8)
idx <- order(x[i, ])
lines(x[i, ][idx], fit_scam_r2_lm$maps[[i]]$res_model$fitted.values[idx], col = "red", lwd = 2)

# residual plot for loess, which is not monotonous
log_res2_loess <- log(fit_scam_r2_loess$maps[[i]]$model$residuals^2 + 1e-5)
plot(x[i, ], log_res2_loess, main = paste0("R^2: ", sprintf("%.3e", R2_r2_lm[i])), cex.main = 0.8)
idx <- order(x[i, ])
lines(x[i, ][idx], fit_scam_r2_loess$maps[[i]]$res_model$fitted[idx], col = "red", lwd = 2)



# loess as residual model is much suitable 
scatter(i, x_scam, y_scam, pred = pred_scam_r2_loess)
scatter(i, x_scam, y_scam, pred = pred_scam_r2_lm)

#


# draw plots 
pdf(file = "residual_plots_vs_plots_PI_6genes.pdf")
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



load_all("../../array2rnaseq")


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








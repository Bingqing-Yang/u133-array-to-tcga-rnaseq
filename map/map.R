# Learn mapping from microarray to RNA-seq

library(io)
library(ggplot2)
library(dplyr)
library(scam)
library(foreach)
library(doParallel)

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
x_scam <- marr.f[which(models =='scam'), ][1:200, ]
y_scam <- rseq.f[which(models =='scam'), ][1:200, ]
probes_scam <- probes[which(models =='scam'), ][1:200, ]
models_scam <- select_models(x_scam, y_scam);
table(models_scam) # lm:0  scam: 10298
saveRDS(x_scam, "./data/x_500.rds")
saveRDS(y_scam, "./data/y_500.rds")




# # parallel computing
# detectCores()
# cl <- makeCluster(6)
# registerDoParallel(cl)
# # stop parallel
# stopCluster(cl)

# source("../../array2rnaseq/R/predict.R")
# source("../../array2rnaseq/R/models.R")
# source("../../array2rnaseq/R/functions.R")
# source("../../array2rnaseq/R/array2rnaseq.R")
# export_func <- c("residual_model", "array2rnaseq", "lm_map", "mutual_info", "predict_interval",
#                  "residual_plot", "scam_map", "scam_select_k", "scatter", "select_models", "res_fit",
#                  "loess_model", "r_square", "aicc", "predict.residual_model"
#                  )
# 
# export_vari <- c("x_scam", "y_scam")



# fit models
fit_scam_r2_lm <- array2rnaseq(X = x_scam, Y = y_scam, models = models_scam)
fit_scam_r2_loess <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res2", model_type = "loess", log = TRUE, export_vari, export_func)
fit_scam_r2_lm_noexp <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res2", model_type = "lm", log = FALSE, export_vari, export_func)
fit_scam_r2_loess_noexp <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res2", model_type = "loess", log = FALSE, export_vari, export_func)


fit_scam_r_lm <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "lm", log = TRUE, export_vari, export_func)
fit_scam_r_loess <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "loess", log = TRUE, export_vari, export_func)
fit_scam_r_lm_noexp <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "lm", log = FALSE, export_vari, export_func)
fit_scam_r_loess_noexp <- array2rnaseq(x_scam, y_scam, models = models_scam, res_type = "res_abs", model_type = "loess", log = FALSE, export_vari, export_func)

# ---- save model ----
saveRDS(fit_scam_r2_lm, "./models/fit_scam_r2_lm.rds")
saveRDS(fit_scam_r2_loess, "./models/fit_scam_r2_loess.rds")
saveRDS(fit_scam_r2_lm_noexp, "./models/fit_scam_r2_lm_noexp.rds")
saveRDS(fit_scam_r2_loess_noexp, "./models/fit_scam_r2_loess_noexp.rds")
saveRDS(fit_scam_r_lm, "./models/fit_scam_r_lm.rds")
saveRDS(fit_scam_r_loess, "./models/fit_scam_r_loess.rds")
saveRDS(fit_scam_r_lm_noexp, "./models/fit_scam_r_lm_noexp.rds")
saveRDS(fit_scam_r_loess_noexp, "./models/fit_scam_r_loess_noexp.rds")





# fit_scam_r2_lm <- readRDS("./models/fit_scam_r2_lm_1000.rds")
# fit_scam_r2_loess <- readRDS("./models/fit_scam_r2_loess_1000.rds")
# fit_scam_r_loess_noexp <- readRDS("./models/fit_scam_r_loess_noexp_1000.rds")
# fit_scam_r_lm <- readRDS("./models/fit_scam_r_lm.rds")
# fit_scam_r_loess <- readRDS("./models/fit_scam_r_loess.rds")

# ----------------

# # predict interval for scam
# pred_scam_r2_lm <- predict.array2rnaseq(fit_scam_r2_lm$maps, new_data = x_scam, level = 0.95)
# pred_scam_r2_loess <- predict.array2rnaseq(fit_scam_r2_loess$maps, new_data = x_scam, res_type = "res2", model_type = "Gam", log = TRUE, level = 0.95)
# pred_scam_r2_lm_noexp <- predict.array2rnaseq(fit_scam_r2_lm_noexp$maps, new_data = x_scam, res_type = "res2", model_type = "lm", log = FALSE, level = 0.95)
# pred_scam_r2_loess_noexp <- predict.array2rnaseq(fit_scam_r2_loess_noexp$maps, new_data = x_scam, res_type = "res2", model_type = "Gam", log = FALSE, level = 0.95)
# 
# pred_scam_r_lm <- predict.array2rnaseq(fit_scam_r_lm$maps, new_data = x_scam, res_type = "res_abs", model_type = "lm", log = TRUE, level = 0.95)
# pred_scam_r_loess <- predict.array2rnaseq(fit_scam_r_loess$maps, new_data = x_scam, res_type = "res_abs", model_type = "Gam", log = TRUE, level = 0.95)
# pred_scam_r_lm_noexp <- predict.array2rnaseq(fit_scam_r_lm_noexp$maps, new_data = x_scam, res_type = "res_abs", model_type = "lm", log = FALSE, level = 0.95)
# pred_scam_r_loess_noexp <- predict.array2rnaseq(fit_scam_r_loess_noexp$maps, new_data = x_scam, res_type = "res_abs", model_type = "Gam", log = FALSE, level = 0.95)
# 
# # ---- save pred ----
# # saveRDS(pred_scam_r2_lm, "./pred/pred_scam_r2_lm.rds")
# # saveRDS(pred_scam_r2_loess, "./pred/pred_scam_r2_loess.rds")
# # saveRDS(pred_scam_r_lm, "./pred/pred_scam_r_lm.rds")
# # saveRDS(pred_scam_r_loess, "./pred/pred_scam_r_loess.rds")
# 
# # pred_scam_r2_lm <- readRDS("./pred/pred_scam_r2_lm.rds")
# # pred_scam_r2_loess <- readRDS("./pred/pred_scam_r2_loess.rds")
# # pred_scam_r_lm <- readRDS("./pred/pred_scam_r_lm.rds")
# # pred_scam_r_loess <- readRDS("./pred/pred_scam_r_loess.rds")
# # ----
# 
# # calculate probability that points fall in the interval for different model
# ratio_r2_lm <- ratio_interval(object = fit_scam_r2_lm$maps, X = x_scam, Y = y_scam, level = 0.95)
# mean(ratio_r2_lm)  # 0.2822109     0.7108163
# 
# ratio_r2_loess <- ratio_interval(object = fit_scam_r2_loess$maps, X = x_scam, Y = y_scam, 
#                               res_type = "res2", model_type = "loess", level = 0.95)
# mean(ratio_r2_loess)  # 0.2819388   0.7147619
# 
# ratio_r2_lm_noexp <- ratio_interval(object = fit_scam_r2_lm_noexp$maps, X = x_scam, Y = y_scam, 
#                               res_type = "res2", model_type = "lm", log = FALSE, level = 0.95)
# mean(ratio_r2_lm_noexp)   # 0.5211565   0.9481633
# 
# ratio_r2_loess_noexp <- ratio_interval(object = fit_scam_r2_loess_noexp$maps, X = x_scam, Y = y_scam, 
#                               res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
# mean(ratio_r2_loess_noexp)  # 0.5203401   0.9515306
# 
# ratio_r_loess <- ratio_interval(object = fit_scam_r_loess$maps, X = x_scam, Y = y_scam, 
#                                        res_type = "res_abs", model_type = "loess", log = TRUE, level = 0.95)
# mean(ratio_r_loess)  # 0.5203401   0.9515306
# 
# ratio_r_loess_noexp <- ratio_interval(object = fit_scam_r_loess_noexp$maps, X = x_scam, Y = y_scam, 
#                                        res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)
# mean(ratio_r_loess_noexp)  # 0.5203401   0.9515306
# 
# 
# 
# 
# 
# #### summary:
# # 1. 不需要关注基因是否存在异质性，因为同方差其实就是异质性的一个特例；
# # 2. 可以focus on 非线性的基因，肉眼看是否存在异质性，并且将其标记，记录；
# 
# 
# # calculate R^2 and AIC for residual model
# 
# x <- x_scam
# y <- y_scam
# 
# R_r2_lm <- R_r2_loess <- R_r2_lm_noexp <- R_r2_loess_noexp <- rep(NULL, nrow(x))
# R_r_lm <- R_r_loess <- R_r_lm_noexp <- R_r_loess_noexp <- rep(NULL, nrow(x))
# 
# aic_r2_lm <- aic_r2_loess <- aic_r2_lm_noexp <- aic_r2_loess_noexp <- rep(NULL, nrow(x))
# aic_r_lm <- aic_r_loess <- aic_r_lm_noexp <- aic_r_loess_noexp <- rep(NULL, nrow(x))
# 
# 
# for (i in 1:nrow(x)) {
#   R_r2_lm[i] <- R2_output(i, y, fit_scam_r2_lm$maps[[i]]$res_model, model_type = "lm")
#   aic_r2_lm[i] <- aic_output(fit_scam_r2_lm$maps[[i]]$res_model, model_type = "lm")
# 
#   
#   R_r2_loess[i] <- R2_output(i, y, fit_scam_r2_loess$maps[[i]]$res_model, model_type = "loess")
#   aic_r2_loess[i] <- aic_output(fit_scam_r2_loess$maps[[i]]$res_model, model_type = "loess")
#   
#   R_r2_lm_noexp[i] <- R2_output(i, y, fit_scam_r2_lm_noexp$maps[[i]]$res_model, model_type = "lm")
#   aic_r2_lm_noexp[i] <- aic_output(fit_scam_r2_lm_noexp$maps[[i]]$res_model, model_type = "lm")
#   
#   R_r2_loess_noexp[i] <- R2_output(i, y, fit_scam_r2_loess_noexp$maps[[i]]$res_model, model_type = "loess")
#   aic_r2_loess_noexp[i] <- aic_output(fit_scam_r2_loess_noexp$maps[[i]]$res_model, model_type = "loess")
#   
#   R_r_lm[i] <- R2_output(i, y, fit_scam_r_lm$maps[[i]]$res_model, model_type = "lm")
#   aic_r_lm[i] <- aic_output(fit_scam_r_lm$maps[[i]]$res_model, model_type = "lm")
# 
#   R_r_loess[i] <- R2_output(i, y, fit_scam_r_loess$maps[[i]]$res_model, model_type = "loess")
#   aic_r_loess[i] <- aic_output(fit_scam_r_loess$maps[[i]]$res_model, model_type = "loess")
# 
#   R_r_lm_noexp[i] <- R2_output(i, y, fit_scam_r_lm_noexp$maps[[i]]$res_model, model_type = "lm")
#   aic_r_lm_noexp[i] <- aic_output(fit_scam_r_lm_noexp$maps[[i]]$res_model, model_type = "lm")
# 
#   R_r_loess_noexp[i] <- R2_output(i, y, fit_scam_r_loess_noexp$maps[[i]]$res_model, model_type = "loess")
#   aic_r_loess_noexp[i] <- aic_output(fit_scam_r_loess_noexp$maps[[i]]$res_model, model_type = "loess")
# 
# }
# 
# 
# R_lst <- c("R_r2_lm", "R_r2_loess", "R_r2_lm_noexp", "R_r2_loess_noexp", "R_r_lm", "R_r_loess", "R_r_lm_noexp", "R_r_loess_noexp")
# R2 <- sapply(R_lst, get)
# R2_data <- data.frame(R2)
# barplot(colMeans(R2_data), main = "Model comparison on Rsquared", ylab = "Rsquared", las=2) 
# colMeans(R2_data)
# # hist(R2_data$R_r2_loess)
# # write.csv(R2_data, file = "./pred/R2_data.csv", row.names = 1)
# 
# aic_lst <- c("aic_r2_lm", "aic_r2_loess", "aic_r2_lm_noexp", "aic_r2_loess_noexp", "aic_r_lm", "aic_r_loess", "aic_r_lm_noexp", "aic_r_loess_noexp" )
# aic <- sapply(aic_lst, get)
# aic_data <- data.frame(aic)
# colMeans(aic_data)
# # barplot(colMeans(aic_data), main = "Model comparison on AICc", ylab = "AICc", las=2) 
# 
# # aic_r2_lm       aic_r2_loess    aic_r2_lm_noexp aic_r2_loess_noexp           aic_r_lm 
# # 1317.817532           2.645420        -687.200294          -4.176701         910.246829 
# # aic_r_loess     aic_r_lm_noexp  aic_r_loess_noexp 
# # 1.259126        -342.089156          -3.002049 
# 
# # summary: 
# # 1. exclude linear model, because of monotonicity;
# # 2. for loess model, r2_loess_noexp better than r_loess_noexp based on AIC;
# 
# 
# 
# 
# # generate scatter plot of diff model
# pdf(file = "fit_scam_r2_loess.pdf") 
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r2_loess, res_type = "res2", model_type = "loess", log = TRUE, level = 0.95)
#   print(pic)
# }
# dev.off()
# 
# pdf(file = "fit_scam_r2_loess_noexp.pdf")
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r2_loess_noexp, res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
#   print(pic)
# }
# dev.off()
#   
# pdf(file = "fit_scam_r_loess.pdf")
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r_loess, res_type = "res_abs", model_type = "loess", log = TRUE, level = 0.95)
#   print(pic)
# }
# dev.off()
# 
# 
# pdf(file = "fit_scam_r_loess_noexp.pdf")
# for (i in 1:nrow(x_scam)) {
#   pic <- combine_scatter(i, x_scam, y_scam, fit_scam_r_loess_noexp, res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)
#   print(pic)
# }
# dev.off()
# 
# 
# 
# # generate scatter plot for special gene
# pic1 <- c(2,15,16)
# pic2 <- c(23,28,36,38)
# 
# pdf(file = "r2_loessVSr2_loess_noexp.pdf")
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
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #############
# a <- c(909,823,406,370,631,878,135,661,659,356,187,405,763,524,819,706,388,764,917,344,225)
# b <- c(356,405,706,763,894,909)
# gene_lst <- union(a, b)
# 
# # look at scatter plot of genes selected by lm, R^2 > 0.1
# pdf(file = "scatter_lm.1.pdf")
# 
# # draw scatter plots 
# #for (i in idx_dec_R0.1) {
# for (i in 1:dim(x)[1]) {
#   pic <- plot(x[i, ], y[i, ], main = paste0(" Scatter plot of Gene: ", rownames(x)[i]) ,
#               xlab = "Log microarray intensity", ylab = "Log RNA-seq data")
#   print(pic)
# }
# dev.off()
# 
# 
# 
# 
# # 3. look at residual plot, loess would be suitable ---------------------
# # # 356 405 706 763 894 909
# # # 909 823 406 370 631 878 135 661 659 356 187 405 763 524 819 706 388 764 917 344 225
# # # 823 388
# # i = 356
# # # residual plot for lm, which is monotonous
# # log_res2_lm <- log(fit_scam_r2_lm$maps[[i]]$model$residuals^2 + 1e-5)
# # plot(x[i, ], log_res2_lm, main = paste0("R^2: ", sprintf("%.3e", R2_r2_lm[i])), cex.main = 0.8)
# # idx <- order(x[i, ])
# # lines(x[i, ][idx], fit_scam_r2_lm$maps[[i]]$res_model$fitted.values[idx], col = "red", lwd = 2)
# # 
# # # residual plot for loess, which is not monotonous
# # log_res2_loess <- log(fit_scam_r2_loess$maps[[i]]$model$residuals^2 + 1e-5)
# # plot(x[i, ], log_res2_loess, main = paste0("R^2: ", sprintf("%.3e", R2_r2_lm[i])), cex.main = 0.8)
# # idx <- order(x[i, ])
# # lines(x[i, ][idx], fit_scam_r2_loess$maps[[i]]$res_model$fitted[idx], col = "red", lwd = 2)
# # 
# # 
# # 
# # # loess as residual model is much suitable 
# # scatter(i, x_scam, y_scam, pred = pred_scam_r2_loess)
# # scatter(i, x_scam, y_scam, pred = pred_scam_r2_lm)
# ###
# pdf(file = "111.pdf")
# for (i in 1:dim(x_scam)[1]) {
#   i = 1
#   y_lab <- log(fit_scam_r2_lm$maps[[i]]$model$residuals^2)
#   x_lab <- x_scam[i, ]
#   
#   plot(x_lab, y_lab, main = paste0("Residual scatter plot of Gene: ", rownames(x)[i]),
#        xlab = "Log Microarray intensity", ylab = "Log residual square in scam", cex.main = 0.8)
#   idx <- order(x_lab)
#   c <- fit_scam_r2_lm$maps[[i]]$C
#   
#   # 红色lm(log(res2) ~ x)
#   res_model <- fit_scam_r2_lm$maps[[i]]$res_model
#   d <- data.frame(x = x_lab)
#   res_pred <- predict(res_model, d)
#   lines(x_lab[idx], res_pred[idx], col = "red")
#   
#   # 绿色loess(log(res2) ~ x)
#   res_model <- fit_scam_r2_loess$maps[[i]]$res_model
#   res_pred <- predict(res_model, d)
#   lines(x_lab[idx], res_pred[idx], col = "green")
#   
#   # 蓝色loess(res2 ~ x)
#   res_model <- fit_scam_r2_loess_noexp$maps[[i]]$res_model
#   res_pred <- predict(res_model, d)
#   lines(x_lab[idx], log(res_pred[idx] + c) , col = "blue")
#   
#   
#   legend("topright", legend = c("lm(log(res2) ~ x)" , "loess(log(res2) + c ~ x)", "loess(res2 ~ x)"),
#           col = c("red", "green", "blue"), lty = 1, title = "Model type", cex = 0.8)
#   
#   scatter(i, x_scam, y_scam, pred = pred_scam_r2_lm, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by linear model"))
#   scatter(i, x_scam, y_scam, pred = pred_scam_r2_loess, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by linear model"))
#   scatter(i, x_scam, y_scam, pred = pred_scam_r2_loess_noexp, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by linear model"))
#  
#   
#   #readline("Enter...") 
# }
# dev.off() 
# 
# 
# 
# 
# 
# 
# 
# 
# pdf(file = "1233.pdf")
# # for (i in idx_dec_R0.2) {
# hist(R2_r2_lm, breaks = 100, main = "Histogram of R^2 of lm", xlab = "R^2", cex.main = 0.8)
# hist(R2_r2_loess, breaks = 100, main = "Histogram of R^2 of loess", xlab = "R^2", cex.main = 0.8)
# hist(R2_r_loess_noexp, breaks = 100, main = "Histogram of R^2 of loess without exp", xlab = "R^2", cex.main = 0.8)
# for (i in 1:dim(x)[1]) {
#   C <- fit_scam_r2_lm$maps[[i]]$C
#   log_res2_lm <- log(fit_scam_r2_lm$maps[[i]]$model$residuals^2 + C)
#   log_res2_loess <- log(fit_scam_r2_loess$maps[[i]]$model$residuals^2 + C)
#   res_abs_pred <- fit_scam_r_loess_noexp$maps[[i]]$res_model$fitted
#   
#   
#   
#   
# 
#   idx <- order(x[i, ])
#   # log_lm
#   lines(x[i, ][idx], fit_scam_r2_lm$maps[[i]]$res_model$fitted.values[idx], col = "blue", lwd = 2)
#   # log_loess
#   lines(x[i, ][idx], fit_scam_r2_loess$maps[[i]]$res_model$fitted[idx], col = "red", lwd = 2)
#   # log
#   lines(x[i, ][idx], log(res_abs_pred^2 + C)[idx], col = "orange", lwd = 2)
#   legend("topright", legend = c(paste0("lm", "  R^2: ", sprintf("%.3e", R2_r2_lm[i])),  
#                                 paste0("loess", "  R^2: ",sprintf("%.3e", R2_r2_loess[i])),
#                                 paste0("loess_no_exp", "  R^2: ",sprintf("%.3e", R2_r_loess_noexp[i]))),
#          col = c("blue", "red", "orange"), lty = 1, title = "Model type", cex = 0.8)
# 
#   
#   plot(x[i, ], abs(fit_scam_r2_lm$maps[[i]]$model$residuals), main = paste0("Residual scatter plot of Gene: ", rownames(x)[i]),
#        xlab = "Log Microarray intensity", ylab = "Absolute value of residual in scam", cex.main = 0.8)
#   lines(x[i, ][idx], res_abs_pred[idx], col = "orange", lwd = 2)
#   
#   scatter(i, x_scam, y_scam, pred = pred_scam_r2_lm, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by linear model"))
#     
#   scatter(i, x_scam, y_scam, pred = pred_scam_r2_loess, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by loess model"))
# 
#   scatter(i, x_scam, y_scam, pred = pred_scam_r_loess_noexp, tittle = paste0("Scatter plot of gene: ", rownames(x)[i], "\n Prediction interval was generate by loess model"))   
# }
# dev.off()
# 
# rmse <- aic <- NULL
# for (i in 1:100) {
#   rmse <-  append(fit_scam_r2_loess$maps[[i]]$span_rmse, rmse)
#   aic <- append(fit_scam_r2_loess$maps[[i]]$span_aic, aic)
# }
# 
# 
# cond <- which(aic < rmse)
# cond
# length(cond)


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








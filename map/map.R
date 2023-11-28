# Learn mapping from microarray to RNA-seq

library(io)
library(ggplot2)
library(dplyr)
library(scam)

# for production use
#library(array2rnaseq)

# for development only
library(devtools)


load_all("../../array2rnaseq/R")

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
# subset data-star --------------------------------------------------------

# load_all("../../array2rnaseq/R")
# # models <- models[1:1000]
# marr.f <- marr.f[1:1000, ]
# rseq.f <- rseq.f[1:1000, ]
# probes <- probes[1:1000, ]

# 只筛选出来适合lm模型的数据
x <- marr.f[which(models =='lm'), ]
y <- rseq.f[which(models =='lm'), ]
probes_lm <- probes[which(models =='lm'), ]
models_lm <- select_models(x, y);

# save model that the residual model based on lm
fit <- array2rnaseq(x, y, models = models_lm)
saveRDS(fit, "./fit_ls.rds")

# save model that the residual model based on loess
fit_loess <- array2rnaseq(x, y, models = models_lm)
saveRDS(fit_loess, "./fit_loess.rds")

pred <- predict(fit$maps, x)
pred_loess <- predict(fit_loess$maps, x)

############ 将异方差数据plot， 观察是log(res^2) 和log(|res|)哪个好？
hetero_idx <- fit$hetero_idx
for (i in hetero_idx){
  res <- lm(y[i, ] ~ x[i, ])$residuals
  log_res_abs <- log(abs(res) + 1e-5)
  log_res2 <- log(res^2 + 1e-5)
  
  par(mfrow = c(1, 2))
  
  residual_model_log_2 <- lm(log_res2 ~ x[i, ])
  R2 <- summary(residual_model_log_2)$r.squared
  p_value <- anova(residual_model_log_2)$"Pr(>F)"[1]
  plot(x[i, ], log_res2, main = paste0("R^2: ", round(R2,8), " p: ", round(p_value, 8)))
  abline(residual_model_log_2, col = "blue")
  
  
  residual_model_log_abs <- lm(log_res_abs ~ x[i, ])
  R2 <- summary(residual_model_log_abs)$r.squared
  p_value <- anova(residual_model_log_abs)$"Pr(>F)"[1]
  # y_range <- par()$usr[3:4]
  min_y <- min(range(log_res2), range(log_res_abs))
  max_y <- max(range(log_res2), range(log_res_abs))
  y_range <- c(min_y, max_y)
  
  plot(x[i, ], log_res_abs,ylim = y_range, main = paste0("R^2: ", round(R2,8), " p: ", round(p_value, 8)))
  abline(residual_model_log_abs, col = "red")
  
  print(paste0("This is: ", i))
  readline("Enter to contunue: ")
}
## summary:
## 1. log(|res|) 比 log(res^2) 好
## 2. problem: there are log(|res|) who appear nonlinear trend; 
# index of genes:5, 12, 14, 129, 132, 193, 204

########## loess VS lm : lm(log(|res|) ~ x) VS loess(log(|res|) ~ x)
hetero_idx <- fit$hetero_idx
for (i in hetero_idx){
  model <- lm(y[i, ] ~ x[i, ])
  res <- model$residuals
  log_res_abs <- log(abs(res) + 1e-5)
  
  res_model <- fit$maps[[i]]$loess$res_model
  
  par(mfrow = c(1, 2))
  plot(x[i, ], log_res_abs)
  idx <- order(x[i, ])
  y_hat_loess <- predict(res_model, x[i, ][idx])
  lines(x[i, ][idx], y_hat_loess, col = "blue")
  
  
  plot(x[i, ], log_res_abs)
  data_res_array <- data.frame(arraydata = x[i, ], log_res_abs = log_res_abs)
  residual_model_log_abs <- lm(log_res_abs ~ arraydata, data = data_res_array)
  abline(residual_model_log_abs, col = "red")

  readline("Press to continue: ")

  }


######

## 看一下散点图based on lm vs loess of residual model 是否有差异
pred <- predict(fit$maps, x)
library(gridExtra)
for (i in 1:dim(x)[1]){
  pic1 <- scatter(i, x, y, pred = pred)
  pic2 <- scatter(i, x, y, pred = pred_loess)
  grid.arrange(pic1, pic2, ncol = 2)
  readline("Enter to contunue: ")
}
## summary:
## 实际上两者肉眼看没有什么区别，但是数据有差异
## 生成的PI也是直线，不是曲线
## 曲线好还是直线好呢？


models_types <- select_models(marr.f, rseq.f);

# Fit models
fits <- array2rnaseq(marr.f, rseq.f, models = models_types)

# generate prediction interval
pred <- predict(fits$maps, marr.f)


# subset data-end ---------------------------------------------------------



x_scam <- marr.f[which(models =='scam'), ][1:10, ]
y_scam <- rseq.f[which(models =='scam'), ][1:10, ]
probes_scam <- probes[which(models =='scam'), ][1:10, ]
models_scam <- select_models(x_scam, y_scam);

# save model that the residual model based on lm
fit_scam <- array2rnaseq(x_scam, y_scam, models = models_scam)
saveRDS(fit, "./fit_ls.rds")

# save model that the residual model based on loess
fit_loess <- array2rnaseq(x, y, models = models_lm)
saveRDS(fit_loess, "./fit_loess.rds")

pred_scam <- predict(fit_scam$maps, x_scam)
pred_loess <- predict(fit_loess$maps, x)




















# predictions
plot(marr.f, preds$mean, pch=".")

# calibration plot
plot(preds$mean, rseq.f, pch=".", col=models)
abline(a=0, b=1, col="grey")
cor(c(preds$mean), c(rseq.f), use="complete.obs")

# examine specific genes
gene <- 64;
plot(marr.f[gene, ], preds$mean[gene, ], pch=".")
points(marr.f[gene, ], preds$lower[gene, ], pch=".", col="blue")
points(marr.f[gene, ], preds$upper[gene, ], pch=".", col="blue")
points(marr.f[gene, ], rseq.f[gene, ], pch=".", col="red")

# modeify J ——> 100
rs <- unlist(lapply(1:100,
  function(j) cor(rseq.f[j, ], preds$mean[j, ])
));
summary(rs)

fves <- unlist(lapply(1:100,
  function(j) fve(rseq.f[j, ], preds$mean[j, ])
));
summary(fves)


saveRDS(maps, "../map/mapping.rds")
write.table(preds, "../map/prediction_interval.txt", sep = "\t")
write.table(models, "../map/model_type.txt", sep = "\t")







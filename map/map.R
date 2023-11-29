# Learn mapping from microarray to RNA-seq

library(io)
library(ggplot2)
library(dplyr)
library(scam)

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

# select data fitted on linear model -------------------------------------

x_lm <- marr.f[which(models =='lm'), ]
y_lm <- rseq.f[which(models =='lm'), ]
probes_lm <- probes[which(models =='lm'), ]
models_lm <- select_models(x_lm, y_lm);

# save model that the residual model based on lm
fit_lm <- array2rnaseq(x_lm, y_lm, models = models_lm, residual_model = "lm")
saveRDS(fit_lm, "./model/fit_lm.rds")

# # save model that the residual model based on loess
fit_loess <- array2rnaseq(x_lm, y_lm, models = models_lm, residual_model = "loess")
saveRDS(fit_loess, "./model/fit_loess.rds")

# read model
# fit_lm <- readRDS("./model/fit_lm.rds")
# fit_loess <- readRDS("./model/fit_loess.rds")

# predict interval for linear model
pred_lm <- predict.array2rnaseq(fit_lm$maps, X = x_lm, new_data = x_lm)
pred_loess <- predict.array2rnaseq(fit_loess$maps, x_lm, new_data = x_lm)


# check: which is better: log(res^2) vs log(|res|) --------------------------------------------------------

hetero_idx <- fit_lm$hetero_idx

for (i in hetero_idx){
  res <- lm(y_lm[i, ] ~ x_lm[i, ])$residuals
  log_res_abs <- log(abs(res) + 1e-5)
  log_res2 <- log(res^2 + 1e-5)
  
  par(mfrow = c(1, 2))
  
  residual_model_log_2 <- lm(log_res2 ~ x_lm[i, ])
  R2 <- summary(residual_model_log_2)$r.squared
  p_value <- anova(residual_model_log_2)$"Pr(>F)"[1]
  plot(x_lm[i, ], log_res2, main = paste0("R^2: ", round(R2,8), " p: ", round(p_value, 8)))
  abline(residual_model_log_2, col = "blue")
  
  
  residual_model_log_abs <- lm(log_res_abs ~ x_lm[i, ])
  R2 <- summary(residual_model_log_abs)$r.squared
  p_value <- anova(residual_model_log_abs)$"Pr(>F)"[1]
  
  min_y <- min(range(log_res2), range(log_res_abs))
  max_y <- max(range(log_res2), range(log_res_abs))
  y_range <- c(min_y, max_y)
  
  plot(x_lm[i, ], log_res_abs,ylim = y_range, main = paste0("R^2: ", round(R2,8), " p: ", round(p_value, 8)))
  abline(residual_model_log_abs, col = "red")

  print(paste0("This is: ", i))
  readline("Enter to contunue: ")
}



## summary:
## 1. log(|res|) is better than log(res^2);
## 2. problem: there are log(|res|) who appear nonlinear trend; 
## 3. negative slope
# index of genes:5, 12, 14, 129, 132, 193, 204

# check: lm(log(|res|) ~ x) VS loess(log(|res|) ~ x) -------------------------------------
hetero_idx <- fit_lm$hetero_idx
for (i in hetero_idx){
  i = 2
  x <- x_lm
  y <- y_lm
  
  model <- lm(y[i, ] ~ x[i, ])
  res <- model$residuals
  log_res_abs <- log(abs(res) + 1e-5)
  
  par(mfrow = c(1, 2))
  
  # residual model based on linear model
  plot(x[i, ], log_res_abs)
  res_model_lm <- fit_lm$maps[[i]]$res_model$res_model  ###  改成res
  idx <- order(x[i, ])
  d <- data.frame(x = x[i, ][idx])
  y_hat_lm <- predict(res_model_lm, d)
  lines(d$x, y_hat_lm, col = "blue")
  
  # residual model based on loess
  plot(x[i, ], log_res_abs)
  res_model_loess <- fit_loess$maps[[i]]$res_model$res_model
  y_hat_loess <- predict(fit_loess$maps[[i]]$res_model$res_model, d)
  lines(d$x, y_hat_loess, col = "red")

  readline("Press to continue: ")

}


# log(|res|) scatter plot based on lm or loess -----------------------------------------------

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




# scam -----------------------------------------------
load_all("../../array2rnaseq")
x_scam <- marr.f[which(models =='scam'), ][1:200, ]
y_scam <- rseq.f[which(models =='scam'), ][1:200, ]
probes_scam <- probes[which(models =='scam'), ][1:200, ]
models_scam <- select_models(x_scam, y_scam);

# save model that the residual model based on lm
fit_scam_lm <- array2rnaseq(x_scam, y_scam, models = models_scam,  residual_model = "lm")
saveRDS(fit_scam_lm, "./model/fit_scam_lm.rds")

# save model that the residual model based on loess
fit_scam_loess <- array2rnaseq(x_scam, y_scam, models = models_scam, residual_model = "loess")
saveRDS(fit_scam_loess, "./model/fit_scam_loess.rds")

pred_scam_lm <- predict.array2rnaseq(fit_scam_lm$maps, X = x_scam, new_data = x_scam)
pred_scam_loess <- predict.array2rnaseq(fit_scam_loess$maps, X = x_scam, new_data = x_scam)


for (i in 1:dim(x_scam)[1]){
  pic1 <- scatter(i, x_scam, y_scam, pred = pred_scam_lm)
  # pic2 <- scatter(i, x_scam, y_scam, pred = pred_scam_loess)
  # grid.arrange(pic1, pic2, ncol = 2)
  readline("Enter to contunue: ")
}
















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








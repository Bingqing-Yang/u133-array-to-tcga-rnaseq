library(io)
library(ggplot2)
library(dplyr)
library(scam)
library(foreach)


# install.packages(c("io", "ggplot2", "scam","devtools"), dependencies = TRUE)

# for production use
# library(array2rnaseq)

# for development only
library(devtools)



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


# select gene fitted on lm
x_lm <- marr.f[which(models =='lm'), ][1:10, ]
y_lm <- rseq.f[which(models =='lm'), ][1:10, ]
probes_lm <- probes[which(models =='lm'), ][1:10, ]
models_lm <- select_models(x_lm, y_lm);
table(models_lm)



# fit model
lm_r2_loess <- array2rnaseq(x_lm, y_lm, models = models_lm, 
                            res_type = "res2", model_type = "loess", 
                            log = TRUE)
lm_r2_loess_noexp <- array2rnaseq(x_lm, y_lm, models = models_lm, 
                                  res_type = "res2", model_type = "loess", 
                                  log = FALSE)
lm_r_loess <- array2rnaseq(x_lm, y_lm, models = models_lm, 
                           res_type = "res_abs", model_type = "loess", 
                           log = TRUE)
lm_r_loess_noexp <- array2rnaseq(x_lm, y_lm, models = models_lm, 
                                 res_type = "res_abs", model_type = "loess", 
                                 log = FALSE)


# prediction interval
pred_lm_r2_loess <- predict.array2rnaseq(lm_r2_loess$maps, new_data = x_lm, 
                                         res_type = "res2", model_type = "loess", 
                                         log = TRUE, level = 0.95)

pred_lm_r2_loess_noexp <- predict.array2rnaseq(lm_r2_loess_noexp$maps, new_data = x_lm, 
                                               res_type = "res2", model_type = "loess", 
                                               log = FALSE, level = 0.95)

pred_lm_r_loess <- predict.array2rnaseq(lm_r_loess$maps, new_data = x_lm, 
                                        res_type = "res_abs", model_type = "loess", 
                                        log = TRUE, level = 0.95)

pred_lm_r_loess_noexp <- predict.array2rnaseq(lm_r_loess_noexp$maps, new_data = x_lm, 
                                              res_type = "res_abs", model_type = "loess", 
                                              log = FALSE, level = 0.95)









# calculate the rate of points fall in PI
ratio_r2_loess <- rate_interval(object = lm_r2_loess$maps, X = x_lm, Y = y_lm,
                                res_type = "res2", model_type = "loess", level = 0.95)
mean(ratio_r2_loess)  # 0.3742358   0.7992767

ratio_r2_loess_noexp <- rate_interval(object = lm_r2_loess_noexp$maps, X = x_lm, Y = y_lm,
                                      res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
mean(ratio_r2_loess_noexp)  # 0.5271104   0.95373

ratio_r_loess <- rate_interval(object = lm_r_loess$maps, X = x_lm, Y = y_lm,
                               res_type = "res_abs", model_type = "loess", log = TRUE, level = 0.95)
mean(ratio_r_loess)  # 0.3742358   0.7992767

ratio_r_loess_noexp <- rate_interval(object = lm_r_loess_noexp$maps, X = x_lm, Y = y_lm,
                                     res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)

# generate scatter plot of diff model
dir.create("./plots/")
pdf(file = "./plots/lm_r2_loess.pdf")
for (i in 1:nrow(x_lm)) {
  pic <- combine_scatter(i, x_lm, y_lm, lm_r2_loess, res_type = "res2", model_type = "loess", log = TRUE, level = 0.95)
  print(pic)
}
dev.off()


pdf(file = "./plots/lm_r2_loess_noexp.pdf")
for (i in 1:nrow(x_lm)) {
  pic <- combine_scatter(i, x_lm, y_lm, lm_r2_loess_noexp, res_type = "res2", model_type = "loess", log = FALSE, level = 0.95)
  print(pic)
}
dev.off()

pdf(file = "./plots/lm_r_loess.pdf")
for (i in 1:nrow(x_lm)) {
  pic <- combine_scatter(i, x_lm, y_lm, lm_r_loess, res_type = "res_abs", model_type = "loess", log = TRUE, level = 0.95)
  print(pic)
}
dev.off()

pdf(file = "lm_r_loess_noexp.pdf")
for (i in 1:nrow(x_lm)) {
  pic <- combine_scatter(i, x_lm, y_lm, lm_r_loess_noexp, res_type = "res_abs", model_type = "loess", log = FALSE, level = 0.95)
  print(pic)
}
dev.off()




# models selection
x <- x_lm
y <- y_lm

aic_r2_loess <- aic_r2_loess_noexp <- aic_r_loess <- aic_r_loess_noexp <-  rep(NULL, nrow(x))

for (i in 1:nrow(x)) {
  if (isTRUE(lm_r2_loess$maps[[i]]$hetero)) {
    aic_r2_loess[i] <- aic_output(lm_r2_loess$maps[[i]]$res_model, model_type = "loess")
    aic_r2_loess_noexp[i] <- aic_output(lm_r2_loess_noexp$maps[[i]]$res_model, model_type = "loess")
    aic_r_loess[i] <- aic_output(lm_r_loess$maps[[i]]$res_model, model_type = "loess")
    aic_r_loess_noexp[i] <- aic_output(lm_r_loess_noexp$maps[[i]]$res_model, model_type = "loess")
    
  }
}

# mean: r2_loess_noexp is best
aic_lst <- c("aic_r2_loess",  "aic_r2_loess_noexp",  "aic_r_loess", "aic_r_loess_noexp" )
aic_data <- data.frame(sapply(aic_lst, get))
colMeans(aic_data, na.rm = TRUE)
barplot(colMeans(aic_data, na.rm = TRUE), main = "Model selection on AICc", ylab = "AICc", las=2)
# aic_r2_loess aic_r2_loess_noexp        aic_r_loess  aic_r_loess_noexp 
# 2.628583          -4.196380           1.242288          -3.149621 



# count: r2_loess_noexp is best
na_idx <- which(is.na(aic_data$aic_r2_loess)) # there are 83 points satify homoscedasticity
aic_data <- aic_data[-na_idx, ]
dim(aic_data) # 154   4
table(apply(aic_data, 1, which.min))
# 2   4 
# 141  13





# prepare data
aic_all <- c(aic_r2_loess, aic_r2_loess_noexp, aic_r_loess, aic_r_loess_noexp)
model_name <- rep(aic_lst, rep(nrow(aic_data), 4))
aic_plot <- data.frame(AIC = unlist(aic_data), Model_name = model_name)

# boxplot
anova_result <- ggplot(data = aic_plot, mapping = aes(
  x = Model_name, y = AIC)) +
  geom_boxplot() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 20))
print(anova_result)


# analysis of variance (ANOVA)
result <- aov(AIC ~ Model_name, data = aic_plot)
summary(result)
TukeyHSD(result)

# do not satisfy normality and equal variance
# can not use anova analysis
# https://statsandr.com/blog/anova-in-r/#whats-next
shapiro.test(aic_plot$AIC)  # not normality
bartlett.test(AIC ~ Model_name, data = aic_plot) # not equal variance

# non-parametric alternative to one-way ANOVA test
# https://bookdown.org/thomas_pernet/Tuto/non-parametric-tests.html
# but do not know which pairs of models are different
kruskal.test(AIC ~ Model_name, data = aic_plot)  # p-value < alpha

# calculate pairwise comparisons between residual models
pairwise.wilcox.test(aic_plot$AIC, aic_plot$Model_name,
                     p.adjust.method = "BH")



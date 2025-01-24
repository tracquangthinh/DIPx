library(ggplot2)
library(reshape2)
library(ggprism)
library(RColorBrewer)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(svglite)
library(ggsci)
library(patchwork)
library(magrittr)

#figure 5a
load("data/challenge_data_with_drugcomb.RData")
target_genes =readRDS("data/target_genes_DREAM.RDS")

n_overlapped_target = sapply(seq(nrow(use.drug)), function(i) {
  drugA = use.drug$COMPOUND_A[i]
  drugB = use.drug$COMPOUND_B[i]
  length(intersect(target_genes[[drugA]], target_genes[[drugB]]))
})
is_overlapped = sapply(n_overlapped_target, function(s) if(s) 1 else 0)

synergy_per_group = split(use.drug$SYNERGY_SCORE, is_overlapped)

df = melt(synergy_per_group)
colnames(df) = c("synergy", "is_overlap")

df = df[-which(df$synergy > 230), ]
t.test(synergy_per_group[[1]], synergy_per_group[[2]])

t_test_df <- data.frame(
  group1 = 1, group2 = 2,
  label = c("p-val > 0.72"),
  y.position = c(150)
)


plt_5a <- ggplot(df) + 
  geom_violin(aes(x = is_overlap, y = synergy, fill = is_overlap), trim = FALSE,  outlier.shape = NA) +
  stat_summary(data = df, aes(x = is_overlap, y = synergy), fun=median, geom="point", size=2, color="black") +
  theme_classic() +
  scale_fill_brewer(palette="Dark2") +
  theme(legend.position = "none") +
  labs(y = "Drug synergy - Loewe score", x = "Target genes between drug A and drug B") +
  scale_x_discrete(labels = list("Non-overlapped (9659)", "Overlapped (495)")) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  add_pvalue(t_test_df, tip.length = 0, label.size = 7) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ylim(-150, 150)

plt_5a

#####################################################################################################################
# figure 5b
load("data/use.drug_DREAM_full_16May2022.Rdata")
z_score = readRDS("data/drug_vs_drug_zscore.RDS")
target_genes =readRDS("data/target_genes_DREAM.RDS")

drugs = rownames(z_score)

x = sapply(seq(nrow(use.drug)), function(i) {
  d1 = use.drug$COMPOUND_A[i]
  d2 = use.drug$COMPOUND_B[i]
  if(d1 %in% drugs & d2 %in% drugs){
    # return(drug_sim[d1, d2])
    return(z_score[d1, d2])
  } else{
    return(-100)
  }
})

use.drug$z_score = x
pick = which(use.drug$z_score == -100)
use.drug = use.drug[-pick, ]

n_overlapped_target = sapply(seq(nrow(use.drug)), function(i) {
  drugA = use.drug$COMPOUND_A[i]
  drugB = use.drug$COMPOUND_B[i]
  length(intersect(target_genes[[drugA]], target_genes[[drugB]]))
})
is_overlapped = sapply(n_overlapped_target, function(s) if(s) 1 else 0)
use.drug$is_overlap = is_overlapped


th = quantile(use.drug$z_score, c(0, 0.25, 0.5, 0.75, 1))
cut_sim <- cut(use.drug$z_score, th, include.lowest = TRUE)
t = split(use.drug$SYNERGY_SCORE, cut_sim)

df = melt(t)
colnames(df) = c("synergy", "group")
df$group = factor(df$group, c("[-2,-0.33]", "(-0.33,-0.141]", "(-0.141,2.97]", "(2.97,29.3]"))

t_test_df <- data.frame(
  group1 = c(1, 2, 3, 2), group2 = c(2, 3, 4, 4),
  label = c("1.7 * 10^-7", "0.9", "5.7 * 10^-8", 
            paste0("p-value < ", "7.8 * 10^-7")),
  y.position = c(23, 25, 27, 29)
)
# t_test_df <- df  %>% t_test(synergy ~ group) %>% add_xy_position(x = "group", dodge = 0.8)


plt_5b <- ggplot(df) + 
  geom_boxplot(aes(x = group, y = synergy, fill = group), trim = FALSE,  outlier.shape = NA) +
  theme_classic() +
  scale_fill_brewer(palette="Set2") +
  theme(legend.position = "none") +
  labs(y = "Drug synergy - Loewe score", x = "z-score between target genes of drug A and drug B") +
  # scale_x_discrete(labels = list("Non-overlapped (9945)", "Overlapped (516)")) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  theme(plot.caption = element_text(size = 22)) +
  # add_pvalue(t_test_df, tip.length = 0, label.size = 7, paste = TRUE) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  coord_cartesian(ylim=c(-5, 30))+
  labs(caption = parse(text = "Lowest~vs~Highest~group:~p-val < 7.8 %*% 10^-7"))

plt_5b

################################################################################
# figure 5c
x = readRDS("data/PAS_functional_inter.RDS")
t1 = tapply(x$pred, x$COMBINATION_ID, median)
t2 = tapply(x$SYNERGY_SCORE, x$COMBINATION_ID, median)

df = data.frame(pred = t1, obs = t2)

cor.test(t1, t2, method = "spearman")

plt_5c <- ggplot(df, aes(
    x = pred,
    y = obs)) +
    geom_point(shape = 21, alpha = 0.9, size = 3, fill = "#EBF5D8", color = "#5A7B29", stroke = 2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Median(Prediction) per combination", y = "Median(Observation) per combination") +
    theme(plot.title = element_text(size = 26)) +
    theme(axis.title = element_text(size = 26)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    geom_abline(slope = 1, intercept = 0) +
    annotate(geom="text", x = 14, y = -1, label=parse(text = "r == 0.46 (p-val < 2.2 %*% 10^-16)"), size = 9) +
    ggtitle("Functional interaction z-score ~ PAS")
plt_5c


plt_5 <- ggarrange(plt_5a, plt_5b, plt_5c,
  nrow = 2, ncol = 2, common.legend = F,
  labels = c("a", "b", "c") , font.label = list(size = 30)
)
ggsave(plt_5, file = "figure_5.pdf", width = 20, height = 16)

############################################################################
# figure 2a

load("data/pred_without_bootstrap_without_drug_sim.RData")

test_info$pred = pred
test_info$obs = test_info$SYNERGY_SCORE

df = test_info[test_info$cohort == "test", ]

taji = read.csv("data/pred_ori1.csv")
taji$exp = paste0(taji$CELL_LINE, ".", taji$COMBINATION_ID)
taji$exp = toupper(taji$exp)
pick = match(df$experiment, taji$exp)

df$taji_pred = taji[pick, ]$PREDICTION
c1 = cor(df$pred, df$obs, method = "spearman")
c2 = cor(df$taji_pred, df$obs, method = "spearman", use="complete.obs")

plt_2a <- ggplot(df, aes(
    x = pred,
    y = obs)) +
    geom_point(shape = 21, alpha = 0.85, fill = "#F0E2B6", size = 4, colour = "#966729", stroke = 2) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Predicted", y = "Observed") +
    theme(plot.title = element_text(size = 26)) +
    theme(axis.title = element_text(size = 26)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    geom_abline(slope = 1, intercept = 0) +
    annotate(geom="text", x = 30, y = -120, label=parse(text="r == 0.50(p-val < 2.2 %*% 10^-16)"), size = 8) +
    ggtitle("Test Set 1 - DIPx")
plt_2a

##### for bootstrap #############################
# figure 2b
taji_bpred = read.csv("data/guan_bootstrap_pred.csv")
rownames(taji_bpred) = toupper(paste0(taji_bpred$CELL_LINE, ".", taji_bpred$COMBINATION_ID))
bpred = read.csv("data/bootstrap_pred.csv")
rownames(bpred) = toupper(paste0(bpred$CELL_LINE, ".", bpred$COMBINATION_ID))

taji_bpred = taji_bpred[rownames(df), ]
bpred = bpred[rownames(df), ]

cor_per_bootstrap = sapply(1:100, function(i){
  x = bpred[, paste0("bpred_", i)] 
  cor(x, df$obs, method = "spearman")
})

taji_cor_per_bootstrap = sapply(1:100, function(i){
  x = taji_bpred[, paste0("bpred_", i)] 
  cor(x, df$obs, method = "spearman", use="complete.obs")
})
d_cor = cor_per_bootstrap - taji_cor_per_bootstrap
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

cors = data.frame(cor_per_bootstrap, taji_cor_per_bootstrap)

df = melt(cors)
colnames(df) = c("group", "cor")
df$group = as.factor(df$group)

# df = df[-which(df$synergy > 230), ]
t.test(cors[[1]], cors[[2]])

t_test_df <- data.frame(
  group1 = 1, group2 = 2,
  label = c("p-val < 0.02"),
  y.position = c(0.485)
)

library(ggsci)
colors = pal_nejm(alpha = 0.9)(7)
colors = colors[c(3, 4)]

plt_2b <- ggplot(df) + 
  geom_boxplot(aes(x = group, y = cor, fill = group), trim = FALSE,  outlier.shape = NA) +
  add_pvalue(t_test_df, label = "label", tip.length = 0, label.size = 8, parse = TRUE)+
  theme_classic() +
  scale_fill_brewer(palette="Dark2") +
  # scale_fill_manual(values = colors) +
  theme(legend.position = "none") +
  labs(y = "Boostrap - Spearman cor(pred,obs)", x = "") +
  scale_x_discrete(labels = list("DIPx", "TAIJI-M")) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  theme(plot.title= element_text(size = 26)) +
  ylim(c(0.35, 0.5)) 
  # ggtitle("")
plt_2b



#########################################################
# figure 2c
load("data/pred_without_bootstrap_without_drug_sim.RData")

test_info$pred = pred
test_info$obs = test_info$SYNERGY_SCORE

df = test_info[test_info$cohort == "test", ]
taji_bpred = read.csv("data/guan_bootstrap_pred.csv")
rownames(taji_bpred) = toupper(paste0(taji_bpred$CELL_LINE, ".", taji_bpred$COMBINATION_ID))
bpred = read.csv("data/bootstrap_pred.csv")
rownames(bpred) = toupper(paste0(bpred$CELL_LINE, ".", bpred$COMBINATION_ID))

taji_bpred = taji_bpred[rownames(df), ]
bpred = bpred[rownames(df), ]


# cls =  unique(df$CELL_LINE)
cls = table(df$CELL_LINE)
cls = names(cls[cls >= 7])



cor_per_cl = sapply(cls, function(cl) {
  x = df[df$CELL_LINE ==  cl, ]
  bx = bpred[rownames(x), ]
  median(sapply(1:100, function(i) {
    cor(bx[, paste0("bpred_", i)] ,x$obs, method = "spearman")

  }))
})

taji_cor_per_cl = sapply(cls, function(cl) {
  x = df[df$CELL_LINE ==  cl, ]
  bx = taji_bpred[rownames(x), ]
  median(sapply(1:100, function(i) {
    cor(bx[, paste0("bpred_", i)] ,x$obs, method = "spearman")

  }))
})

df = data.frame(x = cor_per_cl, y = taji_cor_per_cl)
pick = which(df$x > df$y)
df$color = 0
df[pick, ]$color = 1

library(ggrepel)
library(RColorBrewer)
library(ggpubr)

colors <- brewer.pal(9, "Set1")

plt_2c <- ggplot(df, aes(
    x = x,
    y = y)) +
    geom_point(alpha = 0.8, size = 6, aes(col = factor(color))) +
    scale_colour_manual(values = colors[c(9, 5)]) +
    # geom_point(aes(
      # color = as.factor(color),
      # shape = as.factor(shape), size = size
    # ), alpha = 0.3) +
    # scale_size_identity()+
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "DIPx", y = "TAIJI-M") +
    theme(plot.title = element_text(size = 26)) +
    theme(axis.title = element_text(size = 26)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    geom_abline(slope = 1, intercept = 0) +
    ggtitle("cor(pred, obs) per cell line") +
    xlim(c(-0.55, 1)) + 
    ylim(c(-0.55, 1)) +
    annotate(geom="text", x=0.75, y=-0.25, label="63%", size = 26/.pt) + 
    annotate(geom="text", x=-0.2, y=0.75, label="37%", size = 26/.pt)

plt_2abc <- ggarrange(plt_2a, plt_2b, plt_2c,
  nrow = 1, ncol = 3, common.legend = F,
  labels = c("a", "b", "c") , font.label = list(size = 30)
)
plt_2abc

########################################################################################################
#figure 2d

load("data/pred_without_bootstrap_without_drug_sim.RData")

test_info$pred = pred
test_info$obs = test_info$SYNERGY_SCORE

df = test_info[test_info$cohort == "ext", ]

taji_bpred = read.csv("data/guan_bootstrap_pred.csv")
rownames(taji_bpred) = toupper(paste0(taji_bpred$CELL_LINE, ".", taji_bpred$COMBINATION_ID))
bpred = read.csv("data/bootstrap_pred.csv")
rownames(bpred) = toupper(paste0(bpred$CELL_LINE, ".", bpred$COMBINATION_ID))

taji_bpred = taji_bpred[rownames(df), ]
bpred = bpred[rownames(df), ]

cor_per_bootstrap = sapply(1:100, function(i){
  x = bpred[, paste0("bpred_", i)] 
  cor(x, df$obs, method = "spearman")
})

taji_cor_per_bootstrap = sapply(1:100, function(i){
  x = taji_bpred[, paste0("bpred_", i)] 
  cor(x, df$obs, method = "spearman", use="complete.obs")
})

d_cor = cor_per_bootstrap - taji_cor_per_bootstrap
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

cors = data.frame(cor_per_bootstrap, taji_cor_per_bootstrap)

df = melt(cors)
colnames(df) = c("group", "cor")
df$group = as.factor(df$group)

t.test(cors[[1]], cors[[2]])

t_test_df <- data.frame(
  group1 = 1, group2 = 2,
  label = c("0.07"),
  y.position = c(0.263)
)

plt_2d <- ggplot(df) + 
  geom_boxplot(aes(x = group, y = cor, fill = group), trim = FALSE,  outlier.shape = NA) +
  theme_classic() +
  scale_fill_brewer(palette="Dark2") +
  theme(legend.position = "none") +
  labs(y = "Boostrap - Spearman cor(pred, obs)", x = "") +
  scale_x_discrete(labels = list("DIPx", "TAIJI-M")) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  theme(plot.title= element_text(size = 26)) +
  add_pvalue(t_test_df, tip.length = 0, label.size = 8, parse = T)+
  ylim(c(0.15, 0.27)) +
  ggtitle("Test Set 2")
plt_2d

#################################
# figure 2e
load("data/pred_without_bootstrap_without_drug_sim.RData")

test_info$pred = pred
test_info$obs = test_info$SYNERGY_SCORE

df = test_info[test_info$cohort == "ext", ]

taji_bpred = read.csv("data/guan_bootstrap_pred.csv")
rownames(taji_bpred) = toupper(paste0(taji_bpred$CELL_LINE, ".", taji_bpred$COMBINATION_ID))
bpred = read.csv("data/bootstrap_pred.csv")
rownames(bpred) = toupper(paste0(bpred$CELL_LINE, ".", bpred$COMBINATION_ID))

taji_bpred = taji_bpred[rownames(df), ]
bpred = bpred[rownames(df), ]


load("data/use.drug_challenge.RData")
use.drug = use.drug[which(use.drug$cohort == "train"), ]
train_drugs = unique(c(use.drug$COMPOUND_A, use.drug$COMPOUND_B))
length(train_drugs)

pick2 = which(df$COMPOUND_A %in% train_drugs & df$COMPOUND_B %in% train_drugs)
length(pick2)

pick11 = which(df$COMPOUND_A %in% train_drugs & !(df$COMPOUND_B %in% train_drugs))
pick12 = which(!(df$COMPOUND_A %in% train_drugs) & (df$COMPOUND_B %in% train_drugs))
length(pick11)
length(pick12)
pick1 = c(pick11, pick12)

pick0 = which(!(df$COMPOUND_A %in% train_drugs) & !(df$COMPOUND_B %in% train_drugs))
length(pick0)

# 2 drugs
cor_per_b_2 = sapply(1:100, function(i) {
  y = df[pick2, ]$obs
  x = bpred[pick2, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})

taji_cor_per_b_2 = sapply(1:100, function(i) {
  y = df[pick2, ]$obs
  x = taji_bpred[pick2, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})

d_cor = cor_per_b_2 - taji_cor_per_b_2
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)


# 1 drug
cor_per_b_1 = sapply(1:100, function(i) {
  y = df[pick1, ]$obs
  x = bpred[pick1, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})

taji_cor_per_b_1 = sapply(1:100, function(i) {
  y = df[pick1, ]$obs
  x = taji_bpred[pick1, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})

d_cor = cor_per_b_1 - taji_cor_per_b_1
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

# 0 drug
cor_per_b_0 = sapply(1:100, function(i) {
  y = df[pick0, ]$obs
  x = bpred[pick0, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})

taji_cor_per_b_0 = sapply(1:100, function(i) {
  y = df[pick0, ]$obs
  x = taji_bpred[pick0, paste0("bpred_", i)]
  cor(x, y, method = "spearman",  use="complete.obs")
})
d_cor = cor_per_b_0 - taji_cor_per_b_0
d_cor = -d_cor
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

cors = list(cor_per_b_0, taji_cor_per_b_0, cor_per_b_1, taji_cor_per_b_1, cor_per_b_2, taji_cor_per_b_2)


df = melt(cors)
colnames(df) = c("cor", "group")

df$color = 1
df[which(df$group %in% c(1, 3, 5)), ]$color = 0
df[which(df$group %in% c(2, 4, 6)), ]$group = df[which(df$group %in% c(2, 4, 6)), ]$group - 1

library(rstatix)
t_test_df <- df %>% group_by(group) %>% t_test(cor ~ color) %>% add_xy_position(x = "group", dodge = 0.8)
colors <- brewer.pal(3, "Dark2")

t_test_df$p = c("0.02", "0.25", "6 %*% 10^-4")


plt_2e <- ggplot(df) + 
  geom_boxplot(aes(x = factor(group), y = cor, fill = factor(color)), trim = FALSE,  outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(name = "Method", labels=c('DIPx', 'TAIJI-M'), values = colors[1:2]) +
  theme(legend.position = "none") +
  labs(y = "Boostrap - Spearman cor(pred, obs)", x = "In the training set") +
  scale_x_discrete(labels = c("No drug", "One drug", "Two drugs")) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  theme(plot.title= element_text(size = 26)) +
  add_pvalue(t_test_df, tip.length = 0, label.size = 8, xmin = "xmin", 
               xmax = "xmax", parse = T)+
  ylim(c(0.05, 0.37)) +
  theme(legend.position = "none")+
  ggtitle("Test Set 2")
plt_2e


plt_2de <- ggarrange(plt_2d, plt_2e + rremove("ylab"),
  nrow = 1, ncol = 2, common.legend = F,
  labels = c("d", "e") , font.label = list(size = 30)
)
plt_2de

##############################################################################################
# figure 2fg
load("data/pred_without_bootstrap_without_drug_sim.RData")

test_info$pred = pred
test_info$obs = test_info$SYNERGY_SCORE

ic50_threshold = readRDS("data/IC50_threshold")


pick = match(test_info$COMPOUND_A, names(ic50_threshold))
test_info$IC50_A_threshold = ic50_threshold[pick]

pick = match(test_info$COMPOUND_B, names(ic50_threshold))
test_info$IC50_B_threshold = ic50_threshold[pick]

# testing set 2
plts = list()
set = c("test", "ext")
for(k in 1:2) {
  df = test_info[test_info$cohort == set[k], ]

  taji_bpred = read.csv("data/guan_bootstrap_pred.csv")
  rownames(taji_bpred) = toupper(paste0(taji_bpred$CELL_LINE, ".", taji_bpred$COMBINATION_ID))
  bpred = read.csv("data/bootstrap_pred.csv")
  rownames(bpred) = toupper(paste0(bpred$CELL_LINE, ".", bpred$COMBINATION_ID))

  taji_bpred = taji_bpred[rownames(df), ]
  bpred = bpred[rownames(df), ]

  pick2 = df$IC50_A <= df$IC50_A_threshold & df$IC50_B <= df$IC50_B_threshold
  pick0 = df$IC50_A > df$IC50_A_threshold & df$IC50_B > df$IC50_B_threshold

  pick1 = ((df$IC50_A <= df$IC50_A_threshold) &
          (df$IC50_B > df$IC50_B_threshold)) |
          ((df$IC50_A > df$IC50_A_threshold) &
          (df$IC50_B <= df$IC50_B_threshold))

  sum(pick2)
  sum(pick1)
  sum(pick0)

  # 2 drugs work
  cor_per_b_2 = sapply(1:100, function(i) {
    y = df[pick2, ]$obs
    x = bpred[pick2, paste0("bpred_", i)]
    cor(x, y, method = "spearman", use = "complete.obs")
  })

  taji_cor_per_b_2 = sapply(1:100, function(i) {
    y = df[pick2, ]$obs
    x = taji_bpred[pick2, paste0("bpred_", i)]
    cor(x, y, method = "spearman", use = "complete.obs")
  })
  d_cor = cor_per_b_2 - taji_cor_per_b_2
  z = (0-mean(d_cor))/sd(d_cor)
  2*pnorm(z)



  # 1 drug
  cor_per_b_1 = sapply(1:100, function(i) {
    y = df[pick1, ]$obs
    x = bpred[pick1, paste0("bpred_", i)]
    cor(x, y, method = "spearman", use = "complete.obs")
  })

  taji_cor_per_b_1 = sapply(1:100, function(i) {
    y = df[pick1, ]$obs
    x = taji_bpred[pick1, paste0("bpred_", i)]
    cor(x, y, method = "spearman", use = "complete.obs")
  })
  d_cor = cor_per_b_1 - taji_cor_per_b_1
  z = (0-mean(d_cor))/sd(d_cor)
  2*pnorm(z)

  # 0 drug
  cor_per_b_0 = sapply(1:100, function(i) {
    y = df[pick0, ]$obs
    x = bpred[pick0, paste0("bpred_", i)]
    cor(x, y, method = "spearman", use = "complete.obs")
  })

  taji_cor_per_b_0 = sapply(1:100, function(i) {
    y = df[pick0, ]$obs
    x = taji_bpred[pick0, paste0("bpred_", i)]
    cor(x, y, method = "spearman",  use="complete.obs")
  })

  d_cor = cor_per_b_0 - taji_cor_per_b_0
  z = (0-mean(d_cor))/sd(d_cor)
  2*pnorm(z)

  cors = list(cor_per_b_0, taji_cor_per_b_0, cor_per_b_1, taji_cor_per_b_1, cor_per_b_2, taji_cor_per_b_2)


  df = melt(cors)
  colnames(df) = c("cor", "group")
  # df$group = as.factor(df$group)

  df$color = 1
  df[which(df$group %in% c(1, 3, 5)), ]$color = 0
  df[which(df$group %in% c(2, 4, 6)), ]$group = df[which(df$group %in% c(2, 4, 6)), ]$group - 1

  t_test_df <- df %>% group_by(group) %>% t_test(cor ~ color) %>% add_xy_position(x = "group", dodge = 0.8)
  if(k == 1) {
    t_test_df$p = c("0.09", "0.07", "0.20")
  } else {
    t_test_df$p = c("0.68", "0.15", "1.4 %*% 10^-4")
  }

  colors <- brewer.pal(3, "Dark2")


  plt <- ggplot(df) + 
    geom_boxplot(aes(x = factor(group), y = cor, fill = factor(color)), trim = FALSE,  outlier.shape = NA) +
    theme_classic() +
    scale_fill_manual(name = "Method", labels=c('DIPx', 'TAIJI-M'), values = colors[1:2]) +
    theme(legend.position = "none") +
    labs(y = "Boostrap - Spearman cor(pred, obs)", x = "Monotherapy sensitivity") +
    scale_x_discrete(labels = c("No drug", "One drug", "Two drugs")) +
    theme(axis.text.y = element_text(size = 22)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.title = element_text(size = 26)) +
    theme(plot.title= element_text(size = 26)) +
    add_pvalue(t_test_df, tip.length = 0, label.size = 8, xmin = "xmin", 
               xmax = "xmax", parse = T)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    ggtitle(paste0("Test Set ", k))

  plt <- plt + theme(legend.position = "none")
  if(k == 2) {
    plt <- plt +
    theme(legend.position = c(0.6, 0.05)) +
    theme(legend.title = element_text(size = 26),
            legend.text = element_text(size = 26)) +
    theme(legend.direction="horizontal",
          legend.background=element_blank()) +
    theme(legend.key.size = unit(1.5, 'cm'))
  }
  plts[[k]] = plt
}


plt_2fg <- ggarrange(plts[[1]], plts[[2]] + rremove("ylab"),
  nrow = 1, ncol = 2, common.legend = F,
  labels = c("f", "g"), font.label = list(size = 30)
)
plt_2fg

plt_2 <- ggarrange(plt_2abc, plt_2de, plt_2fg,
  nrow = 3, ncol = 1, common.legend = T,
  labels = NULL , font.label = list(size = 30)
)
ggsave(plt_2, file = "figure_2.png", width = 20, height = 22)


########################################################
# figure 3b

load("data/result_trainquant_merged.RData")
bpred = merckPred
rm(merckPred)

rownames(bpred) = test_info$Folder
colnames(bpred) = paste0("bpred_", 1:100)

pred = readRDS("data/pred")

test_info$obs = test_info$SYNERGY_SCORE
test_info$pred = pred


cor_per_b = sapply(1:100, function(i) {
  y = test_info$obs
  x = bpred[, paste0("bpred_", i)]
  cor(x, y, method = "spearman", use = "complete.obs")
})

merck_cl = read.csv("data/merck_cellline.csv", sep = "\t")
pick = match(test_info$CELL_LINE, merck_cl$cell_line)
test_info$tissue = merck_cl[pick, ]$tissue

cor_per_cl = lapply(unique(test_info$CELL_LINE), function(cl) {
  pick = which(test_info$CELL_LINE == cl)
  x = bpred[pick, ]
  y = test_info[pick, ]$obs
  sapply(1:100, function(i) {
    cor(x[, paste0("bpred_", i)], y, method = "spearman")
  })
})
names(cor_per_cl) = unique(test_info$CELL_LINE)

median_cor_per_cl = sapply(cor_per_cl, median)
pick = order(median_cor_per_cl, decreasing = T)
cor_per_cl = cor_per_cl[pick]


df = melt(cor_per_cl)

pick = match(df[, 2], merck_cl$cell_line)
df$tissue = merck_cl[pick, ]$tissue

colnames(df) = c("cor", "cl", "tissue")
df$text = substr(df$tissue, 1, 1)

colors <- brewer.pal(6, "Set3")

df$tissue = paste0(df$tissue, " (", df$text, ")")

plt_3c <- ggplot(df) + 
  geom_boxplot(aes(x = reorder(factor(cl), cor, fun = median), y = cor, fill = factor(tissue)),  outlier.shape = NA) +
  geom_text(data = data.frame(), aes(x = df$cl, y=0.3, label=df$text), size=6,
            check_overlap = T) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  theme(legend.position = "none") +
  labs(y = "Boostrap - Spearman cor(pred, obs)", x = "Cell line") +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1)) +
  theme(axis.title = element_text(size = 26)) +
  theme(plot.title= element_text(size = 26)) +
  ylim(c(-0.2, 0.3)) +
  theme(legend.position = c(0.8, 0.1)) +
  theme(legend.title = element_blank(),
          legend.text = element_text(size = 20)) +
  theme(legend.direction="horizontal",
        legend.background=element_blank()) +
  theme(legend.key.size = unit(2, 'cm')) +
  # ggtitle("ONEIL") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")
plt_3c

#############################################
# figure 3a

ic50_threshold = readRDS("data/merck_IC50_threshold")


pick = match(test_info$COMPOUND_A, names(ic50_threshold))
test_info$IC50_A_threshold = ic50_threshold[pick]

pick = match(test_info$COMPOUND_B, names(ic50_threshold))
test_info$IC50_B_threshold = ic50_threshold[pick]

# testing set 2
df = test_info

pick2 = df$IC50_A <= df$IC50_A_threshold & df$IC50_B <= df$IC50_B_threshold
pick0 = df$IC50_A > df$IC50_A_threshold & df$IC50_B > df$IC50_B_threshold

pick1 = ((df$IC50_A <= df$IC50_A_threshold) &
        (df$IC50_B > df$IC50_B_threshold)) |
        ((df$IC50_A > df$IC50_A_threshold) &
        (df$IC50_B <= df$IC50_B_threshold))

sum(pick2)
sum(pick1)
sum(pick0)

# 2 drugs work
cor_per_b_2 = sapply(1:100, function(i) {
  y = df[pick2, ]$obs
  x = bpred[pick2, paste0("bpred_", i)]
  cor(x, y, method = "spearman", use = "complete.obs")
})

# 1 drug
cor_per_b_1 = sapply(1:100, function(i) {
  y = df[pick1, ]$obs
  x = bpred[pick1, paste0("bpred_", i)]
  cor(x, y, method = "spearman", use = "complete.obs")
})

# 0 drug
cor_per_b_0 = sapply(1:100, function(i) {
  y = df[pick0, ]$obs
  x = bpred[pick0, paste0("bpred_", i)]
  cor(x, y, method = "spearman", use = "complete.obs")
})

d_cor = cor_per_b_1 - cor_per_b_0
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

d_cor = cor_per_b_2 - cor_per_b_0
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

d_cor = cor_per_b_2 - cor_per_b_1
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

cors = list(cor_per_b_0, cor_per_b_1, cor_per_b_2)


df = melt(cors)
colnames(df) = c("cor", "group")


t_test_df <- df %>% t_test(cor ~ group) %>% add_xy_position(x = "group", dodge = 2)
t_test_df$y.position[2] = t_test_df$y.position[1] + 0.02
t_test_df$y.position[3] = t_test_df$y.position[2] + 0.02

colors <- brewer.pal(3, "Dark2")
t_test_df$p = c("0.54", "0.12", "0.04")


plt_3a <- ggplot(df) + 
  geom_boxplot(aes(x = factor(group), y = cor, fill = factor(group)),  outlier.shape = NA) +
  theme_classic() +
  scale_fill_brewer(palette="Greens") +
  theme(legend.position = "none") +
  labs(y = "Boostrap - Spearman cor(pred, obs)", x = "Monotherapy sensitivity") +
  scale_x_discrete(labels = c("No drug", "One drug", "Two drugs")) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  theme(plot.title= element_text(size = 26)) +
  add_pvalue(t_test_df, label = "p", tip.length = 0, label.size = 8, 
             xmin = "xmin", xmax = "xmax",
             parse = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "none")+
  ggtitle("ONEIL")
plt_3a


######################################################################################################
# supp figure 3
load("data/DreamChallenge_drugInfo_final_refined_09Dec2021.RData")

merckDrugs=unique(c(unique(test_info$COMPOUND_A),unique(test_info$COMPOUND_B)))
merckDrugs=toupper(merckDrugs)


drugInfo=drugInfo[drugInfo$drugNameUpdate!="",]
syno=drugInfo$drugNameUpdate
names(syno)=drugInfo$ChallengeName
comDrug_dream_merck=sapply(syno,function(x){
    x1=strsplit(x,",")[[1]]
    x1=trimws(x1)
    x1=toupper(x1)
    x1=gsub("-","",x1)
    x1=gsub(" ","",x1)
    x1=x1[x1!=""]
    x1=x1[x1 %in% merckDrugs]
    return(x1)
  })
pick=lengths(comDrug_dream_merck)>0
table(pick)

train_drugs = unlist(comDrug_dream_merck[pick])

test_info$COMPOUND_A = toupper(test_info$COMPOUND_A)
test_info$COMPOUND_B = toupper(test_info$COMPOUND_B)

df = test_info

pick2 = which(df$COMPOUND_A %in% train_drugs & df$COMPOUND_B %in% train_drugs)
length(pick2)

pick11 = which(df$COMPOUND_A %in% train_drugs & !(df$COMPOUND_B %in% train_drugs))
pick12 = which(!(df$COMPOUND_A %in% train_drugs) & (df$COMPOUND_B %in% train_drugs))
length(pick11)
length(pick12)
pick1 = c(pick11, pick12)

pick0 = which(!(df$COMPOUND_A %in% train_drugs) & !(df$COMPOUND_B %in% train_drugs))
length(pick0)

# 2 drugs
cor_per_b_2 = sapply(1:100, function(i) {
  y = df[pick2, ]$obs
  x = bpred[pick2, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})


# 1 drug
cor_per_b_1 = sapply(1:100, function(i) {
  y = df[pick1, ]$obs
  x = bpred[pick1, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})


# 0 drug
cor_per_b_0 = sapply(1:100, function(i) {
  y = df[pick0, ]$obs
  x = bpred[pick0, paste0("bpred_", i)]
  cor(x, y, method = "spearman")
})
d_cor = cor_per_b_1 - cor_per_b_0
d_cor = -d_cor
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

d_cor = cor_per_b_2 - cor_per_b_0
d_cor = -d_cor
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

d_cor = cor_per_b_2 - cor_per_b_1
d_cor = -d_cor
z = (0-mean(d_cor))/sd(d_cor)
2*pnorm(z)

cors = list(cor_per_b_0, cor_per_b_1, cor_per_b_2)

df = melt(cors)
colnames(df) = c("cor", "group")

t_test_df <- df %>% t_test(cor ~ group) %>% add_xy_position(x = "group", dodge = 2)
t_test_df$y.position[2] = t_test_df$y.position[1] + 0.02
t_test_df$y.position[3] = t_test_df$y.position[2] + 0.02

t_test_df$p = c("0.93", "0.19", "0.09")

plt_3b <- ggplot(df) + 
  geom_boxplot(aes(x = factor(group), y = cor, fill = factor(group)),  outlier.shape = NA) +
  theme_classic() +
  scale_fill_brewer(palette="Greens") +
  theme(legend.position = "none") +
  labs(y = "Boostrap - Spearman cor(pred, obs)", x = "In the training set") +
  # scale_x_continuous(breaks = c(1.5, 3.5), labels = c("PAMSYN", "TAJI-M")) +
  scale_x_discrete(labels = c("No drug", "One drug", "Two drugs")) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  theme(plot.title= element_text(size = 26)) +
  add_pvalue(t_test_df, label = "p", tip.length = 0, label.size = 8,
             xmin = "xmin", xmax = "xmax", parse = T)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ylim(c(-0.1, 0.25)) +
  theme(legend.position = "none")
plt_3b

plt_3 <- ggarrange(plt_3a, plt_3c, widths = c(1, 3),
                   nrow = 1, ncol = 2, common.legend = F,
  labels = c("a", "b") , font.label = list(size = 30)
)
plt_3
ggsave(plt_3, file = "figure_3.pdf", width = 22, height = 14)
ggsave(plt_3b, file = "figure_3b.pdf", width = 14, height = 14)

#################
# pas vs imp
load("data/AKT_1.ERBB_AKT.ERBB_imp.RData")
load("data/AKT.ERBB_breast_imp.RData")

quant_imp=quantile(imp,probs=0.95)
quant_imp

quant_pas =quantile(median_pas, probs=0.95)
quant_pas

plot(imp, median_pas, log="y", xlab="Feature importance", ylab="PAS expression (median)", cex.lab=1.5)
abline(v= quant_imp, col="blue")
abline(h= quant_pas, col="red")

pick= imp > quant_imp & median_pas > quant_pas

pick_imp = imp[pick]
pick_pas = median_pas[pick]

pick = order(pick_imp, decreasing = T)
pick_imp = pick_imp[pick]
pick_pas= pick_pas[pick]

pick_imp[1:10]


df = data.frame(imp = imp, pas = median_pas)
df$color = 1
pick = which(names(imp) %in% names(pick_imp)[1:6])
df[pick, ]$color = 2

plt <- ggplot(df, aes( x = imp,
    y = pas)) +
    geom_point(alpha = 0.6, size = 7, col = 'gray') +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL) +
    theme(plot.title = element_text(size = 26)) +
    theme(axis.title = element_text(size = 26)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    geom_hline(yintercept = quant_pas) +
    geom_vline(xintercept = quant_imp) +
    xlim(c(-2.5, 4)) + ylim(c(0, 13))
plt

ggsave(plt, file="AKT_ERBB.png", height = 10, width = 10, dpi = 300)

#######
# heatmap
library(reshape2)
library(ggplot2)
tbl = read.csv("data/top_pathways_BCL2L1_IAP.csv")
tbl = read.csv("data/top_pathways_Doxorubicin.TOP2.csv")
pathway = rownames(tbl)
pathway = sapply(pathway, function(s) {
  temp = unlist(strsplit(s, "_"))
  paste(temp[2:length(temp)], collapse = "_")
})
tbl$pathway = pathway
imp = tbl$Imp_score
tbl = tbl[, -1]

df = melt(tbl, id = c("pathway"))
table(df$variable)

df$value = round(df$value, 2)

plt <- ggplot(df, aes(variable, pathway)) + 
  geom_tile(aes(fill = value), color = "white", lwd = 0.75, linetype = 1) +
  geom_text(aes(label = value), size = 5) +
  scale_fill_gradient(low = "white", high = "red", limits = c(-10, 30))+
  ylim(rev(tbl$pathway)) +
  coord_fixed() +
  labs(title = "BCL2L1 + AZD5582, cell line: SW900", x = "Functional interaction", y = "Pathway") +
  # theme(axis.title=element_blank(),
  #       axis.text=element_blank(),
  #       axis.ticks=element_blank()) +
  guides(fill = guide_colourbar(title = "z-score", size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text=element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, size = 15, hjust = 1, vjust = 1)) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.position="top")
plt

ggsave(plt, file="BCL2L1_IAP_heatmap.png", height = 16, width = 16, dpi = 150)
ggsave(plt, file="Doxorubicin_TOP2_heatmap.png", height = 16, width = 16, dpi = 150)
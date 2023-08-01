library(doParallel)
library(randomForestSRC)
source("functions.R")

# number of cores for the pipeline
n_core <- 1
registerDoParallel(n_core)

### input ###
# data includes a toy example
target_updown_path <- "data/target_updown.RData"
synergy_path <- "data/synergy_score.RData"
gex_path <- "data/gex.RData"
geneset_path <- "data/geneset.RData"
driver_gene_path <- "data/driver_gene.RData"

### generate pas ###
pas_path <- generate_pas(gex_path, synergy_path, target_updown_path, geneset_path, driver_gene_path)


### load X, Y for train/test###
xy <- get_pas_matrix(pas_path, synergy_path, geneset_path)

# Here is the toy example so we use 90% for train and 10% for test
total_sample <- nrow(xy)
n_train <- round(total_sample * 0.9)

xy_train <- xy[1:n_train, ]
xy_test <- xy[(n_train + 1):total_sample, ]

model <- rfsrc(SYNERGY_SCORE ~ ., data = xy_train, importance = FALSE)

x_test <- xy_test[, -which(colnames(xy_test) == "SYNERGY_SCORE")]
pred <- predict.rfsrc(model, x_test)$predicted

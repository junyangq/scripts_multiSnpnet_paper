library(multiSnpnet)
library(snpnet)
library(glmnetPlus)
library(reshape2)
library(ggplot2)
library(data.table)
library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)

rank <- 20  # can be e.g. 5, 10, 20, 35
prev_iter <- 0  # -1 if starting from the latest iteration

genotype_file <- "#PATH/TO/GENOTYPE/DATA"   # ukb genotype data

phe_list_in <- "phe_list_biomarkers.csv"
phenotype_file <- "#PATH/TO/PHENOTYPE/DATA"   # path to the phenotype file
phe_list <- fread(phe_list_in, sep = ",")  # read list of phenotypes to analyze, two columns: (GBE_ID, phenotype)
phenotype_names <- phe_list[["GBE_ID"]]

covariate_names <- c("age", "sex", paste0("PC", 1:10))

standardize_response <- TRUE  # should we standardize response beforehand (rescale when prediction for sure)
save <- TRUE  # should the program save intermediate results?
weight <- NULL

# extend the lambda sequence where lambda_100/lambda_1 = 0.01 to length 150
nlambda <- 150
compute_lambda_min_ratio <- function(nlambda.new, nlambda = 100, ratio = 0.01) {
  exp((nlambda.new-1)/(nlambda-1)*log(ratio))
}
lambda.min.ratio <- compute_lambda_min_ratio(nlambda)

batch_size <- 4000  # size of candidate variant batch in the screening
max.iter <- 30  # maximum BAY iteration
thresh <- 1e-5  # convergence threshold
validation <- TRUE  # is validation set provided
early_stopping <- TRUE  # should we adopt early stopping

# other computational configurations
configs <- list(
  missing.rate = 0.2,  # variants above this missing rate are discarded
  MAF.thresh = 0.001,  # MAF threshold
  nCores = 16,  # number of cores to be used
  standardize.variant = FALSE,  # standardize predictors or not
  meta.dir = "meta/",
  results.dir = file.path(getwd(), paste0("results_rank_",rank)), # parent results directory to save intermediate results
  thresh = thresh,
  glmnet.thresh = thresh
)

mem <- 300000 # Should be specified if using a job scheduler

split.col <- "split"

fit <- multisnpnet(genotype_file = genotype_file,
                    phenotype_file = phenotype_file,
                    phenotype_names = phenotype_names,
                    covariate_names = covariate_names,
                    rank = rank,
                    batch_size = batch_size,
                    split_col = split.col,
                    max.iter = max.iter,
                    configs = configs,
                    lambda.min.ratio = lambda.min.ratio,
                    nlambda = nlambda,
                    standardize_response = standardize_response,
                    save = save,
                    validation = validation,
                    prev_iter = prev_iter,
                    early_stopping = early_stopping,
                    weight = weight,
                    mem = mem)

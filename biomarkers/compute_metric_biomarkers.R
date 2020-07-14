library(multiSnpnet)

fit = NULL

args <- commandArgs(trailingOnly = TRUE)

cat("Start computing for biomarkers rank", args[1], "\n")

saved_path = file.path("/oak/stanford/groups/mrivas/users/mrivas/repos/multiresponse-ukbb/results_biomarkers_unweighted/",
                       paste0("results_rank_", args[1]), "output_lambda_")  # Prefix of the saved R Data files

new_genotype_file = "/oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv"
new_phenotype_file = "/oak/stanford/groups/mrivas/projects/biomarkers/snpnet/biomarkers/biomarkers_covar.w_split.sorted.phe"
covariate_names = c("age", "sex", paste0("PC", 1:10))
split_col = "split"
split_name = c("train", "val", "test")
binary_phenotypes = NULL

phe_dir <- file.path("/oak/stanford/groups/mrivas/users/mrivas/repos/multiresponse-ukbb/results_biomarkers_unweighted/",
                     paste0("results_rank_", args[1]))

files_in_dir <- list.files(phe_dir)
result_files <- files_in_dir[startsWith(files_in_dir, "output_lambda_") & endsWith(files_in_dir, ".RData")]
max_iter <- max(as.numeric(gsub(".RData", "", gsub(pattern = "output_lambda_", "", result_files))))
idx <- 1:max_iter

save_to_dir <- file.path("/scratch/users/junyangq/multiresponse/code/multisnpnet/biomarkers/", paste0("results_rank_", args[1]))
dir.create(save_to_dir, recursive = T)

out <- predict_multisnpnet(fit = NULL, saved_path = saved_path, new_genotype_file = new_genotype_file, new_phenotype_file = new_phenotype_file,
                           idx = idx, covariate_names = covariate_names, split_col = split_col, split_name = split_name,
                           binary_phenotypes = binary_phenotypes, zstdcat_path = "zstdcat")

metric_train = t(out$R2[["train"]])
metric_val = t(out$R2[["val"]])
metric_test = t(out$R2[["test"]])
save(metric_train, metric_val, metric_test,
     file = file.path(save_to_dir, "output_lambda_1.RData"))
save(out, file = file.path(save_to_dir, "eval_results.RData"))

cat("End computing for biomarkers rank", args[1], "\n")

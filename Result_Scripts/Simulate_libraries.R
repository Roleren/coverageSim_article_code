#TODO: 
# 1. Analyse variance to justify bumps
# 2. Add MCC of simple features

library(coverageSim); library(ggplot2)


# Define output dirs
test_dir <- "/export/valenfs/projects/Hakon/coverageSim/"
data_output_dir <- file.path(test_dir, "data_output")
all_real_result_objects_path <- file.path(data_output_dir, "all_saved_objects.rdata")
# Setup for simulations
sensitivity_analysis_dir <- file.path(test_dir, "sensitivity_analysis")
all_sims <- c("no_features", "shape", "periodicity", "periodicity_true",
              "codon", "codonAndperiodicity", "default","optimized") 
all_sim_dirs <- file.path(sensitivity_analysis_dir, all_sims)
names(all_sim_dirs) <- all_sims
for (dir in all_sim_dirs) dir.create(dir, showWarnings = FALSE, recursive = TRUE)

# Load required data
if (file.exists(all_real_result_objects_path)) { 
  load(all_real_result_objects_path)
  source("/export/valenfs/projects/Hakon/coverageSim/Create_simulated_sample_data.R")
  merged_gr <- import.ofst(merged_gr_path, seqinfo = seqinfo(cds))
  merged_gr_covRle <- fimport(merged_gr_covRle_path)
} else stop("You must have run the real data validation and saved in to correct path!")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Sensitivity analysis
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Setup settings for simulations
sims_to_use <- all_sims; run_index <- seq_along(sims_to_use)
sim_dirs_to_use <- all_sim_dirs[sims_to_use]; 
exp_names <- paste0("SA_", sims_to_use); exp_obj_names <- paste0("df_", sims_to_use)
needs_forloop <- rep(FALSE, length(sims_to_use)); names(needs_forloop) <- sims_to_use
use_codons <- use_periodicity <- use_shape <-  needs_forloop

use_periodicity[c("periodicity" ,"periodicity_true", "codonAndperiodicity","default", 
                  "optimized")] <- c(-1,1,1,-1,1) 
use_codons[c("codon", "codonAndperiodicity",
             "optimized")] <- c(1,1,1)
needs_forloop <- (use_periodicity == 1) | (use_codons == 1) 
use_shape[c("shape", "default", "optimized")] <- c(9, 9, 9)
stopifnot(all(lengths(list(use_codons, use_periodicity,
               needs_forloop, use_shape, needs_forloop)) == length(sims_to_use)))
average_frame_use <- all_frame_usage[,mean(alpha), by = frame]$V1
frame4 <- average_frame_use[3]/3
average_frame_use <- c(frame4, average_frame_use[3:2],
                       average_frame_use, frame4)
res_region <- simCountTablesRegions(res, regionsToSample = "cds")
force_remake <- T
for (s in run_index[8]) {
  print(sims_to_use[s])
  exp_not_created <- !(exp_names[s] %in% all_exp$name)
  if (exp_not_created | force_remake) {
    out_dir <- sim_dirs_to_use[s]
    # Define settings
    ideal_cov <- list(cds = list(RFP = quote(rep.int(c(1, 0, 0), length.out = x))))
    frame_sub <-   if (use_periodicity[s] == -1) {average_frame_use} else NULL
    if (use_periodicity[s] == 0)  ideal_cov <- NULL
    tai_all <- NULL
    shape <- list(cds = list(RFP = shapes(use_shape[s])))
    
    print("Settings:")
    print(paste("-- Periodicity:", use_periodicity[s]))
    print(paste("-- Codon:", use_codons[s]))
    print(paste("-- Shape:", use_shape[s]))
    if (needs_forloop[s]) {
      files <- c()
      for (i in seq(1, ncol(res))) {
        if (use_periodicity[s] == 1) {
          frame_sub <- all_frame_usage[variable == unique(variable)[i],]$alpha
          frame44 <- frame_sub[3]/3
          frame_sub <- c(frame44, frame_sub[2:3],
                         frame_sub, frame44)
        }
        if (use_codons[s]) {
          tai_all <- AA_estimators[as.integer(variable) == i, c(2, 3), with = F]
        }
        files <- c(files, file.path(out_dir, paste0(colnames(assay(res[,i])), ".ofst")))
        simNGScoverage(genome_real, count_table = res_region[,i], 
                       out_dir = out_dir, 
                       exp_name = exp_names[s], 
                       ideal_coverage = ideal_cov,
                       seq_bias = tai_all, 
                       auto_correlation = shape,
                       rnase_bias = list(RFP = frame_sub),
                       validate = FALSE,
                       debug_coverage = F)
      }
      create.experiment(out_dir,
                        exp_names[s],
                        txdb = genome_real["txdb"],
                        fa = genome_real["genome"],
                        organism = "Homo sapiens", author = "Simulated by ORFikSim",
                        libtype = as.character(res$libtype),
                        condition = as.character(res$condition),
                        rep = res$replicate, files = files)
    } else {
      simNGScoverage(genome_real, count_table = res_region, 
                      out_dir = out_dir, 
                      exp_name = exp_names[s], 
                     ideal_coverage = ideal_cov,
                      seq_bias = tai_all, 
                      auto_correlation = shape, 
                      rnase_bias = list(RFP = frame_sub),
                      debug_coverage = F)
    }
  }
  assign(exp_obj_names[s], 
         read.experiment(exp_names[s])) #[c(1, 3:10, 2),]
}

# Create coverage table lists
# If only want to update single
dt_all_list_named <- readRDS(file.path(data_output_dir, "all_coverages.rds"))
X <- seq_along(exp_obj_names)[8]
exp_id <- toupper(substr(sims_to_use,1 , 1)); names(exp_id) <- sims_to_use
exp_id[c("periodicity_true", "codonAndperiodicity")] <- c("p", "T")
stopifnot(length(exp_id) == length(unique(exp_id)))
dt_all_list <- 
  lapply(X, function(x) coverage_all_cds_all_samples(get(exp_obj_names[x]), cds, exp_id[x]))
dt_all_list_named["optimized"] <- dt_all_list
saveRDS(dt_all_list_named, file.path(data_output_dir, "all_coverages.rds"))

# IF update all
X <- seq_along(exp_obj_names)[8]
exp_id <- toupper(substr(sims_to_use,1 , 1)); names(exp_id) <- sims_to_use
exp_id[c("periodicity_true", "codonAndperiodicity")] <- c("p", "T")
stopifnot(length(exp_id) == length(unique(exp_id)))
dt_all_list <- 
  lapply(X, function(x) coverage_all_cds_all_samples(get(exp_obj_names[x]), cds, exp_id[x]))
dt_all_list_named["optimized"] <- dt_all_list
saveRDS(dt_all_list_named, file.path(data_output_dir, "all_coverages.rds"))
dt_all_list_named <- dt_all_list_named[c("real", "optimized")]
dt_all_list_r <- c(list(dt_all_real), dt_all_list)
dt_all_list_named <- copy(dt_all_list_r)
names(dt_all_list_named) <- c("real", sims_to_use)
saveRDS(dt_all_list_named, file.path(data_output_dir, "all_coverages.rds"))
dt_all_list_named <- readRDS(file.path(data_output_dir, "all_coverages.rds"))

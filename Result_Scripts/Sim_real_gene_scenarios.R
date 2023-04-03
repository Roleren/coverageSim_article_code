# Simulate 4 example genes for the paper
library(coverageSim); library(ORFik); library(data.table); library(RiboCrypt); library(ggplot2)

# Define output dirs
test_dir <- "/export/valenfs/projects/Hakon/coverageSim/"
data_output_dir <- file.path(test_dir, "data_output")
all_real_result_objects_path <- file.path(data_output_dir, "all_saved_objects.rdata")
# Load required data
if (file.exists(all_real_result_objects_path)) { 
  load(all_real_result_objects_path)
  rm(list=lsf.str())
  genome_real <- c(gtf = ORFik:::getGtfPathFromTxdb(loadTxdb(df)),
                   genome = df@fafile, txdb = df@txdb)
  merged_gr <- import.ofst(merged_gr_path, seqinfo = seqinfo(cds))
  merged_gr_covRle <- fimport(merged_gr_covRle_path)
} else stop("You must have run the real data validation and saved in to correct path!")

# uORF output dirs
SA_genes <- file.path(test_dir, "sensitivity_analysis_genes")
all_sims_g <- c("ATF4", "ADAM10", "C11", "uORF_gene")
all_sim_dirs_g <- file.path(SA_genes, all_sims_g)
names(all_sim_dirs_g) <- all_sims_g
for (dir in all_sim_dirs_g) dir.create(dir, showWarnings = FALSE, recursive = TRUE)
exp_names_g <- paste0("SA_sim_", all_sims_g); exp_obj_names_g <- paste0("df_", all_sims_g)
uORF_genes <- c("ENST00000674920", "ENST00000260408", "ENST00000228136", "ENSTTEST10001"); names(uORF_genes) <- all_sims_g
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
## uORF scenario simulation
# 1. ATF4
# 2. ADAM10
# 3. C11orf58
# 4. Complex_sim
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# First make genome for Complex_sim (2 uORFs overlapping CDS)
sim_genome <- simGenome(n = 1, out_dir = all_sim_dirs_g["uORF_gene"], 
                        max_uorfs = 2,
                        uorf_max_length = 60,
                        uorfs_can_overlap_cds = 2,
                        uorfs_can_overlap = TRUE, debug_on = FALSE)
sim_genome <- c(txdb = "/export/valenfs/projects/Hakon/coverageSim/sensitivity_analysis_genes/uORF_gene/Homo_sapiens_dummy.gtf.db",
                gtf = "/export/valenfs/projects/Hakon/coverageSim/sensitivity_analysis_genes/uORF_gene/Homo_sapiens_dummy.gtf", 
                genome = "/export/valenfs/projects/Hakon/coverageSim//sensitivity_analysis_genes/uORF_gene/Homo_sapiens_dummy.fasta", 
                uorfs = "/export/valenfs/projects/Hakon/coverageSim//sensitivity_analysis_genes/uORF_gene/true_uORFs.rds" )
sim_uORFs <- readRDS(sim_genome["uorfs"]); sim_uORFs
count_table_sim <- simCountTables(loadRegion(sim_genome["txdb"], "cds"),
                                  interceptMean = 10, libtypes = c("RFP"), 
                                  replicates = 1, conditions = "WT")

# Get genes and uORFs to plot
upstream <- c(270, 75, 110, 50); downstream <- c(150, 100, 100, 200)
start_codons <- c(rep("ATG", 2), "GTG")
uORFs_to_use <- list(seq(3), 1, 2)
regionsToSample <- c("leader","cds", "uorf")
# shape_genes <- list(leader = list(RFP = 1),
#                     cds = list(RFP = shapes(3, v = 1)),
#                     uorf = list(RFP = shapes(3, v = 1)))
R_r <- list(list(leader =  list(RFP = 0.05),
                 cds =     list(RFP = 0.4),
                 trailer = list(RFP = 0),
                 uorf =    list(RFP = 0.55)), 
            list(leader =  list(RFP = 0.015),
                 cds =     list(RFP = 0.98),
                 trailer = list(RFP = 0),
                 uorf =    list(RFP = 0.005)),
            list(leader =  list(RFP = 0.05),
                 cds =     list(RFP = 0.75),
                 trailer = list(RFP = 0),
                 uorf =    list(RFP = 0.2)),
            list(leader =  list(RFP = 0.05),
                 cds =     list(RFP = 0.65),
                 trailer = list(RFP = 0),
                 uorf =    list(RFP = 0.3))) 
uorf_prop <- list(c(1, 0.8, 4), "uniform", "uniform", "uniform")
periodicity <- unlist(frame_usage(data.table(R_M = rowSums(dt_all_real)))[, .(alpha)], use.names = FALSE)
periodicity44 <- periodicity[3]/3
periodicity <- c(periodicity44, periodicity[2:3],
                 periodicity, periodicity44)
# tai_merged <- coverageSim:::AA_score(cds_all[uORF_genes[1:3]], merged_gr_covRle,
#                        seqs = translate_orf_seq(cds_all[uORF_genes[1:3]],
#                                                 df, as = "AA", start.as.hash = T,
#                                                 stopm1.as.amp = T, startp1.as.per = T,
#                                                 return.as.list = T))
tai_merged <- load_seq_bias(); tai_merged$variable <- NULL
plot_y_breaks <- c(500, 50, 50, 20)
plot_prop <- list(c(10,6), c(10,4), c(5,3), c(6,7))
x_labels <- c(rep("", length(all_sims_g)-1), "position")
up_cds_show <- c(0, 12, 0, 0); down_cds_show <-c(170, 120, 100, 199)
plot_breaks <- c("     ", "    ", "    ", "    ")
plots_genes <- list()
dt_all_list_named_g <- data.table()
make_new_experiments <- FALSE

for (s in seq_along(uORF_genes)) {
  print(paste("-- Importing required data..."))
  tx_name <- uORF_genes[s]
  is_human_genome <- tx_name != uORF_genes["uORF_gene"]
  is_sim_genome <- !is_human_genome
  genome_to_use <- if (is_sim_genome) {sim_genome} else genome_real
  txdb_to_use <- genome_to_use["txdb"]
  print(paste("Gene: ", names(tx_name)))
  print(paste("- Importing gene regions"))
  #assign(paste0(names(tx_name), "_tx_name","ENST00000260408"))
  mrna_gene <- loadRegion(txdb_to_use, "mrna", names.keep = tx_name)
  cds_gene <- loadRegion(txdb_to_use, "cds", names.keep = tx_name)
  window_show <- windowPerGroup(startSites(cds_gene, T, T, T), mrna_gene,
                                upstream[s], downstream[s])
  uorf_prop_this <- uorf_prop[[s]]
  names(uorf_prop_this) <- rep(tx_name, length(uorf_prop_this))
  
  if (is_sim_genome) { 
    uORFs_gene <- sim_uORFs
    count_table_gene <- count_table_sim
  } else { # human genome uORF
    uORFs_gene <- findMapORFs(window_show, txSeqsFromFa(window_show, df, T, T),
                              longestORF = F, startCodon = start_codons[s])[uORFs_to_use[[s]]]
    temp <- coverageSim:::overlap_props(tx_name, uORFs_gene, merged_gr)
    print(temp[[1]])
    count_table_gene <- temp[[2]]
    if (tx_name == uORF_genes["C11"]) assay(count_table_gene) <- round(assay(count_table_gene)/2, 0)
  }
  
  if (make_new_experiments) {
    count_table_gene_region <- 
      simCountTablesRegions(count_table_gene, regionsToSample = c("leader", "cds", "uorf"), 
                            region_proportion =  R_r[[s]])
    print(paste("-- Simulating..."))
    uORFs_gene_tx_names <- uORFs_gene
    names(uORFs_gene_tx_names) <- txNames(uORFs_gene_tx_names)
    uORFs_gene_tx_names@unlistData$names <- NULL
    periodicity_low <- periodicity * c(1, 1, 1, 3, 1, 1, 1)
    tai_merged_used <- copy(tai_merged)
    tai_merged_used[, alpha := alpha*5]
    if (s == 1) {
      # ATF4 has lower start bias than usual
      tai_merged_used[1, ]$alpha <- tai_merged_used[1, ]$alpha*2/3
    }
     
    simNGScoverage(genome_to_use, count_table_gene_region, 
                   exp_name = exp_names_g[s],
                   out_dir = all_sim_dirs_g[s], 
                   rnase_bias = list(RFP = periodicity_low),
                   seq_bias = tai_merged_used, 
                   true_uorf_ranges = uORFs_gene_tx_names,
                   uorf_prop_within_gene = uorf_prop_this)
  }
  
  print(paste("-- Coverage..."))
  df_exp <- read.experiment(exp_names_g[s])
  sim_reads <- fimport(filepath(df_exp[1,], "pshifted"))
  cov_sim <- coveragePerTiling(window_show, sim_reads, is.sorted = T, as.data.table = T, withFrames = T, fraction = "Sim")
  if (is_human_genome) { 
    cov_real_gene <- coveragePerTiling(window_show, merged_gr_covRle, is.sorted = T,
                                      as.data.table = T, withFrames = T, fraction = "Real")
    dt_all_list_named_g <- rbindlist(list(dt_all_list_named_g, 
                                          data.table(genes = rep.int(names(tx_name), nrow(cov_sim)),
                                                     real = cov_real_gene$count, sim = cov_sim$count)))
  } else cov_real_gene <- data.table()
  print(paste("-- Plotting..."))
  plot <- suppressMessages(coverageSim:::ORF_frame_plot_window(tx_window = window_show,
                                                 sim_reads = sim_reads,
                                                 cov_real = cov_real_gene, cov_sim = cov_sim, uORFs = uORFs_gene,
                                                 cds_region = startRegion(cds_gene, upstream = up_cds_show[s],
                                                                          downstream = down_cds_show[s]), 
                                                 y_breaks = plot_y_breaks[s], 
                                                 break_string = plot_breaks[s],
                                                 prop = plot_prop[[s]], x_lab = x_labels[s]))
  plot(plot)
  assign(paste0("plot_", s), plot)
  plots_genes <- c(plots_genes, plot)
}

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merged plot
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
all_scenarios_plot <- cowplot::plot_grid(plot_2, plot_3, plot_1, plot_4, ncol = 2, 
                                         labels = "AUTO", label_size = 9)
ggsave(file.path(test_dir, "ORF_scenarios.png"), all_scenarios_plot, width = 6, heigh = 5, dpi = 600)
plot(all_scenarios_plot)
# all_scenarios_plot <- cowplot::plot_grid(cowplot::plot_grid(plot_2, plot_3, ncol = 2, 
#                                                             labels = "AUTO", label_size = 9, 
#                                                             rel_widths = c(1,1)),
#                                          cowplot::plot_grid(plot_1, plot_4, ncol = 2, 
#                                                             labels = c("C", "D"), label_size = 9, 
#                                                             rel_widths = c(2,1)),
#                                          ncol = 1)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Regression analysis
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
old_var <- dt_all_list_named_g[, .(var = sd(real)/mean(real), var_sim = sd(sim)/mean(sim)),
                    by = genes]

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Ribocrypt browser analysis
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load done experiments
df_simATF <- read.experiment(exp_names_g[1], output.env = new.env())
df_simAdam <- read.experiment("sim_test_adam", output.env = new.env())
df_simC11 <- read.experiment("sim_test_C11orf58", output.env = new.env())
# Special for complex sim
df_C <- read.experiment("Complex_sim", output.env = new.env())
cds_C <- loadRegion(df_C, "cds")
mrna_C <- loadRegion(df_C, "mrna")
window_show_C <- windowPerGroup(startSites(cds_C, T, T, T), mrna_C, 50, 200)

RiboCrypt::multiOmicsPlot_ORFikExp(window_show, df = df[1,],
                                   reads = list(RFP = merged_gr_covRle), 
                                   annotation = uORFs, 
                                   frames_type = "columns",
                                   BPPARAM = BiocParallel::SerialParam())
RiboCrypt::multiOmicsPlot_ORFikExp(window_show_adam, df = df[1,],
                                   reads = list(RFP = merged_gr_covRle), 
                                   annotation = uORFs_adam, 
                                   frames_type = "columns",
                                   BPPARAM = BiocParallel::SerialParam())
RiboCrypt::multiOmicsPlot_ORFikExp(mrna_RPL, df = df[1,],
                                   reads = list(RFP = merged_gr_covRle), 
                                   annotation = cds_RPL, 
                                   frames_type = "columns",
                                   BPPARAM = BiocParallel::SerialParam())
RiboCrypt::multiOmicsPlot_ORFikExp("ENSTTEST10001", df = df_C,
                                   annotation = cds_C, 
                                   frames_type = "columns",
                                   custom_regions = sim_uORFs,
                                   BPPARAM = BiocParallel::SerialParam())
# OLD
# # ATF4
# ATF_tx_name <- "ENST00000674920"
# mrna_ATF4 <- loadRegion(df, "mrna", names.keep = ATF_tx_name)
# cds_ATF4 <- loadRegion(df, "cds", names.keep = ATF_tx_name)
# 
# RiboCrypt::multiOmicsPlot_ORFikExp(window_show, df = df[1,],
#                                    reads = list(RFP = merged_gr_covRle), 
#                                    annotation = uORFs, 
#                                    frames_type = "columns",
#                                    BPPARAM = BiocParallel::SerialParam())
# 
# window_show <- windowPerGroup(startSites(cds_ATF4, T, T, T), mrna_ATF4, 300, 150)
# uORFs <- findMapORFs(window_show, ORFik::txSeqsFromFa(window_show, df, T, T), longestORF = F, startCodon = "ATG")
# 
# if (!("sim_test_ATF4" %in% all_exp$name)) {
#   temp <- overlap_props(ATF_tx_name, uORFs, merged_gr)
#   temp[[1]]; count_table_ATF4 <- temp[[2]]
#   uorf_prop_ATF4 <- c(1, 0.8, 4); names(uorf_prop_ATF4) <- rep(ATF_tx_name, length(uORFs))
#   
#   tai_uORF <- copy(tai)
#   tai_uORF[trans == "A", usage_score := 30]
#   
#   df_simATF <- simNGScoverage(simGenome = genome_real, count_table = count_table_ATF4, exp_name = "sim_test_ATF4",
#                               regionsToSample = c("leader","cds", "uorf"), out_dir = cor_sim_dir_ATF4, 
#                               region_proportion = 
#                                 list(leader =  list(RFP = 0.05),
#                                      cds =     list(RFP = 0.4),
#                                      trailer = list(RFP = 0),
#                                      uorf =    list(RFP = 0.55)), 
#                               skew_covs = list(leader =  list(RFP = 1),
#                                                cds =     list(RFP = c(5, 2, 1)),
#                                                trailer = list(RFP = 1),
#                                                uorf = list(RFP = c(5, 2, 1))),
#                               tAI = tai_uORF, true_uorf_ranges = uORFs, 
#                               uorf_prop_within_gene = uorf_prop_ATF4)
# }
# 
# df_simATF <- read.experiment("sim_test_ATF4")
# envExp(df_simATF) <- new.env()
# 
# cov_real_ATF4 <- coveragePerTiling(window_show, merged_gr_covRle, is.sorted = T, as.data.table = T, withFrames = T, fraction = "Real")
# 
# # ADAM10
# ADAM10_tx_name <- "ENST00000260408"
# mrna_ADAM10 <- loadRegion(df, "mrna", names.keep = ADAM10_tx_name)
# cds_ADAM10 <- loadRegion(df, "cds", names.keep = ADAM10_tx_name)
# window_show_adam <- windowPerGroup(startSites(cds_ADAM10, T, T, T), mrna_ADAM10, 75, 100)
# uORFs_adam <- findMapORFs(window_show_adam, ORFik::txSeqsFromFa(window_show_adam, df, T, T), longestORF = F, startCodon = "ATG")[1]
# RiboCrypt::multiOmicsPlot_ORFikExp(window_show_adam, df = df[1,],
#                                    reads = list(RFP = merged_gr_covRle), 
#                                    annotation = uORFs_adam, 
#                                    frames_type = "columns",
#                                    BPPARAM = BiocParallel::SerialParam())
# temp_adam <- overlap_props(ADAM10_tx_name, uORFs_adam, merged_gr)
# temp_adam[[1]]; temp_adam[[2]]
# tai_uORF_adam <- copy(tai)
# tai_uORF_adam[trans == "C", usage_score := 10]; tai_uORF_adam[trans == "L", usage_score := 20]
# if (!("sim_test_adam" %in% all_exp$name)) {
#   df_simAdam <- simNGScoverage(simGenome = genome_real, count_table = temp_adam[[2]],
#                                exp_name = "sim_test_adam",
#                                regionsToSample = c("leader","cds", "uorf"), out_dir = cor_sim_dir_adam, 
#                                region_proportion = 
#                                  list(leader =  list(RFP = 0.015),
#                                       cds =     list(RFP = 0.98),
#                                       trailer = list(RFP = 0),
#                                       uorf =    list(RFP = 0.005)), 
#                                skew_covs = list(leader =  list(RFP = 1),
#                                                 cds =     list(RFP = c(5, 2, 1)),
#                                                 trailer = list(RFP = 1),
#                                                 uorf = list(RFP = c(5, 2, 1))),
#                                tAI = tai_uORF_adam, true_uorf_ranges = uORFs_adam)
# }
# df_simAdam <- read.experiment("sim_test_adam")
# cov_real_adam <- coveragePerTiling(window_show_adam, merged_gr, is.sorted = T, as.data.table = T, withFrames = T, fraction = "Real")
# 
# # C11orf58
# RPL_tx_name <- "ENST00000228136"
# mrna_RPL <- loadRegion(df, "mrna", names.keep = RPL_tx_name)
# cds_RPL <- loadRegion(df, "cds", names.keep = RPL_tx_name)
# window_show_RPL <- windowPerGroup(startSites(cds_RPL, T, T, T), mrna_RPL, 110, 100)
# uORFs_RPL <- findMapORFs(window_show_RPL, ORFik::txSeqsFromFa(window_show_RPL, df, T, T), longestORF = F, startCodon = "GTG")[2]
# RiboCrypt::multiOmicsPlot_ORFikExp(mrna_RPL, df = df[1,],
#                                    reads = list(RFP = merged_gr_covRle), 
#                                    annotation = cds_RPL, 
#                                    frames_type = "columns",
#                                    BPPARAM = BiocParallel::SerialParam())
# seqs_RPL <- txSeqsFromFa(uORFs_RPL, df, is.sorted = TRUE)
# seqs_RPL <- translate(seqs_RPL)
# subseq(seqs_RPL, 1, 1) <- "#"
# seqs_RPL <- strsplit(as.character(unlist(seqs_RPL, use.names = FALSE)), split = "")
# RPL_AA <- AA_score(cds_RPL, reads = merged_gr_covRle, seqs = seqs_RPL)
# temp_RPL <- overlap_props(RPL_tx_name, uORFs_RPL, merged_gr)
# temp_RPL[[1]]; temp_RPL[[2]]
# # Fix double CDS bias
# summarized_RPL <- temp_RPL[[2]]
# assay(summarized_RPL) <- round(assay(summarized_RPL)/2, 0)
# tai_uORF_RPL <- copy(tai)
# tai_uORF_RPL <- RPL_AA; tai_uORF_RPL[, trans := AA]; tai_uORF_RPL[, usage_score := AA_ratio_ccn]
# tai_uORF_RPL[chmatch(c("#","C", "L", "K", "A", "P"), tai_uORF_RPL$trans),  usage_score := c(80,10, 30, 5, 5, 30)]
# tai_uORF_RPL[trans != "#", usage_score := usage_score/1.2]
# if (!("sim_test_C11orf58" %in% all_exp$name)) {
#   df_simRPL <- simNGScoverage(simGenome = genome_real, count_table = summarized_RPL,
#                               exp_name = "sim_test_C11orf58",
#                               regionsToSample = c("leader","cds", "uorf"), out_dir = cor_sim_dir_C, 
#                               region_proportion = 
#                                 list(leader =  list(RFP = 0.05),
#                                      cds =     list(RFP = 0.75),
#                                      trailer = list(RFP = 0),
#                                      uorf =    list(RFP = 0.2)), 
#                               skew_covs = list(leader =  list(RFP = 1),
#                                                cds =     list(RFP = c(5, 2, 1)),
#                                                trailer = list(RFP = 1),
#                                                uorf = list(RFP = c(5, 2, 1))),
#                               tAI = tai_uORF_RPL, true_uorf_ranges = uORFs_RPL)
# }
# df_simRPL <- read.experiment("sim_test_C11orf58")
# C11orf58_plot <- ORF_frame_plot_window(tx_window = window_show_RPL,
#                                        sim_reads = fimport(filepath(df_simRPL[1,], "pshifted")),
#                                        cov_real = cov_real_RPL,
#                                        uORFs = uORFs_RPL, cds_region = windowPerGroup(startSites(cds_RPL, T, T, T), mrna_RPL, 0, 100), 
#                                        break_string = "    ", prop = c(5,3), x_lab = "")
# plot(C11orf58_plot)
# 
# 
# 
# # Simulated double overlapping
# if (!("Complex_sim" %in% all_exp$name)) {
#   sim_genome <- simGenome(n = 1, out_dir = cor_sim_dir_simulated, 
#                           max_uorfs = 2,
#                           uorf_max_length = 60,
#                           uorfs_can_overlap_cds = 2,
#                           uorfs_can_overlap = TRUE, debug_on = FALSE)
#   sim_uORFs <- readRDS(sim_genome["uorfs"]); sim_uORFs
#   count_table_sim <- simCountTables(loadRegion(sim_genome["txdb"], "cds"),
#                                     interceptMean = 10, libtypes = c("RFP"), 
#                                     replicates = 1, conditions = "WT")
#   df_C <- simNGScoverage(sim_genome, count_table_sim, exp_name = "Complex_sim", region_proportion = 
#                            list(leader =  list(RFP = 0.05),
#                                 cds =     list(RFP = 0.65),
#                                 trailer = list(RFP = 0),
#                                 uorf =    list(RFP = 0.3)), tAI = tai)
# }
# df_C <- read.experiment("Complex_sim")
# envExp(df_C) <- new.env()
# cds_C <- loadRegion(df_C, "cds")
# RiboCrypt::multiOmicsPlot_ORFikExp("ENSTTEST10001", df = df_C,
#                                    annotation = cds_C, 
#                                    frames_type = "columns",
#                                    custom_regions = sim_uORFs,
#                                    BPPARAM = BiocParallel::SerialParam())
# mrna_C <- loadRegion(df_C, "mrna")
# 
# # Plot all
# ADAM_plot <- ORF_frame_plot_window(tx_window = window_show_adam,
#                                    sim_reads = fimport(filepath(df_simAdam[1,], "pshifted")),
#                                    cov_real = cov_real_adam,
#                                    uORFs = uORFs_adam, cds_region = startRegion(cds_ADAM10, upstream = 15, downstream = 120), 
#                                    y_breaks = 50, break_string = "    ", prop = c(10,4), x_lab = "")
# plot(ADAM_plot)
# 
# C11orf58_plot <- ORF_frame_plot_window(tx_window = window_show_RPL,
#                                        sim_reads = fimport(filepath(df_simRPL[1,], "pshifted")),
#                                        cov_real = cov_real_RPL,
#                                        uORFs = uORFs_RPL, cds_region = windowPerGroup(startSites(cds_RPL, T, T, T), mrna_RPL, 0, 100), 
#                                        break_string = "    ", prop = c(5,3), x_lab = "")
# plot(C11orf58_plot)
# 
# ATF4_plot <- ORF_frame_plot_window(tx_window = window_show,
#                                    sim_reads = fimport(filepath(df_simATF[1,], "pshifted")),
#                                    cov_real = cov_real_ATF4,
#                                    uORFs = uORFs, cds_region = startRegion(cds_ATF4, upstream = -1, downstream = 170), 
#                                    y_breaks = 500, break_string = "      ", prop = c(10,6), x_lab = "")
# plot(ATF4_plot)
# 
# sim_double_plot <- ORF_frame_plot_window(tx_window = windowPerGroup(startSites(cds_C, T, T, T), mrna_C, 50, 200),
#                                          sim_reads = fimport(filepath(df_C[1,], "pshifted")),
#                                          cov_real = data.table(),
#                                          uORFs = sim_uORFs, cds_region = windowPerGroup(startSites(cds_C, T, T, T), mrna_C, 0, 199), 
#                                          y_breaks = 20, prop = c(6,7), x_lab = "position")
# plot(sim_double_plot)
# all_scenarios_plot <- gridExtra::arrangeGrob(ADAM_plot, C11orf58_plot, ATF4_plot, sim_double_plot,
#                                              layout_matrix = t(matrix(c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 5)), nrow = 1)))
# ggsave(file.path(test_dir, "ORF_scenarios.png"), all_scenarios_plot, width = 6, heigh = 7, dpi = 600)
# plot(all_scenarios_plot)

library(ORFik); library(data.table); library(RiboCrypt); library(ggplot2)
source("/export/valenfs/projects/Hakon/coverageSim/Create_simulated_sample_data.R")

## Define output paths
# Set this path as main directory ->
test_dir <- "/export/valenfs/projects/Hakon/coverageSim/"

# Output data files
data_output_dir <- file.path(test_dir, "data_output"); dir.create(data_output_dir, F, T)
all_real_result_objects_path <- file.path(data_output_dir, "all_saved_objects.rdata")
if (file.exists(all_real_result_objects_path)) { 
  load(all_real_result_objects_path)
} else {
merged_gr_path <- file.path(data_output_dir, "real_merged_gr.ofst")
merged_gr_covRle_path <- file.path(data_output_dir, "real_merged_gr.covrds")
gene_count_top_500_cds_path <- file.path(data_output_dir, "genecount_top_500_cds.rds")
study_info_path <- file.path(test_dir, "Study_info.csv")

# # Load predefined tai tables
# tai_table <- file.path(test_dir, "tAI_table/tAI_human.csv")
# tai_table_start <- file.path(test_dir, "tAI_table/tAI_human_with_start.csv")
# tai_table_start_m1s <- file.path(test_dir, "tAI_table/tAI_human_with_start_and_m1stop.csv")

# Read real data experiments
all_exp <- list.experiments(validate = FALSE)
b <- all_exp[organism == "Homo sapiens",]; b <- b[grep("RFP", (paste(libtypes, sep = ""))),]
b <- b[grep("PRJ", (paste(name, sep = ""))),]
study_names <- b$name[c(1, 3:11)] 
real_lib_names <- paste0("R", seq(length(study_names)))



# Pick high count subset of cds to use for simulation
df <- read.experiment("PRJEB26593")
genome_real <- c(gtf = ORFik:::getGtfPathFromTxdb(loadTxdb(df)),
                 genome = df@fafile, txdb = df@txdb)
cds_all <- loadRegion(df, "cds", names.keep = filterTranscripts(df))
cds_all_width <- widthPerGroup(cds_all, FALSE)
cds_all_is_mod3 <- (cds_all_width %% 3) == 0; sum(cds_all_is_mod3); sum(!cds_all_is_mod3)
reads_select <- fimport(filepath(df[1,], "pshifted"))
counts <- countOverlapsW(cds_all, reads_select, "score")
fpkms <- ORFik:::fpkm_calc(counts, widthPerGroup(cds_all,FALSE), sum(counts))
count_filter <- max(600, quantile(counts, 0.9)); count_filter
fpkm_filter <- max(10, quantile(fpkms, 0.9)); fpkm_filter
total_filter <- counts > count_filter & fpkms > fpkm_filter & cds_all_is_mod3
cds_candidates <- cds_all[total_filter];length(cds_candidates)
# Optional: For bad genomes validate start and stop are indeed correct valid codons.
# cds_all_is_start <- startRegionString(cds_candidates, cds_all, df, 0, 2) %in% 
#   c("ATG", "CTG", "TTG", "GTG")
# cds_all_is_stop <- txSeqsFromFa(ORFik:::removeMetaCols(stopCodons(cds_candidates, TRUE)),df, T, F) %in% 
#   c("TGA", "TAG", "TAA")
periodicity_filter <- orfScore(cds_candidates, reads_select, T, "score",
                               overlapGrl = counts[names(cds_candidates)])
cds_candidates <- cds_candidates[periodicity_filter$ORFScores > 0]
# Do not sort here by "best" to avoid only getting the extreme genes
cds <- cds_candidates[seq(min(500,length(cds_candidates)))]


leaders <- loadRegion(df, "leaders", names.keep = names(cds))
trailers <- loadRegion(df, "trailers", names.keep = names(cds))
mrna <- loadRegion(df, "mrna", names.keep = names(cds))
cds_with_A_site <- windowPerGroup(startSites(cds, T, T, T), mrna,
                                  downstream = widthPerGroup(cds,F) - 1, upstream = 3)
stopifnot(all((widthPerGroup(cds, F) + 3) == widthPerGroup(cds_with_A_site, F)))

seqs <- translate_orf_seq(cds, df, as = "AA", start.as.hash = T, stopm1.as.amp = T, 
                          startp1.as.per = T, return.as.list = T)
seqs_codon <- translate_orf_seq(cds, df, as = "codon", start.as.hash = T, stopm1.as.amp = T,
                                startp1.as.per = T, return.as.list = F)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Gene, region, codon and nt counts per library
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
all_counts_cds <- all_counts_leaders <- all_counts_trailers <- study_info <- data.table()
file_names <- c()
cov_all_cds_A_and_P <- data.table()
for (study in study_names) {
  message(study)
  df_real <- read.experiment(study)
  if (study == study_names[1]) {df_merge <- df_real[1,]} else df_merge <- rbind(df_merge, df_real[1,])
  
  file <- filepath(df_real, type = "pshifted")[1]
  study_info <- rbind(study_info, data.table(study, SRR = stringr::str_extract(basename(file), "(SRR|DRR|ERR)[0-9]+")))
  file_names <- c(file_names, file)
  print(basename(file))
  reads <- fimport(file)
  # Region counts
  counts_cds <- countOverlapsW(cds, reads, weight = "score")
  counts_leaders <- countOverlapsW(leaders, reads, weight = "score")
  counts_trailers <- countOverlapsW(trailers, reads, weight = "score")
  all_counts_cds <- cbind(all_counts_cds, counts_cds)
  all_counts_leaders <- cbind(all_counts_leaders, counts_leaders)
  all_counts_trailers <- cbind(all_counts_trailers, counts_trailers)
  # Nucleotide counts
  dt_samp <- coveragePerTiling(cds_with_A_site, reads = reads,
                               is.sorted = TRUE, as.data.table = T)
  cov_all_cds_A_and_P <- cbind(cov_all_cds_A_and_P, dt_samp$count)
  
}
# Save sample info
study_info[, V1 := paste0("R", seq(nrow(study_info)))]
colnames(study_info) <- c("Study ID", "Run ID", "As shown in figures")
fwrite(study_info_path, x = study_info)

all_counts <- all_counts_cds
colnames(cov_all_cds_A_and_P) <- real_lib_names
## Make P and A site tables
cds_lengths <- widthPerGroup(cds, F)
cds_with_A_lengths <- cds_lengths + 3
stopifnot(sum(cds_with_A_lengths) == nrow(cov_all_cds_A_and_P))
gene_split_sites_general <- cumsum(cds_with_A_lengths)
# P sites
gene_split_sites <- c(0, gene_split_sites_general[-length(gene_split_sites_general)])
a_start_codons <- sort(unlist(lapply(c(1,2,3), function(x) x + gene_split_sites)))
dt_all_real <- cov_all_cds_A_and_P[-a_start_codons,]
stopifnot(nrow(dt_all_real) == sum(cds_lengths))
# A sites
gene_split_sites <- gene_split_sites_general
p_stop_codons <- sort(unlist(lapply(c(1,2,3), function(x) gene_split_sites - x)))
dt_all_real_A <- cov_all_cds_A_and_P[-p_stop_codons,]
stopifnot(nrow(dt_all_real_A) == sum(cds_lengths))

# Codon counts
translation_table <- GENETIC_CODE_ORFik(TRUE)
AA_table <- unique(translation_table$AA)
genes_pos_index <- rep.int(seq_along(cds), times = cds_lengths)
AA_real <- seq_usage(dt_all_real, seqs, genes_pos_index, 
                     seqs.order.table = AA_table)
AA_real_Asite <- seq_usage(dt_all_real_A, seqs, genes_pos_index, 
                           seqs.order.table = AA_table)
codon_real <- seq_usage(dt_all_real, seqs_codon, genes_pos_index, 
                        seqs.order.table = translation_table$codon)
codon_real_Asite <- seq_usage(dt_all_real_A, seqs_codon, genes_pos_index, 
                              seqs.order.table = translation_table$codon)
stopifnot(all(AA_real$alpha >= 0 & AA_real_Asite$alpha >= 0))
# Pick seq estimator to use
AA_estimators <- copy(AA_real)
AA_estimators <- AA_estimators[, .(variable, seqs, alpha)]

all_frame_usage <- frame_usage(dt_all_real, relative = T)

all_counts_mrna <- all_counts_cds + all_counts_leaders + all_counts_trailers
colnames(all_counts_cds) <- colnames(all_counts_leaders) <- 
  colnames(all_counts_trailers) <- colnames(all_counts_mrna) <- real_lib_names

# Summarized experiment of cds counts
colData <- DataFrame(libtype = as.factor("RFP"),
                     condition = as.factor("WT"),
                     replicate = as.factor(seq(1, ncol(all_counts))))
colData$SAMPLE <- paste(colData$libtype, colData$condition, colData$replicate, sep = "_")
colnames(all_counts) <- colData$SAMPLE
rownames(colData) <- colnames(all_counts)
res <- SummarizedExperiment(assays=list(counts=as.matrix(all_counts)), rowRanges=cds,
                            colData=colData)
fwrite(all_frame_usage, "/export/valenfs/projects/Hakon/coverageSim/inst/extdata/rnase_bias_human.csv")
fwrite(AA_real[, .(variable, seqs, alpha)], "/export/valenfs/projects/Hakon/coverageSim/inst/extdata/AA_bias_p_site_estimates_human.csv")
fwrite(AA_real_Asite[, .(variable, seqs, alpha)], "/export/valenfs/projects/Hakon/coverageSim/inst/extdata/AA_bias_a_site_estimates_human.csv")
fwrite(codon_real[, .(variable, seqs, alpha)], "/export/valenfs/projects/Hakon/coverageSim/inst/extdata/codon_bias_p_site_estimates_human.csv")
fwrite(codon_real_Asite[, .(variable, seqs, alpha)], "/export/valenfs/projects/Hakon/coverageSim/inst/extdata/codon_bias_a_site_estimates_human.csv")
saveRDS(res, gene_count_top_500_cds_path)

# ECDF of zero inflation for rep1
dt_ecdf_real <- data.table(count = dt_all_real$R1, lib = "Real")
dt_ecdf_real[count > 100, count := 100]
save(list = ls(), file = all_real_result_objects_path)

}
# Merged GRanges object of all real
if (!file.exists(merged_gr_path)) {
  merged_ofst <- ORFik:::ofst_merge(file_names, lib_names = real_lib_names)
  merged_ofst[, c(1,2,3,4,5, 16)] #Columns needed to import
  merged_gr <- GRanges(ORFik:::getGRanges(merged_ofst[, c(1,2,3,4,5, 16)], seqinfo = seqinfo(cds)))
  export.ofst(merged_gr, merged_gr_path)
} else merged_gr <- import.ofst(merged_gr_path, seqinfo = seqinfo(cds))

if (!file.exists(merged_gr_covRle_path)) { 
  merged_gr_covRle <- covRleFromGR(merged_gr)
  ORFik:::export.cov(merged_gr_covRle, merged_gr_covRle_path, seqinfo = seqinfo(merged_gr_covRle))
} else merged_gr_covRle <- fimport(merged_gr_covRle_path)



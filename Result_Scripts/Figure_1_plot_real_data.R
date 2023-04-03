library(ORFik); library(data.table); library(RiboCrypt); library(ggplot2)
source("/export/valenfs/projects/Hakon/coverageSim/Create_simulated_sample_data.R")

# Define output dirs
test_dir <- "/export/valenfs/projects/Hakon/coverageSim/"
data_output_dir <- file.path(test_dir, "data_output")
all_real_result_objects_path <- file.path(data_output_dir, "all_saved_objects.rdata")
all_real_plot_objects_path <- file.path(data_output_dir, "all_saved_objects_plots.rdata")
# Load required data
if (file.exists(all_real_plot_objects_path)) { 
  load(all_real_plot_objects_path)
} else if (file.exists(all_real_result_objects_path)) { 
  load(all_real_result_objects_path)
  source("/export/valenfs/projects/Hakon/coverageSim/Create_simulated_sample_data.R")
  merged_gr <- import.ofst(merged_gr_path, seqinfo = seqinfo(cds))
  merged_gr_covRle <- fimport(merged_gr_covRle_path)
} else stop("You must have run the real data validation and saved in to correct path!")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Figure 1 (Real data features)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# Binned counts over libraries
count_bin <- ggplot(dt_ecdf_real, mapping = aes(x = count)) +
  geom_histogram(bins = 100) + 
  scale_y_continuous(n.breaks = 3) + 
  ylab("Frequency") + xlab("Counts per nt") +
  theme_classic() + coord_cartesian(xlim = c(0, 25)); count_bin 

# Gene region proportions
proportions_cds <- cbind(melt(as.data.table(all_counts_cds / all_counts_mrna)), region = "CDS")
proportions_leader <- cbind(melt(as.data.table(all_counts_leaders /all_counts_mrna)), region = "LEADERS")
proportions_trailer <- cbind(melt(as.data.table(all_counts_trailers / all_counts_mrna)), region = "TRAILERS")
proportions <- rbindlist(list(proportions_cds, proportions_leader, proportions_trailer))
proportions[, region := factor(region, levels = c("LEADERS", "CDS", "TRAILERS"), ordered = T)]
proportions_sum <- proportions[,.(value =  mean(value, na.rm = T)), by = .(variable, region)]
gg_mRNA_regions <- ggplot(proportions, aes(y = round(value*100, 2), x = region)) +
  geom_boxplot(outlier.size = 0.5) +
  ylab("% coverage") +
  xlab("Gene region") +
  scale_y_continuous(breaks=c(25,50,75)) + 
  theme_classic(); gg_mRNA_regions

# Real cor plot
combined_res_subset_real <- coverage_cor(dt_all_real, na.rm = T, rm.selfmatch = T)
combined_res_subset_real[, Comparison := "nt-wise"]

all_counts_proper <- copy(all_counts)
colnames(all_counts_proper) <- paste0("R", seq(1, ncol(all_counts_proper)))
combined_cor_gene_level <- coverage_cor(all_counts_proper, rm.selfmatch = T, na.rm = T)
combined_cor_gene_level[, Comparison := "gene-wise"]
final_cor_real <- rbindlist(list(combined_res_subset_real[, .(Cor, Comparison)], 
                                 combined_cor_gene_level[, .(Cor, Comparison)]))
nt_cor_plot_real <- ggplot(final_cor_real, aes(x = Cor, fill = Comparison)) +
  geom_density(alpha = 0.75) +
  ylab("Density") + 
  xlab("Correlation") + 
  theme_classic()  +
  scale_x_continuous(breaks = c(0.00, 0.50, 1.0), limits = c(-0.1, 1)) +
  scale_fill_manual(values = c("#F7E8D4", "#E8E6E6")) + 
  theme(legend.position="top",
        legend.title=element_text(size=8, margin = margin(r=2,l=0,t=0,b=0)),
        legend.text=element_text(size=7), 
        legend.key.size = unit(0.2, "cm")); nt_cor_plot_real
# Peak analysis
cds_seqs <- translate(txSeqsFromFa(cds, df, TRUE))
subseq(cds_seqs, 1, 1) <- "#"
lt2 <- width(cds_seqs) > 2
subseq(cds_seqs[lt2], width(cds_seqs[lt2]) - 1, width(cds_seqs[lt2]) - 1) <- "&"
cds_codon_usage <- table(strsplit(as.character(unlist(cds_seqs, use.names = FALSE)), split = ""))

peaks_max <- ORFik:::findPeaksPerGene(cds, merged_gr_covRle, top_tx = 0, type = "max")
peaks_max[, gene_length := widthPerGroup(cds)[gene_id]]
peaks_max[, position_rel := round((position / gene_length)*100, 1)]
sum(peaks_max$position == 1); sum(peaks_max$position == peaks_max$gene_length - 5); sum(peaks_max$position %% 3 == 1) / length(cds)
location_of_maxpeaks <- ggplot(peaks_max, aes(x = position_rel)) + 
  geom_histogram(bins = 100) + theme_classic() + xlab("Position (100 bins)") +
  ylab("Density\n(max peaks)"); location_of_maxpeaks
codon_max <- windowPerGroup(startSites(cds[peaks_max$gene_id], T, T, T), mrna[peaks_max$gene_id],
                            upstream = -(peaks_max$position -1),
                            downstream = peaks_max$position + 1); table(widthPerGroup(codon_max))
irl_max <- split(IRanges(start = floor(peaks_max$position/3) + 1, width = 1,
                         names = peaks_max$gene_id), peaks_max$gene_id)
irl_max_A <- split(IRanges(start = floor(peaks_max$position/3) + 2, width = 1,
                           names = peaks_max$gene_id), peaks_max$gene_id)
max_codon_usage <- table(strsplit(as.character(unlist(cds_seqs[irl_max], use.names = FALSE)), split = ""))
max_codon_usage_rel <- (max_codon_usage / cds_codon_usage[names(max_codon_usage)])
max_codon_usage_rel <- round(sort(max_codon_usage_rel / sum(max_codon_usage_rel) * 100))

peaks <- ORFik:::findPeaksPerGene(cds, merged_gr_covRle, top_tx = 0, type = "median")
peaks[, gene_length := widthPerGroup(cds)[gene_id]]
peaks[, position_rel := round((position / gene_length)*100, 1)]
location_of_medianpeaks <- ggplot(peaks, aes(x = position_rel)) + 
  geom_histogram(bins = 100) + theme_classic() + xlab("Position (binned to 100)") + ylab("Max peak coverages"); location_of_medianpeaks
codon_peak <- windowPerGroup(startSites(cds[peaks$gene_id], T, T, T), mrna[peaks$gene_id],
                             upstream = -(peaks$position - 1),
                             downstream = peaks$position + 1); table(widthPerGroup(codon_peak))
seqs_peak <- txSeqsFromFa(codon_peak, df, is.sorted = T)
sort(table(seqs_peak), decreasing = T)
sort(table(translate(seqs_peak)), decreasing = T)

window_size_peak <- 100
cds_ir <- windowPerGroup(startSites(cds[peaks$gene_id], T, T, T), mrna[peaks$gene_id],
                         upstream = -(peaks$position - window_size_peak),
                         downstream = peaks$position + window_size_peak)
summary(widthPerGroup(cds_ir))
cds_ir <- cds_ir[widthPerGroup(cds_ir, FALSE) == (window_size_peak*2 + 1)]; length(cds_ir)
cov_real_ir <- coveragePerTiling(cds_ir, merged_gr_covRle, is.sorted = T, as.data.table = T, withFrames = T)

cov_real_ir_psites <- cov_real_ir[frame == 0,]
cov_real_ir_psites_zscore <- coverageScorings(cov_real_ir_psites, scoring = "zscore")
cov_real_ir_psites_zscore[, rel_pos := position - window_size_peak]
peak_plot <- ggplot(cov_real_ir_psites_zscore, mapping = aes(x = rel_pos, y = score)) + 
  geom_line(color = "#DDA28F", size = 1) + 
  geom_vline(xintercept = c(-27, 27), color = "gray", size = 0.5, alpha = 0.5, linetype = "dotted") + 
  theme_classic() + 
  ylab("Coverage\n(Zscore)") + 
  xlab("Position: peak=0"); peak_plot

# Real frame usage
all_frame_usage[, Data := substring(variable, 1, 1)]
all_frame_usage[, percentage_usage := round((sum/sum(sum))*100, 1), by = variable]
frame_plot_real <- ggplot(all_frame_usage, aes(x = frame, y = percentage_usage, fill = frame)) +
  geom_boxplot(alpha = 0.95) + 
  ylab("% coverage") + 
  xlab("Frame") +
  scale_y_continuous(breaks = c(25,50,75)) +
  theme_classic() +  theme(legend.position="none"); frame_plot_real
# ggheatmap_frame <- ggplot(all_frame_usage_m[Data == "R",],
#                           aes(frame, variable, fill = percentage_usage)) +
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(0,5), space = "Lab", 
#                        name="Rel. Cov.") +
#   theme_classic()+ # minimal theme
#   theme(legend.position="none") +
#   ylab("Ribo-seq library") +
#   xlab("Frame"); ggheatmap_frame
all_AA_to_plot <- copy(AA_real)
all_AA_to_plot[, relative_to_max_score := (mean_txNorm_percentage /
                                             max(mean_txNorm_percentage)*100), by = variable]
# IDEA: clustering
ggheatmap_aa <- ggplot(all_AA_to_plot, aes(seqs, variable, fill = relative_to_max_score)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(0,100), space = "Lab", 
                       name="Rel. Cov.") +
  theme_classic()+ # minimal theme
  theme(legend.position="none") +
  ylab("Ribo-seq library") + xlab("Sequence motif"); ggheatmap_aa

all_codon_real_Asite <- copy(codon_real_Asite)
all_codon_real_Asite[, relative_to_max_score := (mean_txNorm_percentage /
                                                   max(mean_txNorm_percentage)*100), by = variable]
all_codon_real_Asite[, type := "A"]
ggheatmap_codon_A <- ggplot(all_codon_real_Asite, aes(seqs, variable, fill = relative_to_max_score)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(0,100), space = "Lab", 
                       name="Rel. Cov.") +
  theme_classic()+ # minimal theme
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Ribo-seq library"); ggheatmap_codon_A

# real_data_plot <- gridExtra::arrangeGrob(nt_cor_plot_real, gg_mRNA_regions, frame_plot_real, 
#                                          location_of_maxpeaks, count_bin, peak_plot,
#                                          ggheatmap_aa,
#                                          layout_matrix = t(matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3)))
# plot(real_data_plot)
real_data_plot <- cowplot::plot_grid(cowplot::plot_grid(nt_cor_plot_real, gg_mRNA_regions, frame_plot_real, labels = "AUTO", label_size = 12, nrow = 1),
                           cowplot::plot_grid(location_of_maxpeaks, count_bin, peak_plot, labels = c("D", "E", "F"), label_size = 12, nrow = 1),
                           cowplot::plot_grid(ggheatmap_aa, labels = c("G"), label_size = 12),
                           ncol = 1, rel_heights = c(1,1,1.2))
plot(real_data_plot)
ggsave(file.path(test_dir, "figure2.png"), real_data_plot, width = 7.5, heigh = 5.7, dpi = 600)
ggsave(file.path(test_dir, "max_peak_locations.png"), location_of_maxpeaks, width = 6, heigh = 3.5, dpi = 600)
ggsave(file.path(test_dir, "codon_A_site_real_libs.png"), ggheatmap_codon_A, width = 7, heigh = 3.5, dpi = 600)
save(list = ls(), file = all_real_plot_objects_path)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Figure 3 (Simulated data comparison to real)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Now run simulation code
# Then plot ->
ecdf_probs <- seq(0, 1, 0.05)
dt_ecdf_sim <- data.table(count = dt_all_sim$S1, lib = "Sim")
dt_ecdf <- rbindlist(list(dt_ecdf_real, dt_ecdf_sim))
dt_ecdf[count > 100, count := 100]
dt_ecdf <- dt_ecdf[, quantile(count, probs = ecdf_probs), by = lib]
colnames(dt_ecdf) <- c("lib", "count")
dt_ecdf[, quantile := rep(ecdf_probs*100, 2)]

ecdf_plot <- ggplot(dt_ecdf, aes(x = count, y = quantile, color = lib)) +
  geom_line(alpha = 0.7, size = 0.3) + 
  scale_y_continuous(breaks=c(50, 100)) + 
  scale_x_continuous(breaks=c(50, 100)) + 
  scale_color_manual(values = c("gray", "black")) +
  ylab("Quantile") + xlab("Count") +
  theme_classic() +  theme(legend.position="none"); ecdf_plot 
qq_plot <- ggplot(mapping = aes(x = dt_ecdf[lib == "Real",]$count, y = dt_ecdf[lib == "Sim",]$count)) +
  geom_line(alpha = 0.7, size = 1) + 
  scale_y_continuous(breaks=c(50, 100)) + 
  scale_x_continuous(breaks=c(50, 100)) + 
  ylab("Q. Count (S)") + xlab("Q. Count (R)") +
  theme_classic() +  theme(legend.position="none"); qq_plot 

# Example coverage from selected gene

name <- "ENST00000380191"
cov_real <- coveragePerTiling(cds[name], fimport(filepath(read.experiment(study_names[1])[1,], "pshifted")), is.sorted = T, as.data.table = T, withFrames = T, fraction = "Real")
cov_sim <- coveragePerTiling(cds[name], fimport(filepath(df_sim3[1,], "pshifted")), is.sorted = T, as.data.table = T, withFrames = T, fraction = "Sim")
cov_combined <- rbindlist(list(cov_real, cov_sim))
cov_combined[, frame := factor(frame)]
cov_plot <- ggplot(cov_combined, aes(y = count, 
                                     x = position, fill = frame)) + geom_col() + 
  scale_y_continuous(breaks=c(10)) + 
  facet_wrap(~ fraction, ncol = 1, strip.position = "right") + 
  theme_classic() + theme(legend.position = "none"); cov_plot
cov_plot_no_colours <- ggplot(cov_combined, aes(y = count, 
                                                x = position)) + geom_col() + 
  scale_y_continuous(breaks=c(10)) + 
  facet_wrap(~ fraction, ncol = 1, strip.position = "right") + 
  theme_classic() + theme(legend.position = "none"); cov_plot_no_colours

# Frame plot merged

View(data.table(all_frame_usage, O = all_frame_usage_f$frame_usage, S = all_frame_usage_s$frame_usage))
all_frame_usage_f[, Data := substring(variable, 1, 1)]
all_frame_usage_f[, percentage_usage := round((V1/sum(V1))*100, 1), by = variable]
all_frame_usage_m <- rbindlist(list(all_frame_usage, all_frame_usage_f))
colors_frame_m <- factor(c("#F8766D", "#00BA38", "#619CFF"), levels = c("#F8766D", "#00BA38", "#619CFF"),  ordered = T)[all_frame_usage_m$frame]

frame_plot_m <- ggplot(all_frame_usage_m, aes(x = frame, y = percentage_usage, fill = interaction(frame, Data))) +
  scale_fill_manual(values = c("#F8666D", "#00BA38", "#619CFF", "#C8363D", "#007A06", "#517CCF")) + 
  geom_boxplot(alpha = 0.95) + 
  ylab("% coverage") + xlab("Frame") +
  theme_classic() +  theme(legend.position="none"); frame_plot_m


final_plot <- gridExtra::arrangeGrob(nt_cor_plot, qq_plot, frame_plot_m, cov_plot_no_colours, ecdf_plot,
                                     layout_matrix = t(matrix(c(1,2,3,1,5,5,4,4,4), nrow = 3)))
plot(final_plot)
ggsave(file.path(test_dir, "figure3.png"), final_plot, width = 6, height = 5, dpi = 600)



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Supplementary figures
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Test global shapes
cds_150 <- windowPerGroup(startSites(cds, T, T, T), mrna, 0, 151)
cds_150_stop <- windowPerGroup(stopSites(cds, T, T, T), mrna, 150, 9)
table(widthPerGroup(cds_150, F))
table(widthPerGroup(cds_150_stop, F))
dt_150_start <- metaWindow(x = merged_gr_covRle, windows = cds_150, 
                           scoring = "transcriptNormalized", fraction = "start", withFrames = T,
                           zeroPosition = 0)
dt_150_start[, frame := as.factor(frame)]
gg_150_start <- ggplot(dt_150_start, aes(x = position, y = score, fill = frame)) + 
  geom_col() + theme_classic() + ylab("Coverage") + xlab("Position (nt)") + theme(legend.position = "top"); gg_150_start

dt_150_stop <- metaWindow(x = merged_gr_covRle, windows = cds_150_stop,
                          scoring = "transcriptNormalized", fraction = "start",
                          withFrames = T, zeroPosition = 151)
dt_150_stop[, frame := as.factor(frame)]
gg_150_stop <- ggplot(dt_150_stop, aes(x = position, y = score, fill = frame)) + 
  geom_col() + theme_classic() + ylab("Coverage") + xlab("Position (nt)") + theme(legend.position = "top")
plot_150 <- gridExtra::arrangeGrob(gg_150_start, gg_150_stop, nrow = 1); plot(plot_150)
ggsave(file.path(test_dir, "cds_150_each_end.png"), plot_150, width = 6, heigh = 3)

fpkms <- ORFik:::fpkm_calc(counts_cds, widthPerGroup(cds, F), sum(counts_cds))
names(fpkms) <- names(cds)
fpkms <- sort(fpkms, decreasing = T)
i <- 4
#RiboCrypt::multiOmicsPlot_ORFikExp(cds[names(fpkms[i])], df_sim_O2[6:7,], frames_type = "columns", annotation = cds[names(fpkms[i])])
RiboCrypt::multiOmicsPlot_ORFikExp(cds[names(fpkms[i])], df_sim_O2, frames_type = "columns", annotation = cds[names(fpkms[i])])
RiboCrypt::multiOmicsPlot_ORFikExp(cds[names(fpkms[i])], df_merge, frames_type = "columns", annotation = cds[names(fpkms[i])])
A <- (c(5,2,1)*1000)
sd(A)/mean(A)
i <- 4
aab <- coveragePerTiling(cds[names(fpkms[i])], merged_gr_covRle, as.data.table = T)
plot(decompose(ts(aab$count, frequency = 3), type = "multiplicative"))
aab <- cbind(fft_dt, dt_all_real)
aab <- aab[genes == 2 & position < 150,7:9]

aab <- coveragePerTiling(cds, merged_gr_covRle, as.data.table = T)
aab[,position_rev := rev(position), by = genes]
aab_start150NT <- coverageScorings(aab[position < 150,], "transcriptNormalized")
aab_stop150NT <- coverageScorings(aab[position_rev < 150,][, position := position_rev], "transcriptNormalized")
plot(decompose(ts(aab_stop150NT$score, frequency = 3), type = "multiplicative"))
png(file.path(test_dir, "decomposed_signal_150NT_start.png"), width = 6, height = 5, units = "in", res = 300)
plot(decompose(ts(aab_start150NT$score, frequency = 3), type = "multiplicative"))
dev.off()

plot(decompose(ts(aab, frequency = 9), type = "multiplicative"))
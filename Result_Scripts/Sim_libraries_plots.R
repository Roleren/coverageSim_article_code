library(coverageSim);library(ORFik); library(data.table); library(RiboCrypt)
library(ggplot2); library(gridExtra); library(cowplot)

# Define output dirs
test_dir <- "/export/valenfs/projects/Hakon/coverageSim/"
data_output_dir <- file.path(test_dir, "data_output")
dt_all_list_named <- readRDS(file.path(data_output_dir, "all_coverages.rds"))
dt_all_real <- dt_all_list_named[["real"]]
cds_top500 <- rowRanges(readRDS(file.path(data_output_dir, "genecount_top_500_cds.rds")))
cds <- cds_top500
cds_lengths <- widthPerGroup(cds_top500, FALSE)
genes <- genes_pos_index <- rep.int(seq_along(cds_lengths), times = cds_lengths)

theme_plot <- theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(),
                    legend.position = "bottom", strip.background = element_blank(),
                    legend.title=element_text(size=9, margin = margin(r=0,l=0,t=0,b=0)),
                    legend.text=element_text(size=8),
                    legend.key.size = unit(0.3, "cm"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Variance analysis
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
var_dt_all <- data.table(lib_size = colSums(dt_all_real))
var_list <- dt_all_list_named
for (i in seq_along(var_list)) var_dt_all <- cbind(var_dt_all, data.table(unlist(lapply(var_list[i][[1]], function(x) sd(x)/mean(x)))))
colnames(var_dt_all)[-1] <- names(var_list)

var_dt_all_melt <- melt(var_dt_all[,-1])
levels <- var_dt_all_melt[, .(median = median(value)), by = variable][order(median, decreasing = T),]
var_dt_all_melt[, variable := factor(variable, levels = levels$variable, ordered = T)]
var_dt_main <- var_dt_all_melt[variable %in% c("real", "no_features", "optimized"),]
var_plot <- ggplot(var_dt_all_melt, aes(y = log2(value), fill = variable)) + 
  geom_hline(yintercept = log2(levels[variable == "real",]$median)) + 
  geom_boxplot() + theme_classic() + ylab("log2 Relative SD"); var_plot

var_dt_main[variable == "real",]$variable <- "Real"
var_dt_main[variable == "optimized",]$variable <- "Simulated"
var_dt_main[variable == "no_features",]$variable <- "Null model"
var_density_plot <- ggplot(var_dt_main, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.75) +
  ylab("Density") + 
  xlab("Coefficient of variation") + 
  theme_classic()  +
  scale_fill_brewer(palette = 16) + 
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.3, "cm")); var_density_plot



ggsave(file.path(test_dir, "variance_analysis.png"), var_plot, width = 6, heigh = 4, dpi = 600)
cor(var_dt_all[,-1])

# dispersion <- data.table(lib_size = colSums(dt_all_real))
# for (i in seq_along(var_list)) dispersion <- cbind(dispersion, data.table(unlist(lapply(var_list[i][[1]], function(x) (mean(x)^2) / (var(x)-mean(x)) + 0.01))))
# colnames(dispersion)[-1] <- names(var_list)
# dispersion_melt <- melt(dispersion[,-1])
# levels <- dispersion_melt[, .(median = median(value)), by = variable][order(median, decreasing = F),]
# dispersion_melt[, variable := factor(variable, levels = levels$variable, ordered = T)]
# dispersion_plot <- ggplot(dispersion_melt, aes(y = log2(value), fill = variable)) + 
#   geom_hline(yintercept = log2(levels[variable == "real",]$median)) + 
#   geom_boxplot() + theme_classic() + ylab("Dispersion"); dispersion_plot

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Skewness
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
skew_dt <- lapply(dt_all_list_named, function(x) coverageSim:::skewness(x))
skew_dt_final <- data.table(skew = unlist(skew_dt, F), 
                            sample = rep(names(skew_dt), each = unique(lengths(skew_dt))))
levels_skew <- skew_dt_final[, .(median = median(skew)), by = sample][order(median, decreasing = F),]
skew_dt_final[, sample := factor(sample, levels = levels_skew$sample, ordered = T)]
skew_plot <- ggplot(skew_dt_final, aes(y = skew, fill = sample)) + 
  geom_hline(yintercept = log2(levels_skew[sample == "real",]$median)) + 
  geom_boxplot() + theme_classic() + ylab("Skew"); skew_plot
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Frame usage
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
frame_usage_all_sim <- lapply(dt_all_list_named[c("optimized", "real")], function(x) frame_usage(x, normalize_to = 100))
frame_usage_all_sim_dt <- rbindlist(frame_usage_all_sim)
frame_usage_all_sim_dt[, data := substr(variable, 1, 1)]
frame_usage_all_sim_dt_sub <- frame_usage_all_sim_dt[data %in% c("R", "O"),]
frame_plot_all <- ggplot(frame_usage_all_sim_dt_sub, 
                         aes(x = interaction(frame, data), y = frame_usage, fill = frame)) +
  geom_boxplot(alpha = 0.95) + 
  ylab("% coverage") + 
  xlab("Frame") +
  scale_y_continuous(breaks = c(25,50,75)) +
  theme_classic() +  theme(legend.position="none"); frame_plot_all

# all_frame_usage_s <- frame_usage_all_sim[["default"]]
# all_frame_usage_f <- frame_usage_all_sim[["periodicity"]]
# all_frame_usage_o <- frame_usage_all_sim[["optimized"]]

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# auto correlation (codon and nt) (codon to next codon correlation)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#BiocParallel::bp
codon_nt_acf <- coverageSim:::auto_correlation_genes_all(dt_all_list_named[c("real")],
                                           100, by.codon = T, genes = genes,
                                           codon.vs.nt = T, fun = acf)
codon_nt_acf_opt <- coverageSim:::auto_correlation_genes_all(dt_all_list_named[c("optimized")],
                                               100, by.codon = T, genes = genes,
                                               codon.vs.nt = T, fun = acf)
# nt_acf <- auto_correlation_genes_all(dt_all_list_named[ac_names],
#                                      100, by.codon = F, genes = genes,
#                                      codon.vs.nt = T, fun = acf)
# nt_pacf <- auto_correlation_genes_all(dt_all_list_named[c("real")],
#                                       100, by.codon = F, genes = genes,
#                                       codon.vs.nt = T, fun = pacf)

# codon_acf <- auto_correlation_genes_all(dt_all_list_named[ac_names],
#                                         33, by.codon = T, genes = genes,
#                                         codon.vs.nt = F, fun = acf)
# codon_pacf <- auto_correlation_genes_all(dt_all_list_named[c("real")],
#                                          33, by.codon = T, genes = genes,
#                                          codon.vs.nt = F, fun = pacf)
# period_detector(codon_acf[id == "real",], quantile_min = 0.8, line_center = 9)
# period_detector(codon_acf[id == "optimized",], quantile_min = 0.8, line_center = 9)
# period_detector(codon_pacf, quantile_min = 0.9, line_center = 9, xlims = c(3, 30))
# nt_acf_plot  <- ac_boxplot(nt_acf, autocor_name = "nt")
# nt_acf_plot  <- ac_boxplot(nt_acf[id == "real",], autocor_name = "nt")
# nt_acf_plot  <- ac_boxplot(nt_acf[id == "optimized",], autocor_name = "nt")
# nt_pacf_plot <- ac_boxplot(nt_pacf, autocor_name = "nt")
# 
# codon_acf_plot  <- ac_boxplot(codon_acf[sample == 6,])
# codon_pacf_plot <- ac_boxplot(codon_pacf[sample == 1 & mean_Cor > 0.001,])

codon_nt_acf_plot  <- ac_boxplot(codon_nt_acf[mean_Cor > quantile(Cor, 0.92),])
codon_nt_acf_opt_plot  <- ac_boxplot(codon_nt_acf_opt[mean_Cor > quantile(Cor, 0.8),])
codon_to_codon_cor <- codon_nt_acf_plot
ggsave(file.path(test_dir, "codon_to_codon_correlation.png"), codon_to_codon_cor, width = 6, heigh = 3.5, dpi = 600)

# plot(c_c_comp[id == "real", mean(Cor), by = distance]$V1)
# ggplot(c_c_comp[id == "real",], aes(y = Cor, x = distance)) + geom_boxplot() + 
#   xlab(paste("Upstream", autocor_name)) + theme_classic() + facet_wrap(~ Var1_index)
# lapply(dt_all_list_named[["no_features"]], function(x) Box.test(x, lag = 1)$p.value)
# 
# value <- autocor_window(as.vector(
#   extraDistr::rdirmnom(1, 1000, 
#                        sample(c(1, 0.5), 300, replace = T, prob = c(0.1, 0.9))))
#   , 30, padding.rm = T); barplot(value); acf(as.vector(value), lag.max = 30)
# sample <- sample(c(1, 0.5), 300, replace = T, prob = c(0.1, 0.9)); barplot(sample)
# sample_mean <- mean(sample)
# sample <- autocor_window(sample, 30, padding.rm = T); sample <- sample * (sample_mean/ mean(sample)); barplot(sample)
# sample <- as.vector(extraDistr::rdirmnom(1, 1000, sample)); barplot(sample)
# acf(sample)
# 
# sample <- c(5, sample(c(1, 0.5), 300, replace = T, prob = c(0.1, 0.9)))*3; sample_mean <- mean(sample); barplot(sample)
# sample <- as.vector(extraDistr::rdirmnom(1, 1000, sample)); barplot(sample)
# sample <- autocor_window(sample, 30, padding.rm = T); barplot(sample)
# sample <- sample * (sample_mean/ mean(sample));
# sample <- as.vector(extraDistr::rdirmnom(1, 1000, sample)); barplot(sample)
# barplot(acf(sample, lag.max = 40)$acf[-1])
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Spectral analysis (Periodicity)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
fft_dt <- data.table(genes = rep.int(seq(length(cds_top500)), times = widthPerGroup(cds_top500, F)))
fft_dt[, position := seq.int(.N), by = genes]

spectral_analysis <- lapply(dt_all_list_named, function(x) {
  fft <- copy(fft_dt)
  fft[, count := x[, 2]]
  fft <- fft[position < 150,]
  fft <- coverageScorings(fft, "transcriptNormalized")
  fft[, spec.pgram(score, plot = F)[c("spec", "freq")]]
})
fft_dt_all <- rbindlist(spectral_analysis)
fft_dt_all[, id := rep(names(spectral_analysis), each = nrow(spectral_analysis[[1]]))]
fft_dt_3 <- fft_dt_all[id %in% c("real", "optimized", "no_features")]
fft_dt_3[id == "real",]$id <- "Real"
fft_dt_3[id == "optimized",]$id <- "Simulated"
fft_dt_3[id == "no_features",]$id <- "Null model"
fft_dt_3[, id := factor(id, unique(fft_dt_3$id)[c(1,3,2)], ordered = TRUE)]
fft_dt_3[, Period := round(1/freq, 0)]
gg_periods <- ggplot(fft_dt_3, aes(x = 1/freq, y = spec)) +
  geom_line() + xlab("Period") + ylab("Amplitude") + 
  scale_x_continuous(limits = c(0,50), breaks = c(3, 20, 40)) +
  theme_classic() + facet_wrap(~ id, scales = "free_y", ncol = 1) + 
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), 
        legend.position = "bottom", strip.background = element_blank(),
        legend.title=element_text(size=9, margin = margin(r=0,l=0,t=0,b=0)),
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.3, "cm"));  gg_periods
# gg_periods <- ggplot(fft_dt_3, aes(x = as.factor(round(1/freq, 1)), y = id, fill = spec)) +
#   geom_tile(color = "gray") + xlab("Period") + ylab("Data") + 
#   theme_classic() + 
#   scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
#                        midpoint = 3, limit = c(0, 400), space = "Lab", 
#                        name="Rel. Cov.") + 
#   xlim("2","3","4", "6", "12","150") + 
#   theme(legend.position = "none");  gg_periods

ggsave(file.path(test_dir, "FFT_all_analysis.png"), gg_periods, width = 6, heigh = 4, dpi = 600)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TIS / TES
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
binned_cov <- function(count, scaleTo = 100, scoring = "meanPos") {
  count[, max_pos := max(position), by = genes]
  count[, `:=`(scalingFactor, (scaleTo/max_pos))]
  count[, `:=`(position, ceiling(scalingFactor * position))]
  if (length(scaleTo) == 1) {
    count[position > scaleTo, `:=`(position, scaleTo)]
  } else {
    if (any(count[, .(max = max(position)), by = genes]$max > 
            scaleTo)) {
      index <- count[, .(ind = which.max(position)), by = genes]$ind
      count[index, `:=`(position, scaleTo)]
    }
  }
  return(coverageScorings(count, scoring, copy.dt = FALSE))
}

sim_ORFs_full_length <- data.table(count = rowMeans(dt_all_list_named[["optimized"]]), genes)
sim_ORFs_full_length[,position := seq.int(.N), by = genes]
TIS_O <- sim_ORFs_full_length[position < 100,]
TIS_O <- TIS_O[, .(score = sum(count)), by = position]
barplot(TIS_O$score, col = c("red", "green", "blue"))
sim_ORFs_binned <- copy(sim_ORFs_full_length)
sim_ORFs_binned <- binned_cov(sim_ORFs_binned)
sim_ORFs_binned <- coverageScorings(sim_ORFs_binned, "transcriptNormalized")
barplot(sim_ORFs_binned$score)
sim_ORFs_full_length[, max_pos := max(position), by = genes]
sim_ORFs_full_length <- sim_ORFs_full_length[max_pos <= 1000,]
sim_ORFs_full_length <- coverageScorings(sim_ORFs_full_length, "transcriptNormalized")
barplot(sim_ORFs_full_length$score)

real_ORFs_full_length <- data.table(count = rowMeans(dt_all_list_named[["real"]]), genes)
real_ORFs_full_length[,position := seq.int(.N), by = genes]
TIS_R <- real_ORFs_full_length[position < 100,]
TIS_R <- TIS_R[, .(score = sum(count)), by = position]
real_ORFs_binned <- copy(real_ORFs_full_length)
real_ORFs_binned <- binned_cov(real_ORFs_binned)
real_ORFs_binned <- coverageScorings(real_ORFs_binned, "transcriptNormalized")
barplot(real_ORFs_binned$score)
barplot(TIS_R$score, col = c("red", "green", "blue"))
real_ORFs_full_length[, max_pos := max(position), by = genes]
real_ORFs_full_length <- real_ORFs_full_length[max_pos <= 1000,]
real_ORFs_full_length <- coverageScorings(real_ORFs_full_length, "transcriptNormalized")
barplot(real_ORFs_full_length$score)

TIS <- data.table(rowMeans(dt_all_list_named[["real"]]), genes)
TIS[,position := seq.int(.N), by = genes]
TIS <- TIS[position < 100,]
TIS_R <- TIS[, .(score = sum(V1)), by = position]
barplot(TIS_R$score, col = c("red", "green", "blue"))

TIS <- data.table(rowMeans(dt_all_list_named[["no_features"]]), genes)
TIS[,position := seq.int(.N), by = genes]
TIS <- TIS[position < 100,]
TIS <- TIS[, .(score = sum(V1)), by = position]
TIS_N <- copy(TIS)
barplot(TIS_O$score, col = c("red", "green", "blue"))

TIS_merged <- rbindlist(list(data.table(TIS_O, id = "Simulated", frame = rep(seq(3), length.out = nrow(TIS_O))),
                             data.table(TIS_R, id = "Real", frame = rep(seq(3), length.out = nrow(TIS))),
                             data.table(TIS_N, id = "Null model", frame = rep(seq(3), length.out = nrow(TIS)))))
TIS_merged[, frame := as.factor(frame)]
TIS_merged[, id := factor(id, levels = c("Real", "Simulated", "Null model"), ordered = TRUE)]
gg_TIS <- ggplot(TIS_merged, aes(x = position, y = score, fill = frame)) + 
  geom_col() + theme_classic() + facet_wrap(~ id, ncol = 1, scales = "free") +
  ylab("Meta coverage") + xlab("Position") + 
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), 
        legend.position = "right", strip.background = element_blank(),
        legend.title=element_text(size=9, margin = margin(r=0,l=0,t=0,b=0)),
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.3, "cm")) +
  guides(fill=guide_legend(title="Frame")); gg_TIS


nrows <- nrow(sim_ORFs_full_length)
ORFs_full_length <- rbindlist(list(
  data.table(sim_ORFs_full_length, id = "Simulated", length = "Full length", frame = rep(seq(3), length.out = nrows)),
  data.table(real_ORFs_full_length, id = "Real", length = "Full length", frame = rep(seq(3), length.out = nrows)),
  data.table(sim_ORFs_binned, id = "Simulated", length = "Binned length", frame = 4),
  data.table(real_ORFs_binned, id = "Real", length = "Binned length", frame = 4)))
ORFs_full_length[, id := factor(id, levels = c("Real", "Simulated", "Null model"), ordered = TRUE)]
ORFs_full_length[, frame := factor(frame, levels = c("1", "2", "3", ""))]
gg_full_ORF <- ggplot(ORFs_full_length, aes(x = position, y = score, fill = frame)) + 
  geom_col() + theme_classic() + facet_wrap(id ~ length, ncol = 2, scales = "free") +
  ylab("Meta coverage") + xlab("Position") + 
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), 
        legend.position = "right", strip.background = element_blank(),
        legend.title=element_text(size=9, margin = margin(r=0,l=0,t=0,b=0)),
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.3, "cm")) +
  guides(fill=guide_legend(title="Frame")); gg_full_ORF
ggsave(file.path(test_dir, "supplementary_metacoverage_cds.png"), gg_full_ORF, width = 6, heigh = 4, dpi = 600)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Max peak
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
peaks_max_sim <- lapply(split(dt_all_list_named[c("real","optimized", "no_features")], c(1,2,3)), function(x) data.table(rowMeans(x[[1]]), genes, id = names(x[1])))
peaks_max_sim <- rbindlist(peaks_max_sim)
peaks_max_sim[,position := seq.int(.N), by = .(id, genes)]
peaks_max_sim[, gene_length := widthPerGroup(cds)[genes]]
peaks_max_sim[, position_rel := round((position / gene_length)*100, 1)]
peaks_max_sim[, max_peak := max(V1), by = .(id, genes)]
peaks_max_sim <- peaks_max_sim[V1 == max_peak,]
peaks_max_sim[, number_of_hits := sum(.N), by = .(id, genes)]
peaks_max_sim[, id := factor(chmatch(id, unique(id)), labels = c("Real", "Simulated", "Null model"), levels = seq_along(unique(id)), ordered = TRUE)]
sum(peaks_max_sim$position == 1); sum(peaks_max_sim$position == peaks_max_sim$gene_length - 5); sum(peaks_max_sim$position %% 3 == 1) / length(cds)
location_of_maxpeaks_sim <- ggplot(peaks_max_sim, aes(x = position_rel)) + 
  geom_histogram(bins = 100) + theme_classic() + facet_wrap(~ id, scales = "free_y", ncol = 1) + 
  xlab("Position (100 bins)") + ylab("Density\n(max peaks)") +
  theme_plot; location_of_maxpeaks_sim
  

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# AA analysis (Skipped for now)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
df <- read.experiment("PRJEB26593")
AA_table <- unique(coverageSim:::GENETIC_CODE_ORFik(TRUE)$AA)

seqs <- coverageSim:::translate_orf_seq(cds, df, as = "AA", start.as.hash = T, stopm1.as.amp = T, 
                          startp1.as.per = T, return.as.list = T)
AA_opt <- coverageSim:::seq_usage(dt_all_list_named[["optimized"]], seqs, genes_pos_index, 
                    seqs.order.table = AA_table)
AA_opt[, relative_to_max_score := (mean_txNorm_percentage /
                                     max(mean_txNorm_percentage)*100), by = variable]
AA_opt[, variable := factor(variable, labels = unique(gsub("O", "S", variable)))]
ggheatmap_aa_opt <- ggplot(AA_opt, aes(seqs, variable, fill = relative_to_max_score)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(0,100), space = "Lab", 
                       name="Rel. Cov.") +
  theme_classic()+ # minimal theme
  theme(legend.position="none") +
  ylab("Ribo-seq library") + xlab("Codon"); ggheatmap_aa_opt


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# All sim analysis plots
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
final_sim_plot <- plot_grid(plot_grid(gg_TIS, gg_periods, labels = "AUTO", label_size = 11),
                                     plot_grid(location_of_maxpeaks_sim, var_density_plot, labels = c("C", "D"), label_size = 11),
                                     plot_grid(ggheatmap_aa_opt, labels = "E", label_size = 11),
                                     rel_heights = c(3,2,2), ncol = 1)
#layout_matrix = t(matrix(c(1,2,3,1,5,5,4,4,4), nrow = 3))
plot(final_sim_plot)
ggsave(file.path(test_dir, "all_analysis_sim.png"), final_sim_plot, width = 6, heigh = 7, dpi = 600)



# df_test <- df_merge
# df_test <- df_codonAndperiodicity
# envExp(df_test) <- new.env()
# i <- 4
# RiboCrypt::multiOmicsPlot_ORFikExp(cds[i], df = df_test[c(2, 7),],
#                                    frames_type = "columns",annotation = cds[i],
#                                    BPPARAM = BiocParallel::SerialParam())

#tai_table_start <- file.path(test_dir, "tAI_table/tAI_human_with_start.csv")
# df <- read.experiment("PRJEB26593")
# reads <- fimport(filepath(df[1,], "pshifted"))
# cds <- loadRegion(df, "cds", names.keep = filterTranscripts(df))
# cds_hits <- countOverlapsW(cds, reads, "score")
# 
# region_ranges <- cds[which.max(cds_hits)]
# region_ranges_noname <- region_ranges
# names(region_ranges_noname) <- NULL
# tile <- tile1(region_ranges_noname, sort.on.return = TRUE)
# dt_range <- data.table::data.table(seqnames = as.character(unlist(seqnames(tile), use.names = FALSE)),
#                                    start = unlist(start(tile), use.names = FALSE),
#                                    end = unlist(start(tile), use.names = FALSE),
#                                    strand = as.character(unlist(strand(tile), use.names = FALSE)))
# seqs <- translate(txSeqsFromFa(region_ranges, df, is.sorted = TRUE))

# tAI <- fread(tai_table_start)
# start_codon_is_seperate_AA <- "#" %in% tAI$trans
# if (start_codon_is_seperate_AA) subseq(seqs, 1, 1) <- "#"
# dt_range[, AA := ""]
# seqs <- strsplit(as.character(unlist(seqs, use.names = FALSE)), split = "")
# dt_range[, genes := groupings(tile)]
# dt_range[, position := seq_len(.N), by = genes]
# tAI_short <- tAI[, c("trans", "usage_score")]
# tAI_merged <- data.table::merge.data.table(data.table(trans = unlist(seqs, use.names = FALSE)),
#                                            tAI_short, by = "trans", sort = FALSE)
# dt_range[(position %% 3) == 1, proportion := tAI_merged$usage_score]
# dt_range[(position %% 3) == 2, proportion := tAI_merged$usage_score]
# dt_range[(position %% 3) == 0, proportion := tAI_merged$usage_score]
# 
# values <- dt_range$proportion
# values <- values * c(5,2,1)
# cov <- data.table(values, frame = rep(c(1,2,3), length.out = length(values)))
# ggplot(cov, aes(x = seq(length(values)), y = values, fill = frame)) + 
#   geom_col() + theme_classic()
# 
# cov[, values := values * frollmean(cov$values, n = 30, fill = 50)]
# ggplot(cov, aes(x = seq(length(values)), y = values, fill = frame)) + 
#   geom_col() + theme_classic()
# cov_real <- coveragePerTiling(cds[which.max(cds_hits)], reads, as.data.table = T, is.sorted = T, withFrames = T)
# ggplot(cov_real, aes(x = position, y = count, fill = frame)) + 
#   geom_col() + theme_classic()
# cov_real_cutoff <- copy(cov_real)
# cov_real_cutoff[count > mean(count)*2, count := frollmean(count, n = 15, fill = 50)]
# ggplot(cov_real_cutoff, aes(x = position, y = count, fill = frame)) + 
#   geom_col() + theme_classic()
# cov_real[, mean_cov := frollmean(cov$values, n = 15, fill = 50)]
# ggplot(cov_real, aes(x = position, y = mean_cov, fill = frame)) + 
#   geom_col() + theme_classic()
# peaks_index <- order(values, decreasing = T)[1:10]
# peaks_index <- 
# 
# abs(cos(seq(0, pi*(x/30), pi/30)))


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Original cor plot
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Real vs sim cor
# method <- c("pearson", "spearman")[1]
# dt_all_sim <- copy(dt_all_list_named["default"][[1]])
# dt_all_sim_f <- copy(dt_all_list_named["optimized"][[1]])
# combined_res_subset <- coverage_cor(dt_all_real, dt_all_sim, dt_all_sim_f,
#                                     method = method,
#                                     na.rm = T, rm.selfmatch = T)
# combined_res_subset <- combined_res_subset[(Var1_type != Var2_type & Var1_index == Var2_index & Var1_type == "R") |  (Var1_type == "R" & Var2_type == "R"),]
# table(combined_res_subset$Comparison)
# nt_cor_plot <- ggplot(combined_res_subset, aes(x = Cor, fill = Comparison)) +
#   geom_density(alpha = 0.75) +
#   ylab("Density") + 
#   xlab("Correlation per nt") + 
#   theme_classic()  +
#   scale_fill_brewer(palette = 16) + 
#   xlim(-0.1, 0.7) + 
#   theme(legend.position="top",
#         legend.title=element_blank(),
#         legend.text=element_text(size=8), 
#         legend.key.size = unit(0.3, "cm")); nt_cor_plot
# 
# 
# ggplot(combined_res_subset, aes(x = Comparison, y = Cor)) +
#   geom_violin(aes(fill = Comparison)) +
#   geom_boxplot(alpha = 0.9) + 
#   ylab("Correlation of coverage per nt") + 
#   theme_classic() +  theme(legend.position="none")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Distributions (PDF)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# dist_all_temp <- lapply(dt_all_list_named, function(x) suppressWarnings(melt(x)))
# dist_all <- rbindlist(dist_all_temp)
# dist_all[, id := rep(names(dist_all_temp), each = nrow(dist_all_temp[[1]]))]
# dist_all[value > 10, value := 10]
# ecdf_probs <- seq(0, 1, 0.005)
# dist_all <- dist_all[, quantile(value, probs = ecdf_probs), by = .(id, variable)]
# dist_all[, quantile := rep(ecdf_probs*100, length(dist_all_temp)*lengths(dt_all_list_named)[1])]
# colnames(dist_all)[3] <- "value"
# dist_all[, count_dist := value - dist_all[id == "real"]$value]
# dist_all_final <- dist_all[, .(count_dist = sum(abs(count_dist))), by = .(id, variable)]
# 
# levels <- dist_all_final[, .(median = median(count_dist)), by = id][order(median, decreasing = T),]
# dist_all_final[, id := factor(id, levels = levels$id, ordered = T)]
# count_density <- ggplot(dist_all_final, aes(y = count_dist, fill = id)) +
#   geom_boxplot(alpha = 0.7, size = 0.3) + theme_classic() + ylab("Total error of permilles\n (count: 1-10)"); count_density
# #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# # Correlation (nt level)
# #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# dt_all_comb_melt <- coverage_cor(dt_all_list_named, is.listed = T)
# dt_all_comb_melt <- dt_all_comb_melt[Var1_type != Var2_type & Var1_index == Var2_index & Var1_type == "R",]
# comparison_median <- dt_all_comb_melt[, .(median = median(Cor)), by = Var2_type][order(median, decreasing = T),]
# dt_all_comb_melt[, Var2_type := factor(Var2_type, levels = comparison_median$Var2_type, ordered = TRUE)]
# 
# nt_cor_plot_single <- ggplot(dt_all_comb_melt, aes(x = Var2_type, y = Cor, fill = Var2_type)) +
#   geom_boxplot(alpha = 0.35) +
#   ylab("Correlation per nt") + 
#   xlab("Comparison: Model vs Real") + 
#   geom_hline(yintercept = comparison_median[Var2_type == "N",]$median, linetype = "dashed", col = "gray4") + 
#   theme_classic()  +
#   scale_fill_viridis_d() +
#   theme(legend.position="none",
#         legend.title=element_blank(),
#         legend.text=element_text(size=8), 
#         legend.key.size = unit(0.3, "cm")); nt_cor_plot_single
# ggsave(file.path(test_dir, "SA_analysis.png"), nt_cor_plot_single, width = 6, heigh = 5, dpi = 600)

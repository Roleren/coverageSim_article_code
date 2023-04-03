# Simple walkthrough of steps done to create Ribo-seq
# Estimators to simulate replicates from real dataset
# mean, variance, kurtiosis, rnase bias (with auto correlation), 
# codon bias (with auto correlation), Ribo-seq lab protocol bias and sequencer machine bias
library(coverageSim);library(ggplot2)
test_dir <- "/export/valenfs/projects/Hakon/coverageSim/"

gene_size <- 60 # In codons
# Step 1 Ideal coverage
ideal_cov <- rep(c(1, 0, 0), gene_size)
barplot(ideal_cov, col = c("red", "green", "blue"))
# Step 2 Ideal RNAase inacurracy
RNAse_reach <- 3
alpha_roll_RNAse <- coverageSim:::autocor_window(ideal_cov, RNAse_reach, na.rm = TRUE)
barplot(alpha_roll_RNAse, col = c("red", "green", "blue"))

# Step 3 Codon
# codon_alphas <- c(1, sample(c(0.5, 0.2, 0.1, 0.01), size = gene_size - 1, replace = T,
#                             prob = c(0.1, 0.1, 0.1, 0.5)))
seq_biases <- load_seq_bias(bias = "all")[variable == "R7",]$alpha
codon_alphas <- c(seq_biases[1], sample(seq_biases[-1], size = gene_size - 1, replace = T))
barplot(codon_alphas)
quantile <- sample(c(0.85, 9, 9.5), 1, T, prob = c(0.25, 0.5, 0.25))/10
to_be_peaks <- codon_alphas >= quantile(codon_alphas, quantile)
barplot(to_be_peaks)
codon_alphas_scaled <- codon_alphas
codon_alphas_scaled[to_be_peaks] <- codon_alphas_scaled[to_be_peaks] *
  (codon_alphas_scaled[to_be_peaks] / max(codon_alphas)) * 3 * sample(seq(0.5, 3, by = 0.1), size = sum(to_be_peaks), replace = T)
barplot(codon_alphas_scaled)
# Step 4 Add codon autocorrelation
codon_auto_cor_reach <- 30
codon_ac_alphas <- coverageSim:::autocor_window(codon_alphas_scaled, codon_auto_cor_reach, 
                                  padding.rm = T, penalty = 1.25)
barplot(codon_ac_alphas)
acf(codon_ac_alphas)
codon_ac_rnase_alphas_ideal <- rep(codon_ac_alphas, each = 3)
codon_ac_rnase_alphas_ideal[seq(codon_ac_rnase_alphas_ideal) %% 3 %in% c(2,0)] = 0
codon_ac_rnase_alphas_ideal <- c(rep(0, RNAse_reach*2), codon_ac_rnase_alphas_ideal,
                                 rep(0, RNAse_reach*2))
barplot(codon_ac_rnase_alphas_ideal, col = c("red", "green", "blue"))
# alpha_roll_final <- alpha_roll_final * alpha_roll_RNAse
rnase_W <- c(0.5,1,2,10,2,1,0.5)*2
rnase_W <- c(0.07023652/3, 0.21070956/2, 0.19992970/2, 0.78399483,
            0.19992970, 0.21070956, 0.07023652)
final_dist <- frollapply(codon_ac_rnase_alphas_ideal, 
                                       FUN = function(x) sum(x*rnase_W),
                                       n = c(7), align = "center", fill = NA)
final_dist <- final_dist[!is.na(final_dist)]
final_dist <- final_dist * (mean(codon_alphas) / mean(final_dist))
final_dist[final_dist == 0] <- 1e-24
barplot(final_dist[1:100], col = c("red", "green", "blue"))
# High coverage, high average alphas
sample_Ribo_1 <- extraDistr::rdirmnom(n = 1, size = 1000, alpha = final_dist*50)
sample_Ribo_1 <- as.vector(t(sample_Ribo_1)); sample_Ribo_1
barplot(sample_Ribo_1, col = c("red", "green", "blue"))
# High coverage, low average alphas
sample_Ribo_3 <- extraDistr::rdirmnom(n = 1, size = 1000, alpha = final_dist)
sample_Ribo_3 <- as.vector(t(sample_Ribo_3)); sample_Ribo_3
barplot(sample_Ribo_3, col = c("red", "green", "blue"))
coverageSim:::frame_usage(data.table(sample_Ribo_3), normalize_to = 100)
coverageSim:::frame_usage(data.table(sample_Ribo_1, sample_Ribo_3), normalize_to = 100)
# Low coverage, high average alphas
sample_Ribo_2 <- extraDistr::rdirmnom(n = 1, size = 50, alpha = final_dist*50)
sample_Ribo_2 <- as.vector(t(sample_Ribo_2)); sample_Ribo_2
barplot(sample_Ribo_2, col = c("red", "green", "blue"))
# Low coverage, low average alphas
sample_Ribo_4 <- extraDistr::rdirmnom(n = 1, size = 50, alpha = final_dist/20)
sample_Ribo_4 <- as.vector(t(sample_Ribo_4)); sample_Ribo_4
barplot(sample_Ribo_4, col = c("red", "green", "blue"))


all_steps_list <- list("Uniform triplet translocation" = ideal_cov, "RNase bias" = alpha_roll_RNAse, "Codon bias" = codon_alphas, "Codon auto correlation" = codon_ac_alphas,
                       "Probability density function" = final_dist, "Sample 1 (Ribo-seq)" = sample_Ribo_1, "Sample 2 (Ribo-seq)" = sample_Ribo_3)
all_covs <- lapply(seq_along(all_steps_list), function(i) {
  x <- all_steps_list[[i]]; name <- names(all_steps_list[i])
  data.table(count = x, pos = seq(x), id = name, frame = rep(as.factor(seq(3)), length.out = length(x)))
})

all_cov_dt <- rbindlist(all_covs); all_cov_dt[, id := factor(id, levels = unique(id), ordered = TRUE)]
all_cov_dt <- all_cov_dt[-1, ]
all_cov_dt[!(id %in% c("Codon bias", "Codon auto correlation")), pos := pos - 4]
all_cov_dt[(id %in% c("Codon bias", "Codon auto correlation")), pos := pos - 1]
all_cov_dt[id %in% c("Codon bias", "Codon auto correlation"), frame := NA]
# meta_cov <- ggplot(all_cov_dt, aes(x = pos, y = count, fill = frame)) +
#   geom_col() + theme_classic() + ylab("Coverage") + xlab("Position (nt)") +
#   facet_wrap(~ id, ncol = 1, scales = "free") +
#   theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank()); meta_cov
# meta_cov

meta_cov1 <- ggplot(all_cov_dt[id %in% unique(all_cov_dt$id)[c(3,4)]], aes(x = pos, y = count)) + 
  geom_col() + theme_classic() + ylab("Coverage") + xlab("Position (codon)") +
  facet_wrap(~ id, ncol = 1, scales = "free") + 
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), strip.background = element_blank());
meta_cov2 <- ggplot(all_cov_dt[id %in% unique(all_cov_dt$id)[c(1,2)]], aes(x = pos, y = count, fill = frame)) + 
  geom_col() + theme_classic() + ylab("") + xlab("Position (nt)") +
  facet_wrap(~ id, ncol = 1, scales = "free") +
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), legend.position = "none", strip.background = element_blank());
meta_cov3 <- ggplot(all_cov_dt[id %in% unique(all_cov_dt$id)[c(5)]], aes(x = pos, y = count, fill = frame)) + 
  geom_col() + theme_classic() + ylab("Coverage") + xlab("") +
  facet_wrap(~ id, ncol = 1, scales = "free") + guides(fill=guide_legend(title="Frame")) +
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), 
        legend.position = "bottom", strip.background = element_blank(),
        legend.title=element_text(size=9, margin = margin(r=0,l=0,t=0,b=0)),
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.3, "cm"));
meta_cov4 <- ggplot(all_cov_dt[id %in% unique(all_cov_dt$id)[c(6)]], aes(x = pos, y = count, fill = frame)) + 
  geom_col() + theme_classic() + ylab("Coverage") + xlab("Position (nt)") +
  facet_wrap(~ id, ncol = 2, scales = "free") +
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), legend.position = "none", strip.background = element_blank());
meta_cov5 <- ggplot(all_cov_dt[id %in% unique(all_cov_dt$id)[c(7)]], aes(x = pos, y = count, fill = frame)) + 
  geom_col() + theme_classic() + ylab("") + xlab("Position (nt)") +
  facet_wrap(~ id, ncol = 2, scales = "free") +
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(), legend.position = "none", strip.background = element_blank());
# grid <- gridExtra::arrangeGrob(meta_cov1, meta_cov2, meta_cov3, meta_cov4, layout_matrix = matrix(c(1,2,1,2,3,3,3,3,4,4,4,4), nrow = 6, byrow = T))
# plot(grid)
grid <- cowplot::plot_grid(cowplot::plot_grid(meta_cov1, meta_cov2, labels = "AUTO", label_size = 12),
                           cowplot::plot_grid(meta_cov3, labels = "C", label_size = 12),
                           cowplot::plot_grid(meta_cov4, meta_cov5, labels = c("D", "E"), label_size = 12),
                           ncol = 1, rel_heights = c(1.5,1.5,1))
plot(grid)
ggsave(file.path(test_dir, "sim_cov_steps.png"), grid, dpi = 600, width = 6, height = 5.7)
ggsave(file.path(test_dir, "sim_cov_steps.pdf"), grid, width = 6, height = 5.7)

list_dt <- list(alpha = data.table(Cor = alpha_roll_final, distance = seq(length(alpha_roll_final))))
alpha_roll_dt <- ac_list_to_dt(list_dt)
ac_boxplot(alpha_roll_dt, autocor_name = "nt", breaks.by = 27)
dt_FFT <- as.data.table(spec.pgram(x = sample_Ribo_1, plot = F)[c(1,2)])
dt_FFT <- as.data.table(spec.pgram(x = codon_ac_alphas, plot = F)[c(1,2)])
ggplot(dt_FFT, aes(1 / freq, spec)) + geom_line() + theme_classic() + 
   xlim(c(1, 10)) + geom_vline(xintercept = 3, colour = "gray", size = 0.1, linetype = "dashed")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Resampling bias
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
resamples <- lapply(seq(5), function(x) as.vector(t(extraDistr::rdirmnom(n = 1, size = 1000, alpha = final_dist*5))))
resamples <- data.table(count = c(final_dist, unlist((resamples))),
                        sample = c(rep("theoretical (pdf)", length.out = length(final_dist)),
                          paste("sim", rep(seq(5), each = length(final_dist)))))
resamples[, sample := factor(sample, levels = unique(sample), ordered = T)]
resamples[, position := seq(.N) - 4, by = sample]
resamples[, frame := position %% 3]
resamples[, frame := as.factor(frame)]
resample_cov <- ggplot(resamples, aes(x = position, y = count, fill = frame)) + 
  geom_col() + theme_classic() + ylab("Coverage") + xlab("Position (nt)") +
  facet_wrap(~ sample, ncol = 1, scales = "free") +
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank()); resample_cov

ggsave(file.path(test_dir, "resample.png"), resample_cov, dpi = 300, width = 6, height = 8)

# codon_extreme <- max(codon_alphas) / median(codon_alphas)
# codon_variance <- sd(codon_alphas) / median(codon_alphas)
# codon_alphas <- codon_alphas**0.5
# quantile <- sample(seq(6, 9), 1, T, prob = c(0.2, 0.3, 0.4, 0.1))/10
# scaler <- 30*(codon_variance / codon_extreme); scaler
# codon_alphas_scaled <- (codon_alphas/(quantile(codon_alphas, 0.8)))**scaler

#TODO: 
# 1. Analyse variance to justify bumps
# 2. Add MCC of simple features

library(ORFik); library(data.table); library(RiboCrypt); library(ggplot2)
source("/export/valenfs/projects/Hakon/coverageSim/Create_simulated_sample_data.R")

# Define output dirs
test_dir <- "/export/valenfs/projects/Hakon/coverageSim/"
data_output_dir <- file.path(test_dir, "data_output")
all_real_result_objects_path <- file.path(data_output_dir, "all_saved_objects.rdata")
# Load required data
if (file.exists(all_real_result_objects_path)) { 
  load(all_real_result_objects_path)
  merged_gr <- import.ofst(merged_gr_path, seqinfo = seqinfo(cds))
  merged_gr_covRle <- fimport(merged_gr_covRle_path)
} else stop("You must have run the real data validation and saved in to correct path!")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Shape analysis validation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create shape model only
fft_dt <- data.table()
cov_dt <- data.table()
genes <- 1000
length <- sample(seq(29, 1000, by = 12), size = genes, replace = TRUE)
#A <- sample(seq(0, 14), size = 50, replace = TRUE, prob = seq.int(1.4, 0, -0.1)^2)
b <- sample(seq(3, 15, by = 0.1), size = genes, replace = TRUE)
pi_s <- sample(pi*seq(5), prob = 1/seq(5), size = genes, replace = TRUE)
shapy <- shapes(3, v = 1)
for (i in seq(genes)) {
  x <- seq(0, pi_s[i], pi_s[i]/5)
  values <- eval(shapy)
  values <- values*c(5,2,1)
  cov <- data.table(pos = seq(length(values)), values, frame = factor(rep(c(5,2,1), length.out = length(values))), genes = i)
  cov_dt <- rbindlist(list(cov_dt, cov))
  spec <- spec.pgram(values, plot = FALSE)
  
  fft_dt <- rbindlist(list(fft_dt, data.table(read_length = i, amplitude = spec$spec, periods = 1 / spec$freq)))
}
max_period_is_3 <- fft_dt[fft_dt[, .I[amplitude == max(amplitude)], by = read_length]$V1,]$periods
table(max_period_is_3 > 2.9 & max_period_is_3 < 3.1)
sd(cov_dt$values)
fft_dt2 <- fft_dt[, .(amplitude = sum(amplitude, na.rm = T), read_length = 1), by = .(periods)]
gg_periods <- ggplot(fft_dt2, aes(x = periods, y = amplitude)) +
  geom_line() +
  scale_x_continuous(limits = c(0,50), breaks = c(3, 20, 40)) +
  theme_classic(); gg_periods

#ORFik::ribo_fft_plot(fft_dt, period_window = c(0, 50))
# genome <- simGenome(out_dir = ,n = 500, max_uorfs = 0, genome_name = "test_shapes")

cov_final <- cov_dt[, .(count = sum(values) / .N), by = .(pos)]
cov_final[, frame := as.factor(rep(c(0,1,2), length.out = nrow(cov_final)))]
plot(cov_final$count)
meta_cov <- ggplot(cov_final, aes(x = pos, y = count, fill = frame)) + 
  geom_col() + theme_classic() + xlab("Coverage") + ylab("Position (nt)")
meta_plots <- gridExtra::arrangeGrob(meta_cov, gg_periods, nrow = 2); plot(meta_plots)

ggplot(cov_dt[genes == 19,], aes(x = pos, y = values, fill = frame)) + 
  geom_col() + theme_classic()
gg_coverage_shape_all <- ggplot(cov_dt[genes %in% 1:30], aes(x = pos, y = values, fill = frame)) + 
  geom_col() + theme_classic() + xlab("Coverage") + ylab("Position (nt)") + facet_wrap(~ genes, scales = "free_x") + scale_x_continuous(breaks = 300)

ggsave(file.path(test_dir, "Suppl_all_fft_sim.png"), meta_plots, width = 5, height = 3, dpi = 400)
ggsave(file.path(test_dir, "Suppl_all_shapes_coverage_sim.png"), gg_coverage_shape_all,  width = 6, height = 10, dpi = 400)

values <- rep(abs(cos(seq(1, pi*(x/30), pi/30))), length.out = x)
waves <- 40
scaler <- sample(20:100, size = waves, replace = T); scaler2 <- sample(seq(0, 1, 0.1), size = waves, replace = T); const <- 0.1
for (i in waves) {values <- values * rep(abs(cos(seq(scaler2[i], pi*(x/scaler[i]), pi/scaler[i]))), length.out = x) }
values <- values + const

plot(values)
values <- values*c(5,2,1)
values <- rmultinom(1, size = 1000, prob = cov_dt[genes == 1,]$values)
values <- dt_all_real[seq(widthPerGroup(cds[1])),]$R4
values <- dt_all_real[seq(widthPerGroup(cds[1])),]$R5
values <- rnbinom(c(1305), mu = 13, size = 0.1)*c(5,2,1)
values2 <- rnbinom(4*3, mu = 13, size = 0.1)*c(5,2,1)
cov_test <- data.table(count = values, 
                       position = seq(length(values)), library = "1",
                       frame = factor(rep(c(1,2,3), length.out = length(values))))
cov_test2 <- data.table(count = values2, 
                        position = seq(length(values2)), library = "2",
                        frame = factor(rep(c(1,2,3), length.out = length(values2))))
cov_test <- rbindlist(list(cov_test, cov_test2))
ggplot(cov_test, aes(x = position, y = count, fill = frame)) + 
  geom_col() + theme_classic() + facet_wrap(~ library, ncol = 1)
which.max(counts_cds)

spec <- spec.pgram(values, plot = F)
fft_dt <- data.table(read_length = i, amplitude = spec$spec, periods = 1 / spec$freq)
ggplot(fft_dt, aes(x = periods, y = amplitude)) +
  geom_line() +
  scale_x_continuous(limits = c(0,50)) +
  theme_classic()

#quote(abs(cos(seq(0, pi*(x/30), pi/30)) +)

# Using autocorrelation by cosine
lengths <- 30
b <- sample(seq(3,15, by = 0.1), size = length(lengths), replace = TRUE)
pi_s <- sample(pi*seq(5), prob = 1/seq(5), size = length(lengths), replace = TRUE)
periods <-  lapply(seq_along(lengths), function(i) seq(0, pi_s[i], pi_s[i]/lengths[i]))
bump_shape <- lapply(seq_along(lengths), 
                     function(i, l = lengths[i], x = periods[[i]]) 
                       rep_len(eval(shapes(3)), l))
barplot(unlist(bump_shape), col = c("red", "green", "blue"))
barplot(aaa_vec*unlist(bump_shape), col = c("red", "green", "blue"))
# Using autocorrelation by rolling window

###
# aaa_dt <- data.table(count = aaa_vec); aaa_dt$count[1] <- aaa_dt$count[1]*3
# aaa_dt_copy <- copy(aaa_dt)

aaa_dt[rep(c(F, F, T), length.out = nrow(aaa_dt)), count := frollmean(count, 3, align = "left", )]
aaa_dt[rep(c(F, T, F), length.out = nrow(aaa_dt)), count := frollmean(count, 3, align = "left", )]
aaa_dt[rep(c(T, F, F), length.out = nrow(aaa_dt)), count := frollmean(count, 3, align = "left", )]
aaa_dt[is.na(count), count := aaa_dt_copy$count[is.na(aaa_dt$count)]]
auto_cor_test <- auto_correlation(aaa_dt, 10, by.codon = T)
  
barplot(aaa_dt_copy$count, col = c("red", "green", "blue"))
barplot(aaa_dt$count, col = c("red", "green", "blue"))
# Old uORF check
if (is.null(names(region_counts))) names(region_counts) <- rownames(assay)
uorf_tab <- table(txNames(uorf_ranges))
if (uorf_prop_mode == "character") {
  a <- if (uorf_prop_within_gene == "uniform") {
    lapply(names(uorf_tab), function(x) rmultinom(1, region_counts[x], rep(1, uorf_tab[x])))
  } else {
    stop("Not implemented yet")
    lapply(names(uorf_tab), function(x) rmultinom(1, region_counts[x], rep(1, uorf_tab[x])))
  }
} else a <-  lapply(names(uorf_tab), function(x) rmultinom(1, region_counts[x],
                                                           uorf_prop_within_gene[names(uorf_prop_within_gene) %in% x]))
a <- unlist(a); names(a) <- txNames(uorf_ranges)
region_counts <- a
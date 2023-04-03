library(data.table)
res <- data.table()
mu <- rep(c(2,4,6,10), each = 2)
sd <- c(2,3,2,1,2,1,2,4)
for (i in seq_along(mu)) {
  a <- DESeq2::makeExampleDESeqDataSet(n = 50000, interceptMean = mu[i], interceptSD = sd[i])
  res <- rbindlist(list(res, data.table(mu = mu[i], sd = sd[i], 
                                        mean = round(mean(assay(a)), 0), median = round(median(assay(a))),
                                        zero_count = round((sum(assay(a)[,1] == 0) / 50000)*100, 0),
                                        max = max(assay(a)))))
}
fwrite(res, file = "~/Desktop/gene_count_modeling_res.csv")

# Human genome, index (STAR), download data, Align, QC and p-shifting

library(ORFik)
ORFik::config() # Check relative save locations (General)
path_config <- ORFik::config.exper(experiment = "coverageSim_results",
  assembly = "Homo_sapiens_ensembl_GRCh38", "") # Set Relative save locations (Specific)

# Annotation and index
annotation <- getGenomeAndAnnotation("Homo sapiens", path_config["ref"], optimize = TRUE)
index <- STAR.index(annotation)

# Download data
runs <- fread("https://github.com/Roleren/coverageSim_article_code/tree/main/Study_info.csv")
runs <- runs$`Run ID`
ORFik::download.SRA(runs, outdir = path_config["fastq "])

# Align data
STAR.align.folder(path_config["fastq"], path_config["bam "], index)

# Create ORFik experiment
create.experiment(path_config["bam"], "coverageSim_human_data",
                  txdb = paste0(annotation["gtf", ".db"]),
                  fa = annotation["genome"], libtype = "RFP",
                  fraction = paste0("R", seq(10)))
df <- read.experiment("coverageSim_human_data")
# Quality control
ORFikQC(df)
# Shift
shiftFootprintsByExperiment(df)




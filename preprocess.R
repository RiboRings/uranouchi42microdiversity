library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(reshape)
library(ComplexHeatmap)
library(tibble)
library(plotly)
library(iterators)
library(MatrixGenerics)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(doParallel)
library(circlize)

checkm <- read_csv("data/genomeInfo.csv") %>%
  mutate(genome = gsub(".fa", "", genome))
coverm <- read_tsv("data/assays.tsv")
names(coverm) <- gsub(" ", "_",
                 gsub(" \\(%\\)", "",
                 gsub("reference.fa.(UU.*)_1.paired.fastq.gz", "\\1",
                      names(coverm))))
ar_gtdbtk <- read_tsv("data/gtdbtk.ar122.summary.tsv")
bac_gtdbtk <- read_tsv("data/gtdbtk.bac120.summary.tsv")
gtdbtk <- rbind(ar_gtdbtk, bac_gtdbtk) %>%
  select(user_genome, classification) %>%
  separate(classification,
           c("Domain",
             "Phylum",
             "Class",
             "Order",
             "Family",
             "Genus",
             "Species"),
           ";")

diversity_files <- list.files(pattern = 'UU.*_microdiversity_summary.csv',
                              path = "data")

diversity_file_list <- mclapply(paste0("data/", diversity_files),
                                sample_loader,
                                mc.cores = 3)

micro_df <- merge_all(diversity_file_list) %>%
  na.omit() %>%
  dplyr::rename(nucdiv = nucl_diversity) %>%
  mutate(sample = gsub("_microdiversity_summary", "", sample)) %>%
  group_by(genome)

diversity_stats <- micro_df %>%
  summarise(MeanCov = mean(coverage),
            MinCov = min(coverage),
            MaxCov = max(coverage))

diversity_df <- micro_df %>%
  select(genome,
         nucdiv,
         DiSiperMbp,
         NonsynonimousFraction,
         coverage,
         sample) %>%
  pivot_wider(names_from = sample,
              values_from = c(nucdiv, DiSiperMbp, NonsynonimousFraction, coverage),
              values_fill = 0) %>%
  left_join(diversity_stats)

df <- checkm %>%
  right_join(gtdbtk, by = c("genome" = "user_genome")) %>%
  left_join(coverm, by = c("genome" = "Genome")) %>%
  mutate(genome = gsub("matabat2bin.", "", genome)) %>%
  right_join(diversity_df)

names(df) <- gsub("(UU\\w)(\\d{6})", "\\2\\1", names(df))
df <- df[ , order(names(df))]

rpkm_mat <- df %>%
  select(ends_with("RPKM")) %>%
  as.matrix()
rpkm_mat[rpkm_mat == 0] <- min(rpkm_mat[rpkm_mat > 0]) / 2
rpkm_shift_mat <- matrix(nrow = dim(rpkm_mat)[[1]], ncol = dim(rpkm_mat)[[2]] - 1)

for (i in seq(1, dim(rpkm_mat)[[2]] - 1)) {
  
  rpkm_shift_mat[ , i] <- map2_dbl(rpkm_mat[ , i + 1], rpkm_mat[ , i], `/`)
  
}

colnames(rpkm_shift_mat) <- sapply(seq(1, dim(rpkm_shift_mat)[[2]]), function(x) paste0("Ratio_", x))
df <- df %>%
  cbind(as.data.frame(rpkm_shift_mat))

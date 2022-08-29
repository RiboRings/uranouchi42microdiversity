library(readr)
library(reshape)
library(doParallel)

file_dir <- "data/"
file_vec <- list.files(path = file_dir, pattern = "*.contigs.csv")
file_list <- list()

for (n in 1:length(file_vec)){
  file_list[[n]] <- read_csv(paste0(file_dir, file_vec[n]))
}

all_df <- merge_all(file_list) %>%
  group_by(genome) %>%
  summarise(Scaffolds = mean(true_scaffolds)) %>%
  mutate(genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)))

df <- df %>% left_join(all_df)

write_csv2(df, "data/mags.csv")


file_merger <- function(pattern, dir = "data/", method = read_csv2) {
  
  file_vec <- list.files(path = dir, pattern = pattern)
  file_list <- list()
  
  file_list <- foreach (i = file_vec) {
    method(paste0(file_dir, i)) %>%
      na.omit() %>%
      filter(genome %in% selected_genomes) %>%
      mutate(sample = i,
             genome = gsub(".fna", "", gsub("metabat2bin_", "", genome)),
             position = as.factor(position))
  }
  
  full_df <- merge_all(file_list) %>%
    group_by(genome) %>%
    summarise(Scaffolds = mean(true_scaffolds)) %>%
    mutate(genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)))
  
  
}


snvs_files <- list.files(pattern = 'UU.*_snvs_summary.csv',
                         path = "data")
snvs_list <- mclapply(paste0("data/", snvs_files),
                         read_csv2)
snvs_file_list <- list()

sample_loader <- function(sample_name) {
  
  read_csv2(sample_name) %>%
    mutate(sample = gsub(".csv", "", sample_name),
           genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)))
  
}

micro_df <- merge_all(diversity_file_list) %>%
  na.omit() %>%
  dplyr::rename(nucdiv = nucl_diversity) %>%
  group_by(genome)

gc_content <- read_tsv("data/gc.tsv", col_names = c("genome", "GC"))
gc_content <- gc_content %>%
  mutate(genome = gsub("matabat2bin.", "", gsub("_k141_\\d*", "", gsub("_*UU", "UU", genome))))
gc_content <- gc_content %>%
  group_by(genome) %>%
  summarise(GC = mean(GC))

base_content <- read_tsv("data/base_content.tsv", col_names = c("genome", "BaseContent"))
base_content <- base_content %>%
  mutate(genome = gsub("matabat2bin.", "", gsub("_k141_\\d*", "", gsub("_*UU", "UU", genome))))
base_content <- base_content %>%
  group_by(genome) %>%
  summarise(BaseContent = mean(BaseContent))

base_count <- read_tsv("data/base_count.tsv", col_names = c("genome", "BaseCount"))
base_count <- base_count %>%
  mutate(genome = gsub("matabat2bin.", "", gsub("_k141_\\d*", "", gsub("_*UU", "UU", genome))))
base_count <- base_count %>%
  group_by(genome) %>%
  summarise(BaseCount = mean(BaseCount))

seqkit_df <- gc_content %>%
  full_join(base_content) %>%
  full_join(base_count)

seqkit_df <- filtered_df %>%
  left_join(seqkit_df) %>%
  select(genome, GC, BaseContent, BaseCount)

write_csv2(seqkit_df, "data/seqkit.csv")
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
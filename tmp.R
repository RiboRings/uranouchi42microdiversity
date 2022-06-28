library(readr)

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

exp_df <- gene_df %>%
  mutate(genome = gsub("_k141.*", "", gsub("matabat2bin.", "", gene))) %>%
  group_by(genome, gene) %>%
  summarise(gene_length = mean(gene_length))

exp_df2 <- exp_df %>%
  group_by(genome) %>%
  summarise(GeneCount = n(),
            CodingLength = sum(gene_length))

exp_df3 <- exp_df2 %>%
  right_join(filtered_df) %>%
  mutate(CodingDensity = CodingLength / length * 100) %>%
  select(genome, GeneCount, CodingLength, CodingDensity)

write.csv2(exp_df3, "data/coding_density.csv")

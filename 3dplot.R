filtered_df$NucDivMean <- filtered_df %>%
  select(starts_with("nucdiv")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

plot_ly(data = filtered_df,
        x = ~completeness,
        y = ~NucDivMean,
        z = ~contamination,
        type = "scatter3d",
        mode = "markers",
        size = ~AbundMax) %>%
  layout(scene = list(yaxis = list(title = 'Mean ND')))

p1 <- ggplot(filtered_df, aes(x = NucDivMean, y = N50)) +
  geom_point(size = 1) +
  scale_x_continuous(limits = c(0, 0.06),
                     breaks = seq(0, 0.06, by = 0.01)) +
  scale_y_continuous(limits = c(0, 10^6)) +
  labs(x = "Mean Nucleotide Diversity") +
  theme_classic()

tmp <- read_csv2("data/tmp.txt", col_names = FALSE)
tmp1 <- tmp[seq(1, nrow(tmp), by = 2), ]
tmp2 <- tmp[seq(2, nrow(tmp), by = 2), ]
contigs <- data.frame(genome = tmp1$X1, contigs = tmp2$X1)
contigs <- contigs %>%
  mutate(genome = gsub(".fa", "", gsub("matabat2bin.", "", gsub("dereplicated_genomes/", "", genome))),
         contigs = as.numeric(contigs))
filtered_df <- filtered_df %>%
  left_join(contigs)

filtered_df <- filtered_df %>%
  mutate(ContigLength = contigs / length * 1000000)

p2 <- ggplot(filtered_df, aes(x = NucDivMean, y = ContigLength)) +
  geom_point(size = 1) +
  scale_x_continuous(limits = c(0, 0.06),
                     breaks = seq(0, 0.06, by = 0.01)) +
  scale_y_log10() +
  labs(x = "Mean Nucleotide Diversity",
       y = "Contigs per million Base Pairs") +
  theme_classic()

p <- p1 + p2 +
  plot_annotation(tag_levels = 'A')

ggsave("quality.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 5)

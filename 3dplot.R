filtered_df$NucDivMean <- filtered_df %>%
  select(starts_with("nucdiv")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

filtered_df$NucDivMax <- filtered_df %>%
  select(starts_with("nucdiv")) %>%
  as.matrix() %>%
  rowMaxs(na.rm = TRUE)

plot_ly(data = filtered_df,
        x = ~completeness,
        y = ~NucDivMean,
        z = ~contamination,
        type = "scatter3d",
        mode = "markers",
        size = ~AbundMax) %>%
  layout(scene = list(yaxis = list(title = 'Mean ND')))

p1 <- ggplot(filtered_df, aes(x = NucDivMean, y = N50, colour = MaxCov)) +
  geom_point(size = 1) +
  scale_x_continuous(limits = c(0, 0.06),
                     breaks = seq(0, 0.06, by = 0.01)) +
  scale_y_continuous(limits = c(0, 10^6)) +
  labs(x = "Mean Nucleotide Diversity",
       colour = "Max Coverage") +
  scale_colour_gradientn(colours = c("blue", "yellow", "red"),
                         limits = c(0, max(filtered_df$MaxCov)),
                         breaks = seq(0, 1250, by = 250)) +
  theme_classic() +
  stat_cor(method = "pearson",
           label.x.npc = "centre",
           label.y.npc = "top")

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

p2 <- ggplot(filtered_df, aes(x = NucDivMean, y = ContigLength, colour = MaxCov)) +
  geom_point(size = 1) +
  scale_x_continuous(limits = c(0, 0.06),
                     breaks = seq(0, 0.06, by = 0.01)) +
  scale_y_log10() +
  labs(x = "Mean Nucleotide Diversity",
       y = "Number of Contigs per Mbp",
       colour = "Max Coverage") +
  scale_colour_gradientn(colours = c("blue", "yellow", "red"),
                         limits = c(0, max(filtered_df$MaxCov)),
                         breaks = seq(0, 1250, by = 250)) +
  theme_classic() +
  stat_cor(method = "pearson",
           label.x.npc = "centre",
           label.y.npc = "bottom")

p3 <- ggplot(filtered_df, aes(x = MeanCov, y = N50, colour = NucDivMax)) +
  geom_point(size = 1) +
  scale_x_continuous(limits = c(0, 450),
                     breaks = seq(0, 450, by = 50)) +
  scale_y_continuous(limits = c(0, 10^6)) +
  labs(x = "Mean Coverage",
       y = "N50",
       colour = "Max ND") +
  scale_colour_gradientn(colours = c("blue", "yellow", "red"),
                         limits = c(0, max(filtered_df$NucDivMax)),
                         breaks = c(seq(0, 0.06, by = 0.02), 0.075)) +
  theme_classic() +
  stat_cor(method = "pearson",
           label.x.npc = "centre",
           label.y.npc = "top")

p4 <- ggplot(filtered_df, aes(x = MeanCov, y = ContigLength, colour = NucDivMax)) +
  geom_point(size = 1) +
  scale_x_continuous(limits = c(0, 450),
                     breaks = seq(0, 450, by = 50)) +
  scale_y_log10() +
  labs(x = "Mean Coverage",
       y = "Number of Contigs per Mbp",
       colour = "Max ND") +
  scale_colour_gradientn(colours = c("blue", "yellow", "red"),
                         limits = c(0, max(filtered_df$NucDivMax)),
                         breaks = c(seq(0, 0.06, by = 0.02), 0.075)) +
  theme_classic() +
  stat_cor(method = "pearson",
           label.x.npc = "centre",
           label.y.npc = "bottom")

p <- ((p1 | p2) + plot_layout(guides = "collect")) / 
        ((p3 | p4) + plot_layout(guides = "collect")) +
  plot_annotation(tag_levels = 'A')

ggsave("quality.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 15,
       height = 10)

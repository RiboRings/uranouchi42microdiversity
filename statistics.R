mean_nucdiv <- df %>%
  select(starts_with("nucdiv")) %>%
  as.matrix() %>%
  rowMeans()

# mean GC content, stdev GC content, distribution
# mean completeness, stdev completeness, distribution
# mean contamination, stdev contamination, distribution
# mean strain heterogeneity, stdev strain heterogeneity, distribution
# mean length, stdev length, distribution
# mean contigs, stdev contigs, distribution
# mean N50, stdev N50, distribution

stats_df <- read_csv2("mags.csv")
seqkit_df <- read_csv2("data/seqkit.csv")
stats_df <- stats_df %>%
  left_join(seqkit_df)

p1 <- ggplot(stats_df, aes(x = length / 1000000)) +
  geom_histogram(bins = 100) +
  labs(x = "Genome Length (Mbp)",
       y = "Count") +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8)) +
  scale_y_continuous(limits = c(0, 70),
                     breaks = seq(0, 70, by = 10)) +
  theme_classic()

p2 <- ggplot(stats_df, aes(x = contigs)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) +
  scale_y_continuous(limits = c(0, 60),
                     breaks = seq(0, 60, by = 10)) +
  labs(x = "Contigs",
       y = "Count") +
  theme_classic()

p3 <- ggplot(stats_df, aes(x = N50)) +
  geom_histogram(bins = 100) +
  scale_x_log10(labels = scales::number_format()) +
  labs(y = "Count") +
  theme_classic()

pivoted_df <- stats_df %>%
  pivot_longer(cols = c("completeness",
                        "contamination",
                        "strain_heterogeneity",
                        "GC"),
               values_to = "checkm_percent",
               names_to = "checkm_parameter")

pivoted_df$checkm_parameter[pivoted_df$checkm_parameter == "strain_heterogeneity"] <- "Strain Heterogeneity"
pivoted_df$checkm_parameter[pivoted_df$checkm_parameter == "completeness"] <- "Completeness"
pivoted_df$checkm_parameter[pivoted_df$checkm_parameter == "contamination"] <- "Contamination"

p4 <- ggplot(pivoted_df, aes(x = checkm_parameter, y = checkm_percent)) +
  geom_violin() +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 10)) +
  labs(x = "MAG Quality Statistics",
       y = "Percentage (%)") +
  theme_classic()

p <- (p1 / p2 / p3 | p4) +
  plot_annotation(tag_levels = "A")

ggsave("statistics.pdf",
       plot = p,
       width = 20,
       height = 10,
       device = "pdf",
       path = "results")

stats_df <- data.frame("Completeness" = mean(df$completeness),
                       "Contamination" = mean(df$contamination),
                       "MeanCoverage" = mean(df$MeanCov),
                       "MinCoverage" = mean(df$MinCov),
                       "Breadth" = mean(rowMeans(ubiquity_mat)),
                       "DiSiperMbp" = mean(df$MeanDiSiperMbp),
                       "NonSynonymousFraction" = mean(df$MeanNonsynonimousFractionMean),
                       "NucDiv" = mean(mean_nucdiv),
                       "MaxRPKM" = mean(df$AbundMax)) %>%
  t()

colnames(stats_df) <- "U42"

write.csv2(stats_df, "results/statistics.csv")

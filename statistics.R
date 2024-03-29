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

stats_df <- read_csv2("data/mags.csv")
seqkit_df <- read_csv2("data/seqkit.csv")
coding_density_df <- read_csv2("data/coding_density.csv") %>%
  select(-`...1`)

stats_df <- stats_df %>%
  left_join(seqkit_df) %>%
  left_join(coding_density_df)

p1 <- ggplot(stats_df, aes(x = length / 1000000)) +
  geom_histogram(bins = 100,
                 colour = "black",
                 fill = "white") +
  labs(x = "Genome Length (Mbp)",
       y = "Number of MAGs") +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8)) +
  scale_y_continuous(limits = c(0, 70),
                     breaks = seq(0, 70, by = 10)) +
  theme_classic()

p2 <- ggplot(stats_df, aes(x = contigs)) +
  geom_histogram(bins = 100,
                 colour = "black",
                 fill = "white") +
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) +
  scale_y_continuous(limits = c(0, 60),
                     breaks = seq(0, 60, by = 10)) +
  labs(x = "Contigs",
       y = "Number of MAGs") +
  theme_classic()

p3 <- ggplot(stats_df, aes(x = N50)) +
  geom_histogram(bins = 100,
                 colour = "black",
                 fill = "white") +
  scale_x_log10(labels = scales::number_format()) +
  labs(y = "Number of MAGs") +
  theme_classic()

pivoted_df <- stats_df %>%
  pivot_longer(cols = c("completeness",
                        "contamination",
                        "CodingDensity",
                        "GC"),
               values_to = "checkm_percent",
               names_to = "checkm_parameter")

pivoted_df$checkm_parameter[pivoted_df$checkm_parameter == "CodingDensity"] <- "Coding Density"
pivoted_df$checkm_parameter[pivoted_df$checkm_parameter == "completeness"] <- "Completeness"
pivoted_df$checkm_parameter[pivoted_df$checkm_parameter == "contamination"] <- "Contamination"
pivoted_df$checkm_parameter[pivoted_df$checkm_parameter == "GC"] <- "GC Content"

pivoted_df$checkm_parameter <- factor(
  pivoted_df$checkm_parameter, 
  levels = c("Completeness", "Contamination", "GC Content", "Coding Density")
)

p4 <- ggplot(pivoted_df, aes(x = checkm_parameter, y = checkm_percent)) +
  geom_boxplot() +
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

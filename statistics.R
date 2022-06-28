mean_nucdiv <- df %>%
  select(starts_with("nucdiv")) %>%
  as.matrix() %>%
  rowMeans()

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

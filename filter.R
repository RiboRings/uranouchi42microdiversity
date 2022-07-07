#df <- read_csv2("data/mags.csv")

mat <- df %>%
  select(ends_with("RPKM")) %>%
  as.matrix()

df$AbundMean <- rowMeans(mat)
df$AbundMedian <- rowMedians(mat)
df$AbundMin <- rowMins(mat)
df$AbundMax <- rowMaxs(mat)
df$AbundRange <- rowMaxs(mat) - rowMins(mat)
df$AbundSd <- rowSds(mat)

abund_max_idx <- c()
for (i in 1:nrow(rpkm_mat)) abund_max_idx <- append(abund_max_idx, which(rpkm_mat[i, ] == rowMaxs(rpkm_mat, na.rm = TRUE)[[i]])[[1]])
df$AbundMaxIdx <- abund_max_idx

df$AbundRatioMean <- rowMeans(rpkm_shift_mat)
df$AbundRatioMedian <- rowMedians(rpkm_shift_mat)
df$AbundRatioMin <- rowMins(rpkm_shift_mat)
df$AbundRatioMax <- rowMaxs(rpkm_shift_mat)
df$AbundRatioRange <- rowMaxs(rpkm_shift_mat) - rowMins(rpkm_shift_mat)
df$AbundRatioSd <- rowSds(rpkm_shift_mat)

#p1 <- ggplot(df, aes(x = MaxCov, y = Scaffolds)) +
#  geom_point() +
#  geom_smooth() +
#  scale_x_log10() +
#  scale_y_log10() +
#  labs(x = "Maximum Coverage",
#       y = "Number of Contigs") +
#  theme_bw() +
#  theme(panel.grid = element_blank())
#
#p2 <- ggplot(df, aes(x = MaxCov, y = N50)) +
#  geom_point() +
#  geom_smooth() +
#  scale_x_log10() +
#  scale_y_log10() +
#  labs(x = "Maximum Coverage") +
#  theme_bw() +
#  theme(panel.grid = element_blank())

#p3 <- p1 + p2

ggsave("quality_check.pdf",
       plot = p3,
       device = "pdf",
       width = 14,
       height = 7,
       path = "results")

max_coverage_threshold <- 10

filtered_df <- df %>%
  filter(MaxCov > max_coverage_threshold) %>%
  arrange(desc(AbundMax))

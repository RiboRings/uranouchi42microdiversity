df$DiSiperMbpMean <- df %>%
  select(starts_with("DiSiperMbp")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

df$NonsynonimousFractionMean <- df %>%
  select(starts_with("NonsynonimousFraction")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

disi_df <- filtered_df %>%
  select(starts_with("DiSiperMbp_"))
names(disi_df) <- paste0("DiSiperMbp_", seq(1, ncol(disi_df)), ".")
disi_df$AbundMaxIdx <- filtered_df$AbundMaxIdx
disi_array <- sapply(iter(disi_df, by = "row"),
                     max_mapper, pattern = "DiSiperMbp")

nsf_df <- filtered_df %>%
  select(starts_with("NonsynonimousFraction_"))
names(nsf_df) <- paste0("NonsynonimousFraction_", seq(1, ncol(nsf_df)), ".")
nsf_df$AbundMaxIdx <- filtered_df$AbundMaxIdx
nsf_array <- sapply(iter(nsf_df, by = "row"),
                    max_mapper, pattern = "NonsynonimousFraction")

tax_df <- filtered_df %>%
  transmute(Tax = paste(Phylum, Order, Family, gsub("metabat2bin_", "", genome), sep = ";"))

microdiversity_df <- data.frame("DiSiperMbpAtAbundMax" = disi_array,
                                "NonsynonimousFractionAtAbundMax" = nsf_array,
                                "AbundMax" = filtered_df$AbundMax,
                                "Tax" = tax_df$Tax) %>%
  filter(NonsynonimousFractionAtAbundMax != 0) %>%
  mutate(Tax = ifelse(NonsynonimousFractionAtAbundMax > 0.6, Tax, ""))

p <- ggplot(microdiversity_df, aes(x = DiSiperMbpAtAbundMax,
                                   y = NonsynonimousFractionAtAbundMax,
                                   colour = AbundMax)) +
  geom_text_repel(aes(label = Tax),
                  size = 3,
                  colour = "red") +
  geom_point(size = 1) +
  labs(x = "Divergent Sites per Mbp at Max RPKM",
       y = "Nonsynonymous Fraction at Max RPKM",
       colour = "Max RPKM") +
  theme_classic()

ggsave("microdiversity.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)

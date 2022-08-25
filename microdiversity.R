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
  transmute(Tax = paste(Phylum, Genus, sep = ";"))

microdiversity_df <- data.frame("DiSiperMbpAtAbundMax" = disi_array,
                                "NonsynonimousFractionAtAbundMax" = nsf_array,
                                "AbundMax" = filtered_df$AbundMax,
                                "Tax" = tax_df$Tax) %>%
  filter(NonsynonimousFractionAtAbundMax != 0) %>%
  mutate(Tax = ifelse(NonsynonimousFractionAtAbundMax > 0.5 | DiSiperMbpAtAbundMax > 100000 | AbundMax > 50, Tax, ""))

p <- ggplot(microdiversity_df, aes(x = DiSiperMbpAtAbundMax,
                                   y = NonsynonimousFractionAtAbundMax,
                                   colour = AbundMax)) +
  scale_x_continuous(limits = c(0, 150000),
                     breaks = seq(0, 150000, 25000)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  geom_text_repel(aes(label = Tax),
                  size = 1.8,
                  colour = "black",
                  max.overlaps = 1000) +
  geom_point(size = 1) +
  labs(x = "Divergent Sites per Mbp at Max RPKM",
       y = "Nonsynonymous Fraction at Max RPKM",
       colour = "Max RPKM") +
  scale_colour_gradientn(colours = c("blue", "yellow", "red"),
                         limits = c(0, 65),
                         breaks = seq(0, 60, by = 10)) +
  theme_classic()

ggsave("microdiversity.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)

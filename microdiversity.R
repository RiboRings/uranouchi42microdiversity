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

microdiversity_df <- data.frame("genome" = filtered_df$genome, 
                                "DiSiperMbpAtAbundMax" = disi_array,
                                "NonsynonimousFractionAtAbundMax" = nsf_array,
                                "AbundMax" = filtered_df$AbundMax,
                                "Recurrence" = filtered_df$Recurrence,
                                "Tax" = tax_df$Tax) %>%
  filter(NonsynonimousFractionAtAbundMax != 0) %>%
  mutate(Tax = ifelse(NonsynonimousFractionAtAbundMax > 0.5 | DiSiperMbpAtAbundMax > 100000 | AbundMax > 50, Tax, ""))

p <- ggplot(microdiversity_df, aes(x = DiSiperMbpAtAbundMax / 1000,
                                   y = NonsynonimousFractionAtAbundMax,
                                   colour = Recurrence,
                                   size = AbundMax)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 150),
                     breaks = seq(0, 150, 25)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  scale_size_continuous(limits = c(0, 65),
                        breaks = c(1, seq(20, 60, by = 20))) +
  scale_colour_gradientn(colours = c("blue", "yellow", "red"),
                         limits = c(1, 42),
                         breaks = c(1, seq(10, 40, by = 10))) +
  geom_text_repel(aes(label = Tax),
                  size = 1.8,
                  colour = "black",
                  max.overlaps = 1000) +
  labs(x = "Divergent Sites per Kbp at Max RPKM",
       y = "Nonsynonymous Fraction at Max RPKM",
       colour = "Recurrence",
       size = "Max RPKM") +
  theme_classic()

ggsave("microdiversity.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)

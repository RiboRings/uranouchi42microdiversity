selected_genomes <- top_df$genome

time_series_df1 <- top_df %>%
  select(genome,
         Tax,
         starts_with("nucdiv_")) %>%
  pivot_longer(starts_with("nucdiv"),
               values_to = "NucDiv",
               names_to = "Sample") %>%
  mutate(Sample = gsub("nucdiv_", "", Sample))

time_series_df2 <- top_df %>%
  select(genome,
         Tax,
         ends_with("RPKM")) %>%
  pivot_longer(ends_with("RPKM"),
               values_to = "RPKM",
               names_to = "Sample") %>%
  mutate(Sample = gsub("_RPKM", "", Sample))

time_series_df <- time_series_df1 %>%
  full_join(time_series_df2) %>%
  mutate(Sample = gsub("UU\\w", "", Sample)) %>%
  mutate(Sample = paste(substr(Sample, 1, 2), substr(Sample, 3, 4), substr(Sample, 5, 6), sep = "-")) %>%
  mutate(Sample = as.Date(Sample))

nucdiv_range_df <- time_series_df %>%
  group_by(genome) %>%
  summarise(NucDivRange = max(NucDiv, na.rm = TRUE) - min(NucDiv, na.rm = TRUE))

time_series_df <- time_series_df %>%
  left_join(nucdiv_range_df) %>%
  arrange(desc(NucDivRange))

p <- ggplot(time_series_df, aes(x = Sample,
                                group = Tax)) +
  geom_line(aes(y = RPKM), colour = "Black") +
  geom_line(aes(y = NucDiv * 1000), colour = "Dark Gray") +
  scale_y_continuous(name = "RPKM",
                     sec.axis = sec_axis(~./1000,
                                         name = "Nucleotide Diversity")) +
  facet_wrap(~ NucDivRange * Tax,
             nrow = 4) +
  theme_bw() +
  theme(axis.title.y = element_text(colour = "Black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.right = element_text(colour = "Dark Gray"),
        panel.grid = element_blank())

ggsave("multifacets.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 40,
       height = 10)

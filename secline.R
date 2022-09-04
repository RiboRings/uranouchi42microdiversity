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

time_series_df <- time_series_df %>%
  left_join(subgroup_df) %>%
  group_by(Cluster, Sample) %>%
  summarise(NucDiv = mean(NucDiv),
            RPKM = mean(RPKM))

p <- ggplot(time_series_df, aes(x = Sample,
                                group = Cluster)) +
  geom_line(aes(y = RPKM), colour = "Black") +
  geom_line(aes(y = NucDiv * 1000), colour = "Dark Gray") +
  scale_x_date(date_breaks = "month",
               date_labels = "%b %y") +
  scale_y_continuous(name = "RPKM",
                     sec.axis = sec_axis(~./1000,
                                         name = "Nucleotide Diversity")) +
  facet_wrap(~ Cluster,
             nrow = length(unique(time_series_df$Cluster))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(colour = "Black"),
        axis.title.y.right = element_text(colour = "Dark Gray"),
        panel.grid = element_blank())

ggsave("secline.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 7,
       height = 30)

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
  summarise(NucDivStDev = sd(NucDiv),
            RPKMStDev = sd(RPKM),
            NucDiv = mean(NucDiv),
            RPKM = mean(RPKM))

p_clust <- ggplot(time_series_df, aes(x = Sample,
                                      group = Cluster)) +
  geom_line(aes(y = RPKM), colour = "#F8766D") +
  geom_errorbar(aes(ymax = RPKM + RPKMStDev,
                    ymin = RPKM - RPKMStDev),
                alpha = 0.2) +
  geom_line(aes(y = NucDiv * 1000), colour = "#00BFC4") +
  geom_errorbar(aes(ymax = (NucDiv + NucDivStDev) * 1000,
                    ymin = (NucDiv - NucDivStDev) * 1000),
                alpha = 0.2) +
  scale_x_date(date_breaks = "month",
               date_labels = "%b %y") +
  scale_y_continuous(name = "RPKM",
                     limits = c(0, 65),
                     breaks = seq(0, 60, by = 10),
                     sec.axis = sec_axis(~./1000,
                                         name = "Nucleotide Diversity")) +
  facet_wrap(~ Cluster,
             ncol = 3,
             nrow = 2) +
  labs(x = "Time") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(colour = "#F8766D"),
        axis.title.y.right = element_text(colour = "#00BFC4"),
        panel.grid = element_blank())

p <- p_clust / p_max +
  plot_annotation(tag_levels = "A")

ggsave("secline.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 15,
       height = 15)

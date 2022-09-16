lolo_df <- heatmap_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "genome") %>%
  separate(col = genome, into = c("Phylum", "Order", "ID"), sep = ";") %>%
  pivot_longer(cols = c(-Phylum, -Order, -ID),
               values_to = "RPKM",
               names_to = "Sample") %>%
  group_by(Phylum, Order, Sample) %>%
  summarise(RPKM = sum(RPKM))

lolo_df <- lolo_df %>%
  ungroup() %>%
  mutate(Sample = gsub("\\.1", " *", Sample),
         Phylum = gsub("p__", "", Phylum),
         Order = gsub("o__", "", Order)) %>%
  mutate(Time = format(as.Date(Sample), "%d %b %y")) %>%
  mutate(Time = ifelse(str_detect(Sample, "\\*"), paste(Time, "*"), Time)) %>%
  group_by(Sample) %>%
  mutate(Percent = RPKM / sum(RPKM) * 100) %>%
  arrange(desc(Percent))

toto_df <- lolo_df %>%
  ungroup() %>%
  group_by(Phylum) %>%
  summarise(Percent = sum(Percent)) %>%
  arrange(desc(Percent))

lolo_df$Phylum[lolo_df$Phylum %in% toto_df$Phylum[11:nrow(toto_df)]] <- "Other"
toto_df$Phylum[11:nrow(toto_df)] <- "Other"
lolo_df$Phylum <- factor(lolo_df$Phylum, levels = unique(toto_df$Phylum))

p40 <- ggplot(lolo_df, aes(x = Time,
                           y = Percent,
                           fill = Phylum)) +
  geom_col() +
  scale_y_continuous(limits = c(-1, 101),
                     breaks = seq(0, 100, by = 10)) +
  labs(x = "Time (#d #m #y)",
       y = "Percent RPKM (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

momo_df <- lolo_df %>%
  ungroup() %>%
  group_by(Order) %>%
  summarise(Percent = sum(Percent)) %>%
  arrange(desc(Percent))

momo_df$Order[is.na(momo_df$Order)] <- "Other"

lolo_df$Order[lolo_df$Order %in% momo_df$Order[11:nrow(momo_df)]] <- "Other"
momo_df$Order[11:nrow(momo_df)] <- "Other"
lolo_df$Order <- factor(lolo_df$Order, levels = unique(momo_df$Order))

p80 <- ggplot(lolo_df, aes(x = Time,
                           y = Percent,
                           fill = Order)) +
  geom_col() +
  scale_y_continuous(limits = c(-1, 101),
                     breaks = seq(0, 100, by = 10)) +
  labs(x = "Time (#d #m #y)",
       y = "Percent RPKM (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

layout <- "
AA
BB
CC
CC
CC
DD
DD
"

p <- (p40 / p80 / ggtree2 / ggtree1) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")

ggsave("compobar.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 20,
       height = 30)
 
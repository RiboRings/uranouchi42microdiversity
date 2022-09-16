ubiquity_mat <- micro_df %>%
  filter(genome %in% filtered_df$genome) %>%
  select(genome, breadth, sample) %>%
  pivot_wider(names_from = sample,
              values_from = breadth,
              values_fill = 0) %>%
  column_to_rownames(var = "genome") %>%
  as.matrix()

ubiquity <- sapply(iter(ubiquity_mat, by = "row"), function(x) sum(x > 0.50))
ubiquity_df <- data.frame("genome" = rownames(ubiquity_mat), "Recurrence" = ubiquity)
filtered_df <- filtered_df %>%
  full_join(ubiquity_df)

p1 <- ggplot(ubiquity_df, aes(x = Recurrence)) +
  geom_bar(colour = "black",
           fill = "#F8766D") +
  scale_x_continuous(limits = c(0, 43),
                     breaks = c(1, seq(5, 40, by = 5))) +
  scale_y_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, by = 10)) +
  labs(x = "MAG Recurrence",
       y = "Number of MAGs") +
  theme_classic() +
  theme(panel.grid = element_blank())

p7 <- ggplot(ubiquity_df, aes(x = Recurrence)) +
  stat_bin(aes(y = cumsum(..count..)),
           bins = 44,
           colour = "black",
           fill = "#00BFC4") +
  geom_bar(colour = "black",
           fill = "#F8766D") +
  scale_x_continuous(limits = c(0, 43),
                     breaks = c(1, seq(5, 40, by = 5))) +
  scale_y_continuous(limits = c(0, 1200),
                     breaks = seq(0, 1500, by = 100)) +
  labs(x = "MAG Recurrence",
       y = "Cumulative Number of MAGs") +
  theme_classic() +
  theme(panel.grid = element_blank())

ubiquity2 <- lapply(1:nrow(ubiquity_mat), function(x) names(ubiquity_mat[x, ubiquity_mat[x, ] > 0.50]))
tmp_mat <- matrix(nrow = 1132, ncol = 4)
for (i in 1:length(ubiquity2)) {
tmp_mat[i, 1] <- sum(substr(ubiquity2[[i]], 6, 7) %in% c("12", "01", "02"))
tmp_mat[i, 2] <- sum(substr(ubiquity2[[i]], 6, 7) %in% c("06", "07", "08"))
tmp_mat[i, 3] <- sum(substr(ubiquity2[[i]], 6, 7) %in% c("03", "04", "05"))
tmp_mat[i, 4] <- sum(substr(ubiquity2[[i]], 6, 7) %in% c("09", "10", "11"))
}
tmp_df <- as.data.frame(tmp_mat)
colnames(tmp_df) <- c("Winter", "Summer", "Spring", "Autumn")
tmp_df <- tmp_df %>%
  pivot_longer(names_to = "Season",
               values_to = "Count",
               cols = c("Winter", "Summer", "Autumn", "Spring"))

tmp_df$level[tmp_df$Season == "Winter"] <- 1
tmp_df$level[tmp_df$Season == "Spring"] <- 2
tmp_df$level[tmp_df$Season == "Summer"] <- 3
tmp_df$level[tmp_df$Season == "Autumn"] <- 4

p2 <- ggplot(tmp_df, aes(x = reorder(Season, level), y = Count)) +
  geom_boxplot() +
  labs(x = "Season",
       y = "MAG Recurrence") +
  theme_classic()

samples <- df %>%
  select(ends_with("_RPKM")) %>%
  colnames()
samples <- gsub("_RPKM", "", samples)
samples <- gsub("UU\\w", "", samples)
samples <- paste(substr(samples, 1, 2), substr(samples, 3, 4), substr(samples, 5, 6), sep = "-")
samples <- as.Date(samples)
samples_df <- data.frame(Time = samples)

ubiquity_df$ExtractionDate <- gsub("UU\\w(\\d*)_\\d*", "\\1", ubiquity_df$genome)
ubiquity_df$ExtractionDate <- paste(substr(ubiquity_df$ExtractionDate, 1, 2), substr(ubiquity_df$ExtractionDate, 3, 4), substr(ubiquity_df$ExtractionDate, 5, 6), sep = "-")
ubiquity_df$ExtractionDate <- as.Date(ubiquity_df$ExtractionDate)

p3 <- ggplot() +
  geom_bar(data = samples_df,
           stat = "bin",
           aes(x = Time,
               y = ..count..),
           bins = 23,
           colour = "black",
           fill = "white") +
  geom_density(data = ubiquity_df,
               aes(x = ExtractionDate,
                   y = ..count.. * 2),
               colour = "#F8766D") +
  scale_x_date(date_breaks = "month",
               date_labels = "%b %y") +
  scale_y_continuous(limits = c(0, 8),
                     breaks = seq(0, 8),
                     sec.axis = sec_axis(~./2,
                                         name = "MAG Density",
                                         breaks = seq(0, 4, by = 0.5))) +
  labs(x = "Time",
       y = "Number of Samples") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y.right = element_text(colour = "#F8766D"))

p <- ((p1 / p7) | (p2 / p3)) +
  plot_annotation(tag_levels = "A")

ggsave("recurrence.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 14,
       height = 7)

n_top_taxa <- 1132

top_df <- filtered_df %>%
  slice_max(AbundMax, n = n_top_taxa) %>%
  mutate(Tax = paste(Order, genome, sep = ";"))

col_fun1 <- colorRamp2(c(0, 0.1, 1, 5, 10, 20, 40, 60), col = rev(rainbow(8)))
col_fun2 <- colorRamp2(c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075), col = rev(rainbow(8)))

heatmap_mat <- top_df %>%
  select(Tax, ends_with("RPKM")) %>%
  column_to_rownames(var = "Tax") %>%
  as.matrix()

tmp <- gsub("UU\\w", "", gsub("_RPKM", "", colnames(heatmap_mat)))
colnames(heatmap_mat) <- paste(substr(tmp, 1, 2), substr(tmp, 3, 4), substr(tmp, 5, 6), sep = "-")
heatmap_mat <- heatmap_mat[ , order(colnames(heatmap_mat))]
colnames(heatmap_mat) <- format(as.Date(colnames(heatmap_mat)), "%b %y")

h1 <- Heatmap(heatmap_mat,
             heatmap_legend_param = list(at = c(0, 0.1, 1, 5, 10, 20, 40, 60)),
             col = col_fun1,
             name = "RPKM",
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             row_title = "MAG",
             column_title = "Sample",
             row_names_gp = gpar(fontsize = 0.7),
             column_names_gp = gpar(fontsize = 8),
             column_names_rot = 45,
             width = ncol(heatmap_mat)*unit(4.5, "mm"), 
             height = nrow(heatmap_mat)*unit(0.5, "mm"))

nd_mat <- top_df %>%
  select(Tax, starts_with("nucdiv_")) %>%
  column_to_rownames(var = "Tax") %>%
  as.matrix()
colnames(nd_mat) <- paste(substr(tmp, 1, 2), substr(tmp, 3, 4), substr(tmp, 5, 6), sep = "-")
nd_mat <- nd_mat[ , order(colnames(nd_mat))]
colnames(nd_mat) <- format(as.Date(colnames(nd_mat)), "%b %y")

ha <- HeatmapAnnotation(`Max RPKM` = filtered_df$AbundMax,
                        `Max ND` = rowMaxs(top_df %>%
                                             select(starts_with("nucdiv_")) %>%
                                             as.matrix()),
                        which = "row",
                        col = list(`Max RPKM` = col_fun1,
                                   `Max ND` = col_fun2),
                        show_legend = FALSE)

h2 <- Heatmap(nd_mat,
             heatmap_legend_param = list(at = c(0, 0.02, 0.04, 0.06, 0.075)),
             col = col_fun2,
             name = "ND",
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             row_title = "MAG",
             column_title = "Sample",
             row_names_gp = gpar(fontsize = 0.7),
             column_names_gp = gpar(fontsize = 8),
             column_names_rot = 45,
             right_annotation = ha,
             width = ncol(mat)*unit(4.5, "mm"), 
             height = nrow(mat)*unit(0.5, "mm"))

ht_list <- h1 + h2

pdf("results/composition.pdf",
    width = 20,
    height = 30)

draw(ht_list)

dev.off()

max_df <- data.frame(AbundMax = filtered_df$AbundMax,
                     NucDivMax = rowMaxs(top_df %>%
                                           select(starts_with("nucdiv_")) %>%
                                           as.matrix()))

p <- ggplot(max_df, aes(x = AbundMax, NucDivMax)) +
  geom_point() +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 65),
                     breaks = seq(0, 60, by = 10)) +
  labs(x = "Max RPKM",
       y = "Max Nucleotide Diversity") +
  theme_classic()

ggsave("nd_rpkm.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)
                     
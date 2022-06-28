n_top_taxa <- 50

top_df <- filtered_df %>%
  slice_max(AbundMax, n = n_top_taxa) %>%
  mutate(Tax = paste(Order, genome, sep = ";"))

heatmap_mat <- top_df %>%
  select(Tax, ends_with("RPKM")) %>%
  column_to_rownames(var = "Tax") %>%
  as.matrix()

tmp <- gsub("UU\\w", "", gsub("_RPKM", "", colnames(heatmap_mat)))
colnames(heatmap_mat) <- paste(substr(tmp, 1, 2), substr(tmp, 3, 4), substr(tmp, 5, 6), sep = "-")
heatmap_mat <- heatmap_mat[ , order(colnames(heatmap_mat))]

ha <- HeatmapAnnotation(ND = rowMeans(top_df %>% select(starts_with("nucdiv"))),
                         which = "row")

pdf("results/composition.pdf",
    width = 20,
    height = 10)

draw(Heatmap(heatmap_mat,
             heatmap_legend_param = list(at = seq(floor(min(heatmap_mat)), ceiling(max(heatmap_mat)), by = 10)),
             col = c("blue", "cyan", "green", "yellow", "red"),
             name = "RPKM",
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             row_title = "MAG",
             column_title = "Sample",
             row_names_gp = gpar(fontsize = 8),
             left_annotation = ha))

dev.off()

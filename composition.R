n_top_taxa <- 1132

top_df <- filtered_df %>%
  slice_max(AbundMax, n = n_top_taxa) %>%
  mutate(Tax = paste(Phylum, Order, gsub("MAG ", "", ID), sep = ";"))

col_fun1 <- colorRamp2(c(0, 0.1, 1, 5, 10, 20, 40, 60), col = rev(rainbow(8)))
col_fun2 <- colorRamp2(c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075), col = rev(rainbow(8)))

heatmap_mat <- top_df %>%
  select(Tax, ends_with("RPKM")) %>%
  column_to_rownames(var = "Tax") %>%
  as.matrix()

heatmap_mat <- heatmap_mat[order(rownames(heatmap_mat)), ]

tmp <- gsub("UU\\w", "", gsub("_RPKM", "", colnames(heatmap_mat)))
colnames(heatmap_mat) <- paste(substr(tmp, 1, 2), substr(tmp, 3, 4), substr(tmp, 5, 6), sep = "-")
heatmap_mat <- heatmap_mat[ , order(colnames(heatmap_mat))]
colnames(heatmap_mat) <- format(as.Date(colnames(heatmap_mat)), "%b %y")

h1 <- Heatmap(heatmap_mat,
              heatmap_legend_param = list(at = c(0, 0.1, 1, 5, 10, 20, 40, 60),
                                          break_dist = 1,
                                          legend_height = unit(4, "cm")),
              col = col_fun1,
              name = "RPKM",
              cluster_rows = TRUE,
              row_km = 5,
              row_km_repeats = 20,
              cluster_columns = FALSE,
              show_row_dend = FALSE,
              row_title = "MAG",
              row_names_side = "left",
              column_title = "Time Series of Abundance (RPKM)",
              row_names_gp = gpar(fontsize = 0.7),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 45,
              width = ncol(heatmap_mat) * unit(4.8, "mm"), 
              height = nrow(heatmap_mat) * unit(0.6, "mm"))

nd_mat <- top_df %>%
  select(Tax, starts_with("nucdiv_")) %>%
  column_to_rownames(var = "Tax") %>%
  as.matrix()

nd_mat <- nd_mat[order(rownames(nd_mat)), ]

colnames(nd_mat) <- paste(substr(tmp, 1, 2), substr(tmp, 3, 4), substr(tmp, 5, 6), sep = "-")
nd_mat <- nd_mat[ , order(colnames(nd_mat))]
colnames(nd_mat) <- format(as.Date(colnames(nd_mat)), "%b %y")

h2 <- Heatmap(nd_mat,
              heatmap_legend_param = list(at = c(0, 0.02, 0.04, 0.06, 0.075),
                                          break_dist = 1,
                                          legend_height = unit(4, "cm")),
              col = col_fun2,
              name = "ND",
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              row_title = "MAG",
              row_names_side = "left",
              column_title = "Time Series of Nucleotide Diversity (ND)",
              row_names_gp = gpar(fontsize = 0.7),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 45,
              width = ncol(nd_mat) * unit(4.8, "mm"), 
              height = nrow(nd_mat) * unit(0.6, "mm"))

ht_list <- h1 + h2

ht_list <- draw(ht_list)
r.dend <- row_dend(ht_list)
rcl.list <- row_order(ht_list)

# loop to extract genes for each cluster.
for (i in 1:length(row_order(ht_list))){
  if (i == 1) {
    clu <- t(t(row.names(heatmap_mat[row_order(ht_list)[[i]], ])))
    out <- cbind(clu, paste("Cluster", i, sep = " "))
    colnames(out) <- c("Tax", "Cluster")
  } else {
    clu <- t(t(row.names(heatmap_mat[row_order(ht_list)[[i]], ])))
    clu <- cbind(clu, paste("Cluster", i, sep = " "))
    out <- rbind(out, clu)
  }
}

subgroup_df <- as.data.frame(out)

subgroup_df <- top_df %>%
  left_join(subgroup_df) %>%
  select(Tax, Cluster)

pdf("results/composition.pdf",
    width = 20,
    height = 30)

draw(ht_list)

dev.off()
max_df <- data.frame(AbundMax = filtered_df$AbundMax,
                     Recurrence = filtered_df$Recurrence,
                     Tax = subgroup_df$Tax,
                     Cluster = subgroup_df$Cluster,
                     NucDivMax = rowMaxs(top_df %>%
                                           select(starts_with("nucdiv_")) %>%
                                           as.matrix()))

max_df <- max_df %>%
  mutate(Tax = ifelse(AbundMax > 25, Tax, ""))

p_max <- ggplot(max_df, aes(x = AbundMax,
                            y = NucDivMax,
                            colour = Recurrence,
                            shape = Cluster)) +
  geom_point() +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 65),
                     breaks = seq(0, 60, by = 10)) +
  scale_colour_gradientn(colours = c("blue", "yellow", "red"),
                         limits = c(1, 42),
                         breaks = c(1, seq(10, 40, by = 10))) +
  geom_text_repel(aes(label = Tax),
                  size = 1.8,
                  colour = "black",
                  max.overlaps = 1000) +
  labs(x = "Max RPKM",
       y = "Max Nucleotide Diversity") +
  theme_classic()

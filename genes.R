func_annot <- read_tsv("data/special_func_annot.gff",
                       col_names = c("pass", "gene", "KO", "thrshld", "score", "Evalue", "Definition"))
func_cats <- read_tsv("data/KEGG_pathway_ko_uniq.txt")

func_annot_cats <- func_annot %>%
  left_join(func_cats, by = c("KO" = "ko")) %>%
  mutate(gene = gsub("___", "", gene))

abc <- c("UUJ170105", "UUJ170306", "UUJ170412", "UUJ170501", "UUJ170523", "UUJ170612",
         "UUJ170615", "UUF170620", "UUJ170620", "UUF170621", "UUJ170623", "UUJ170626",
         "UUJ170630", "UUJ170706", "UUJ170714", "UUJ170721", "UUJ170726", "UUJ170804",
         "UUJ170818", "UUJ170927", "UUJ171027", "UUJ171110", "UUM171110", "UUJ171113",
         "UUJ171218", "UUJ180110", "UUJ180219", "UUJ180313", "UUJ180410", "UUJ180510",
         "UUJ180528", "UUJ180622", "UUJ180629", "UUJ180718", "UUJ180725", "UUJ180802",
         "UUJ180809", "UUJ180815", "UUJ180827", "UUM180827", "UUJ180903", "UUJ180925")
trick_df <- data.frame(AbundMaxIdx = 1:42, AbundMaxSample = abc)
gen_list <- data.frame(genome = filtered_df$genome, AbundMaxIdx = filtered_df$AbundMaxIdx) %>%
  left_join(trick_df)

mapping <- read_tsv("data/mapping.stb", col_names = c("scaffold", "genome")) %>%
  mutate(genome = gsub("matabat2bin.", "", gsub(".fa", "", genome))) %>%
  right_join(gen_list)

gene_files <- list.files(path = "data/gene_info", pattern = "*_gene_info.tsv")

gene_df_list <- lapply(paste("data", "gene_info", gene_files, sep = "/"),
                       gene_loader,
                       method = read_tsv)

gene_df <- merge_all(gene_df_list)
func_annot_df <- gene_df %>%
  left_join(func_annot_cats) %>%
  mutate(NSF = SNV_N_count / SNV_count)

# func_annot_df <- func_annot_df %>%
#   filter(coverage > 10)
#test <- func_annot_df %>%
#  filter(level2_pathway_name %in% colnames(func_mat))
#  p

func_mat <- func_annot_df %>%
  group_by(genome, level2_pathway_name) %>%
  summarise(NSF = mean(NSF, na.rm = TRUE)) %>%
  pivot_wider(names_from = level2_pathway_name,
              values_from = NSF,
              values_fill = NA) %>%
  column_to_rownames(var = "genome") %>%
  as.matrix()

interesting_guys <- microdiversity_df %>%
  filter(genome %in% rownames(func_mat)) %>%
  arrange(desc(NonsynonimousFractionAtAbundMax), DiSiperMbpAtAbundMax)

func_mat <- read_csv("data/raw_gene_mat.csv") %>%
  inner_join(interesting_guys) %>%
  arrange(desc(NonsynonimousFractionAtAbundMax), DiSiperMbpAtAbundMax) %>%
  select(-NonsynonimousFractionAtAbundMax,
         -DiSiperMbpAtAbundMax,
         -AbundMax,
         -Tax) %>%
  column_to_rownames(var = "genome") %>%
  as.matrix()

func_mat <- func_mat[ , !(colnames(func_mat) %in% c("NA",
                                                    "Not included in regular maps",
                                                    "Poorly characterized",
                                                    "Aging",
                                                    "Immune system",
                                                    "Substance dependence",
                                                    "Cardiovascular disease",
                                                    "Neurodegenerative disease",
                                                    "Endocrine and metabolic disease",
                                                    "Nervous system",
                                                    "Endocrine system",
                                                    "Immune disease",
                                                    "Cancer: overview",
                                                    "Cancer: specific types",
                                                    "Infectious disease: parasitic",
                                                    "Infectious disease: viral",
                                                    "Cellular community - eukaryotes",
                                                    "Viral protein families",
                                                    "Information processing in viruses",
                                                    "Circulatory system",
                                                    "Excretory system",
                                                    "Digestive system",
                                                    "Development and regeneration"))]

#raw_gene_mat <- as.data.frame(func_mat) %>%
#  rownames_to_column(var = "genome")
#write_csv(as.data.frame(raw_gene_mat), "data/raw_gene_mat.csv")
#
#func_mat <- func_mat[rownames(func_mat) %in% unique(interesting_guys$genome), ]

clean_func_mat <- remove_problematic_combs(func_mat)

interesting_guys <- interesting_guys %>%
  filter(genome %in% rownames(clean_func_mat))

ha_row <- rowAnnotation(`Divergent Sites per Kbp` = anno_barplot(interesting_guys$DiSiperMbpAtAbundMax / 1000),
                        `Max RPKM` = anno_barplot(interesting_guys$AbundMax),
                        annotation_name_rot = 90,
                        annotation_name_gp = gpar(fontsize = 8))

clean_func_mat <- clean_func_mat[ , order(colMedians(clean_func_mat, na.rm = TRUE))]
ha_col <- columnAnnotation(Distribution = anno_boxplot(clean_func_mat),
                           annotation_name_gp = gpar(fontsize = 8))

fun_color_range <- colorRampPalette(c("#1b98e0", "red"))
my_colors <- fun_color_range(101)
col_fun <- colorRamp2(breaks = seq(0, 1, by = 0.01), colors = my_colors)

pdf("results/func_hm.pdf",
    width = 10,
    height = 15)

draw(Heatmap(clean_func_mat,
             name = "Mean NSF",
             heatmap_legend_param = list(at = seq(0, 1, by = 0.1)),
             col = col_fun,
             na_col = "grey",
             row_title = "MAG",
             column_title = "Functional Group",
             show_row_names = FALSE,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             column_names_gp = gpar(fontsize = 8),
             column_names_rot = 45,
             right_annotation = ha_row,
             top_annotation = ha_col,
             width = ncol(mat) * unit(4, "mm"), 
             height = nrow(mat) * unit(0.18, "mm")))

dev.off()

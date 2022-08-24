func_annot <- read_tsv("data/special_func_annot.gff",
                       col_names = c("pass", "gene", "KO", "thrshld", "score", "Evalue", "Definition"))
func_cats <- read_tsv("../KEGG_pathway_ko_uniq.txt")

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

func_mat <- func_annot_df %>%
  group_by(genome, level2_pathway_name) %>%
  summarise(NSF = mean(NSF)) %>%
  pivot_wider(names_from = level2_pathway_name,
              values_from = NSF,
              values_fill = NA) %>%
  column_to_rownames(var = "genome") %>%
  as.matrix()

my_tax <- top_df$Tax[top_df$genome %in% rownames(func_mat)]
my_tax_order <- order(gsub(".*\\;", "", gsub("o__", "", my_tax)))
my_mat_order <- order(rownames(func_mat))

my_tax <- my_tax[my_tax_order]
func_mat <- func_mat[my_mat_order, ]

rownames(func_mat) <- my_tax

fun_color_range <- colorRampPalette(c("#1b98e0", "red"))
my_colors <- fun_color_range(101)
col_fun <- colorRamp2(breaks = seq(0, 1, by = 0.01), colors = my_colors)

pdf("results/func_hm.pdf",
    width = 10,
    height = 20)
draw(Heatmap(clean_func_mat,
             name = "Mean NSF",
             heatmap_legend_param = list(at = seq(0, 1, by = 0.1)),
             col = col_fun,
             na_col = "grey",
             row_title = "MAG",
             column_title = "Functional Group",
             row_names_gp = gpar(fontsize = 4),
             column_names_gp = gpar(fontsize = 8),
             column_names_rot = 45))
dev.off()

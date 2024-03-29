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

gene_df <- gene_df %>%
  filter(!is.na(SNV_count)) %>%
  mutate(SNV_count = ifelse(SNV_count == 0, 1, SNV_count))

func_annot_df <- gene_df %>%
  left_join(func_annot_cats) %>%
  mutate(NSF = SNV_N_count / SNV_count)

interesting_guys <- microdiversity_df %>%
  arrange(desc(NonsynonimousFractionAtAbundMax), DiSiperMbpAtAbundMax)

func_bubble <- func_annot_df %>%
  left_join(interesting_guys) %>%
  rename(FunctionalCategory = level2_pathway_name,
         NonsynonymousFraction = NSF) %>%
  filter(!is.na(FunctionalCategory) & !(FunctionalCategory %in% c("NA",
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
                                                                "Development and regeneration",
                                                                "Unclassified: genetic information processing",
                                                                "Unclassified: signaling and cellular processes",
                                                                "Unclassified: metabolism")))

write_csv(func_bubble, "data/raw_gene_mat.csv")
#func_bubble <- read_csv("data/raw_gene_mat.csv")

func_bubble$Group <- ""
func_bubble$Group[func_bubble$DiSiperMbpAtAbundMax > 100000] <- "Number of Divergent Sites per Kbp > 100"
func_bubble$Group[func_bubble$NonsynonimousFractionAtAbundMax > 0.5] <- "Percentage of non-synonymous Mutations at Genome Level > 50%"

shapiro_tests <- by(func_bubble,
                    list(func_bubble$FunctionalCategory),
                    function(x) shapiro.test(x$NonsynonymousFraction[x$Group != ""]))

t_tests <- by(func_bubble,
              list(func_bubble$FunctionalCategory),
              function(x) t.test(x$NonsynonymousFraction[x$Group == "Number of Divergent Sites per Kbp > 100"],
                                 x$NonsynonymousFraction[x$Group == "Percentage of non-synonymous Mutations at Genome Level > 50%"]))

shapiro_test_pvalues <- sapply(shapiro_tests, function(x) x$p.value)
t_test_pvalues <- sapply(t_tests, function(x) x$p.value)

p_value_df <- data.frame(ShapiroTest = shapiro_test_pvalues,
                         tTest <- t_test_pvalues) %>%
  rownames_to_column(var = "FunctionalCategory") %>%
  mutate(Normality = ifelse(ShapiroTest <= 0.05, TRUE, FALSE),
         Significance = ifelse(tTest <= 0.05, TRUE, FALSE))

func_bubble <- func_bubble %>%
  mutate(Significance = ifelse(FunctionalCategory %in% p_value_df$FunctionalCategory[p_value_df$Significance == TRUE], "*", "ns"),
         Normality = ifelse(FunctionalCategory %in% p_value_df$FunctionalCategory[p_value_df$Normality == TRUE], "*", "ns"))

max_nsf <- max(func_bubble$NonsynonimousFractionAtAbundMax)
min_nsf <- min(func_bubble$NonsynonimousFractionAtAbundMax)

pf <- ggplot(func_bubble, aes(x = NonsynonymousFraction * 100,
                              y = reorder(FunctionalCategory, NonsynonymousFraction, median, na.rm = TRUE))) +
  geom_boxplot() +
  geom_vline(aes(xintercept = median(NonsynonymousFraction))) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 10)) +
  labs(x = "Percentage of non-synonymous Mutations at Gene Level (%)",
       y = "Functional Category",
       colour = "Groups") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 45))

ggsave("func_annot_overview.pdf",
       plot = pf,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)

po <- ggplot(subset(func_bubble, Group != ""), aes(x = reorder(FunctionalCategory, NonsynonymousFraction, median, na.rm = TRUE),
                                                   y = NonsynonymousFraction * 100)) +
  geom_boxplot(aes(fill = Group)) +
  stat_compare_means(data = subset(func_bubble, Group != "" & Significance == "*" & Normality == "*"),
                     aes(group = Group),
                     method = "t.test",
                     label = "p.signif",
                     label.y = 100,
                     colour = "blue",
                     size = 5) +
  stat_compare_means(data = subset(func_bubble, Group != "" & Significance == "*" & Normality == "ns"),
                     aes(group = Group),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y = 100,
                     colour = "red",
                     size = 5) +
  labs(x = "Functional Category",
       y = "Percentage of non-synonymous Mutations at Gene Level (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

p <- (pf / po) +
  plot_annotation(tag_levels = "A")

ggsave("microfunc.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 15,
       height = 15)

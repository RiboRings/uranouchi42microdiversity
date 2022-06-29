snvs_files <- list.files(pattern = 'UU.*_snvs_summary.csv_short.csv',
                         path = "data")

snvs_file_list <- mclapply(paste0("data/", snvs_files),
                           sample_loader,
                           mc.cores = detectCores() - 1)

snvs_df <- snvs_file_list %>%
  merge_all() %>%
  filter(genome %in% unique(top_df$genome)) %>%
  transmute(genome,
            scaffold,
            position = as.factor(position),
            con_base,
            A, C, T, G,
            sample = gsub("_snvs_summary.csv_short.csv", "", sample))

for (cur_genome in selected_genomes) {
  
  consensus_base_df <- snvs_df %>%
    filter(genome == cur_genome) %>%
    select(scaffold,
           position,
           con_base,
           sample) %>%
    pivot_wider(names_from = sample,
                values_from = con_base,
                values_fill = NA)
  
  names(consensus_base_df) <- gsub("(UU\\w)(\\d{6})", "\\2\\1", names(consenus_base_df))
  consensus_base_df <- consensus_base_df[ , order(names(consensus_base_df))]

  consensus_base_mat <- consensus_base_df %>%
    mutate(name = paste(scaffold, position, sep = ";")) %>%
    column_to_rownames(var = "name") %>%
    select(ends_with("UUM"),
           ends_with("UUJ"),
           ends_with("UUF")) %>%
    as.matrix()
  
  count_list <- list()
  
  for (cur_sample in 1:42) {
    
    sample_mat <- consensus_base_mat[!is.na(consensus_base_mat[ , cur_sample]), ]
    
    count_list[[cur_sample]] <- sapply(iter(sample_mat, by = "col"),
                                       function(x) sum(sample_mat[ , cur_sample] == x, na.rm = TRUE))
    
  }
  
  count_df <- count_list %>%
    melt() %>%
    transmute(RefSample = as.factor(L1),
              Count = value)
  
  count_df$Time <- as.factor(rep(seq(1, 42), 42))
  
  up_lim <- ceiling(max(count_df$Count) / 10000) * 10000
  down_lim <- floor(min(count_df$Count) / 10000) * 10000
  
  p1 <- ggplot(count_df, aes(x = Time,
                             y = Count,
                             colour = RefSample,
                             group = RefSample)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits = c(down_lim, up_lim),
                       breaks = seq(down_lim, up_lim, by = (up_lim - down_lim) / 5)) +
    labs(y = "Shared Polymorphic Sites",
         colour = "Reference Sample") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  bar_array <- c()
  for (i in 1:nrow(consensus_base_mat)) {
    
    opts <- unique(consensus_base_mat[i, ])
    bar_array <- append(bar_array,
                        max(sapply(opts, function(x) sum(x == consensus_base_mat[i, ], na.rm = TRUE))))
    
  }
  
  bar_df <- as.data.frame(bar_array)
  
  p2 <- ggplot(bar_df, aes(x = bar_array)) +
    geom_bar(aes(y = (..count..) / sum(..count..))) +
    scale_x_reverse(limits = c(43, 0),
                    breaks = seq(42, 1, by = -1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Consensus Base Recurrence",
         y = "Number of SNVs") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  genome_time_series_df <- time_series_df %>%
    filter(genome == cur_genome)
  
  p3 <- ggplot(genome_time_series_df, aes(x = Sample)) +
    geom_line(aes(y = RPKM), colour = "Black") +
    geom_line(aes(y = NucDiv * 1000), colour = "Dark Gray") +
    scale_x_continuous(limits = c(1, 42),
                       breaks = seq(1, 42)) +
    scale_y_continuous(name = "RPKM",
                       limits = c(0, 50),
                       breaks = seq(0, 50, by = 10),
                       sec.axis = sec_axis(~./1000,
                                           name = "Nucleotide Diversity")) +
    theme_bw() +
    theme(axis.title.y = element_text(colour = "Black"),
          axis.title.y.right = element_text(colour = "Dark Gray"),
          panel.grid = element_blank())
  
  tax_df <- top_df %>%
    filter(genome == cur_genome) %>%
    transmute(Tax = paste(Phylum, Order, Family, sep = ";"))
  
  p <- (p2 / p3 | p1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Dynamics of within-species Variation",
                    subtitle = paste0("Genome: ", cur_genome, ". Taxonomy: ", tax_df$Tax))
  
  ggsave(paste0("polymorph", cur_genome, ".pdf"),
         plot = p,
         height = 10,
         width = 20,
         path = "results",
         device = "pdf")
  
}

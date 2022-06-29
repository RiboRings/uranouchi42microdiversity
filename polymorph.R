snvs_df <- read_csv("data/snvs.csv")

snvs_df <- snvs_df %>%
  transmute(genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)),
            scaffold,
            position,
            con_base,
            A, C, T, G,
            sample)

for (cur_genome in selected_genomes) {
  
  message(paste("analysing genome", cur_genome))
  
  consensus_base_df <- snvs_df %>%
    filter(genome == cur_genome) %>%
    select(scaffold,
           position,
           con_base,
           sample) %>%
    pivot_wider(names_from = sample,
                values_from = con_base,
                values_fill = NA)
  
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
    
    message(paste("computing sample", cur_sample))
    
    sample_mat <- consensus_base_mat[!is.na(consensus_base_mat[ , cur_sample]), ]
    
    count_list[[cur_sample]] <- sapply(iter(sample_mat, by = "col"),
                                       function(x) sum(sample_mat[ , cur_sample] == x, na.rm = TRUE))
    
  }
  
  count_df <- count_list %>%
    melt() %>%
    transmute(RefSample = as.factor(L1),
              Count = value)
  
  count_df$Time <- as.factor(rep(seq(1, 42), 42))
  
  count_df <- count_df %>%
    group_by(RefSample) %>%
    summarise(Count,
              Time,
              Percent = Count / max(Count))
  
  message("generating plot 1")
  
  p1 <- ggplot(count_df, aes(x = Time,
                             y = Percent,
                             colour = RefSample,
                             group = RefSample)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.1)) +
    labs(y = "Percent of Shared Polymorphic Sites",
         colour = "Reference Sample") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  message("plot 1 done")
  
  bar_array <- c()
  for (i in 1:nrow(consensus_base_mat)) {
    
    opts <- unique(consensus_base_mat[i, ])
    bar_array <- append(bar_array,
                        max(sapply(opts, function(x) sum(x == consensus_base_mat[i, ], na.rm = TRUE))))
    
  }
  
  bar_df <- as.data.frame(bar_array)
  
  message("generating plot 2")
  
  p2 <- ggplot(bar_df, aes(x = bar_array)) +
    geom_bar(aes(y = (..count..) / sum(..count..))) +
    scale_x_reverse(limits = c(43, 0),
                    breaks = seq(42, 1, by = -1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Consensus Base Recurrence",
         y = "Number of SNVs") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  message("plot 2 done")
  
  genome_time_series_df <- time_series_df %>%
    filter(genome == cur_genome)
  
    message("generating plot 3")
  
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
  
  message("plot 3 done")
  
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
  
  message(paste(cur_genome, "completed"))
  
}

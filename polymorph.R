# import snv data with minimal information on
# - genome: "UUJ150418_matabat2bin.303.fa"
# - scaffold
# - position
# - con_base: "A" / "C" / "T" / "G"
# - A, C, T, G: numerical (0 - ∞)
# - sample: "UUJ150318"
snvs_df <- read_csv("data/snvs.csv")

#alt1 <- read_csv("data/alt1.csv", col_names = FALSE)
#alt2 <- read_csv("data/alt2.csv", col_names = FALSE)
#snvs_df <- rbind(alt1, alt2)
#names(snvs_df) <- c("scaffold",
#                    "position",
#                    "position_coverage",
#                    "allele_count",
#                    "con_base",
#                    "A","C","T","G",
#                    "mutation_type",
#                    "genome",
#                    "sample")
  
# select only relevant columns
snvs_df <- snvs_df %>%
  transmute(genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)),
            scaffold,
            position,
            con_base,
            A, C, T, G,
            sample)

# iterate over all MAGs of interest
for (cur_genome in selected_genomes) {

  message(paste("analysing genome", cur_genome))
  
  # generate data.frame of the consensus bases
  # of the current MAG across the samples
  consensus_base_df <- snvs_df %>%
    filter(genome == cur_genome) %>%
    select(scaffold,
           position,
           con_base,
           sample) %>%
    pivot_wider(names_from = sample,
                values_from = con_base,
                values_fill = NA)
  
  # generate a matrix of the consensus bases
  # where each row is an SNV and each column
  # is a sample, so that  each cell tells you
  # in sample "x" which consensus base was found
  # for SNV "y"
  consensus_base_mat <- consensus_base_df %>%
    mutate(name = paste(scaffold, position, sep = ";")) %>%
    column_to_rownames(var = "name") %>%
    select(ends_with("UUM"),
           ends_with("UUJ"),
           ends_with("UUF")) %>%
    as.matrix()
  
  consensus_base_mat <- consensus_base_mat[ , order(colnames(consensus_base_mat))]
  count_list <- list()
  
  for (cur_sample in 1:ncol(consensus_base_mat)) {
    
    message(paste("computing sample", cur_sample))
    
    sample_mat <- consensus_base_mat[!is.na(consensus_base_mat[ , cur_sample]), ]
    
    count_list[[cur_sample]] <- sapply(iter(sample_mat, by = "col"),
                                       function(x) sum(sample_mat[ , cur_sample] == x, na.rm = TRUE))
    
  }
  
  count_df <- count_list %>%
    melt() %>%
    transmute(RefSample = colnames(consensus_base_mat)[L1],
              Count = value)
  
  count_df$Time <- as.factor(rep(colnames(consensus_base_mat), ncol(consensus_base_mat)))
  
  count_df <- count_df %>%
    group_by(RefSample) %>%
    summarise(Count,
              Time,
              Percent = Count / max(Count))
  
  count_df <- count_df %>%
    mutate(Time = gsub("UU\\w", "", Time)) %>%
    mutate(Time = paste(substr(Time, 1, 2), substr(Time, 3, 4), substr(Time, 5, 6), sep = "-")) %>%
    mutate(Time = as.Date(Time))
  
  winter <- substr(count_df$Time, 6, 7) %in% c("12", "01", "02")
  summer <- substr(count_df$Time, 6, 7) %in% c("06", "07", "08")
  spring <- substr(count_df$Time, 6, 7) %in% c("03", "04", "05")
  autumn <- substr(count_df$Time, 6, 7) %in% c("09", "10", "11")
  
  count_df$Season <- ""
  count_df$Season[winter] <- "Winter"
  count_df$Season[summer] <- "Summer"
  count_df$Season[spring] <- "Spring"
  count_df$Season[autumn] <- "Autumn"
  
  message("generating plot 1")
  
  p1 <- ggplot(count_df, aes(x = Time,
                             y = Percent,
                             colour = Season,
                             group = RefSample)) +
    geom_point() +
    geom_line() +
    scale_x_date(date_breaks = "1 month") +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.1)) +
    labs(y = "Percent of Shared Polymorphic Sites",
         colour = "Reference Sample") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  
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
    geom_bar() +
    scale_x_reverse(limits = c(43, 0),
                    breaks = seq(42, 1, by = -1)) +
    labs(x = "Consensus Base Recurrence",
         y = "Number of Shared Polymorphic Sites") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  message("plot 2 done")
  
  genome_time_series_df <- time_series_df %>%
    filter(genome == cur_genome)
  up_lim <- ceiling(max(genome_time_series_df$RPKM) / 10) * 10
  
  message("generating plot 3")
  
  p3 <- ggplot(genome_time_series_df, aes(x = Sample)) +
    geom_line(aes(y = RPKM), colour = "Black") +
    geom_line(aes(y = NucDiv * 1000), colour = "Dark Gray") +
    scale_x_date(date_breaks = "1 month") +
    scale_y_continuous(name = "RPKM",
                       breaks = seq(0, up_lim, by = 10),
                       sec.axis = sec_axis(~./1000,
                                           name = "Nucleotide Diversity")) +
    theme_bw() +
    labs(x = "Time") +
    theme(axis.title.y = element_text(colour = "Black"),
          axis.title.y.right = element_text(colour = "Dark Gray"),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
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

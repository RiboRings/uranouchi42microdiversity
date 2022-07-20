# import snv data with minimal information on
# - genome: "UUJ150418_matabat2bin.303.fa"
# - scaffold
# - position
# - con_base: "A" / "C" / "T" / "G"
# - A, C, T, G: numerical (0 - ∞)
# - sample: "UUJ150318"
snvs_df <- read_csv("data/snvs.csv")
  
# select only relevant columns
snvs_df <- snvs_df %>%
  transmute(genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)),
            scaffold,
            position,
            con_base,
            A, C, T, G,
            sample)

coverage_correction <- FALSE

# iterate over all MAGs of interest
for (cur_genome in selected_genomes) {

  message(paste("analysing genome", cur_genome))

  tax_df <- top_df %>%
    filter(genome == cur_genome) %>%
    transmute(Tax = paste(Phylum, Order, Family, sep = ";"))
  
  genome_df <- snvs_df %>%
    filter(genome == cur_genome) %>%
    select(scaffold,
           position,
           con_base,
           sample)
  
  # skip MAGs that appear in less than 2 samples
  if (length(unique(genome_df$sample)) < 2) {
    
    message(paste("skipping", cur_genome))
    next
    
  }
  
  # generate data.frame of the consensus bases
  # of the current MAG across the samples
  consensus_base_df <- genome_df %>%
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
  
  # sort columns so that samples are in
  # chronological order
  consensus_base_mat <- consensus_base_mat[ , order(colnames(consensus_base_mat))]
  count_list <- list()
  
  cov_df <- top_df %>%
    filter(genome == cur_genome) %>%
    select(starts_with("coverage")) %>%
    melt() %>%
    mutate(variable = gsub("coverage_", "", variable)) %>%
    filter(variable %in% colnames(consensus_base_mat)) %>%
    arrange(variable)
  
  low_coverage_samples <- cov_df$value < 10
  snv_likelyhood <- sapply(iter(consensus_base_mat, by = "row"), function(x) sum(!is.na(x)) / length(x))
  
  if (coverage_correction) {
    
    for (i in 1:nrow(consensus_base_mat)) {
      
      consensus_base_mat[i, is.na(consensus_base_mat[i, ]) & low_coverage_samples & snv_likelyhood[i] > 0.7]  <- max(consensus_base_mat[i, ], na.rm = TRUE)
      
    }
    
  }
  
  # compute number of shared polymorphic sites
  # each iteration taking a specific sample
  # as reference and comparing it to all the others
  for (cur_sample in 1:ncol(consensus_base_mat)) {
    
    message(paste("computing sample", cur_sample))
    
    sample_mat <- consensus_base_mat[!is.na(consensus_base_mat[ , cur_sample]), ]
    
    count_list[[cur_sample]] <- sapply(iter(sample_mat, by = "col"),
                                       function(x) sum(sample_mat[ , cur_sample] == x, na.rm = TRUE))
    
  }
  
  # melt counts of shared polymorphic sites
  # into 2 columns:
  # - L1 (reference sample identifier)
  # - value (number of shared polymorphic sites)
  count_df <- count_list %>%
    melt() %>%
    transmute(RefSample = colnames(consensus_base_mat)[L1],
              Count = value)
  
  # add column with information on sampling time
  count_df$ObservedSample <- as.factor(rep(colnames(consensus_base_mat), ncol(consensus_base_mat)))
  
  # compute percentage of shared polymorphic sites
  count_df <- count_df %>%
    group_by(RefSample) %>%
    summarise(Count,
              ObservedSample,
              Percent = Count / max(Count))
  
  # convert time to a date object type
  count_df <- count_df %>%
    mutate(Time = gsub("UU\\w", "", ObservedSample)) %>%
    mutate(Time = as.Date(paste(substr(Time, 1, 2),
                                substr(Time, 3, 4),
                                substr(Time, 5, 6),
                                sep = "-")))
  
  # index by season
  winter17 <- substr(count_df$RefSample, 3, 4) %in% c("12", "01", "02") &
    substr(count_df$RefSample, 1, 2) == "17"
  summer17 <- substr(count_df$RefSample, 3, 4) %in% c("06", "07", "08") &
    substr(count_df$RefSample, 1, 2) == "17"
  spring17 <- substr(count_df$RefSample, 3, 4) %in% c("03", "04", "05") &
    substr(count_df$RefSample, 1, 2) == "17"
  autumn17 <- substr(count_df$RefSample, 3, 4) %in% c("09", "10", "11") &
    substr(count_df$RefSample, 1, 2) == "17"
  winter18 <- substr(count_df$RefSample, 3, 4) %in% c("12", "01", "02") &
    substr(count_df$RefSample, 1, 2) == "18"
  summer18 <- substr(count_df$RefSample, 3, 4) %in% c("06", "07", "08") &
    substr(count_df$RefSample, 1, 2) == "18"
  spring18 <- substr(count_df$RefSample, 3, 4) %in% c("03", "04", "05") &
    substr(count_df$RefSample, 1, 2) == "18"
  autumn18 <- substr(count_df$RefSample, 3, 4) %in% c("09", "10", "11") &
    substr(count_df$RefSample, 1, 2) == "18"
  
  # store seasonality information into
  # a new column
  count_df$Season <- ""
  count_df$Season[winter17] <- "Winter 2017"
  count_df$Season[summer17] <- "Summer 2017"
  count_df$Season[spring17] <- "Spring 2017"
  count_df$Season[autumn17] <- "Autumn 2017"
  count_df$Season[winter18] <- "Winter 2018"
  count_df$Season[summer18] <- "Summer 2018"
  count_df$Season[spring18] <- "Spring 2018"
  count_df$Season[autumn18] <- "Autumn 2018"
  
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
  
  bar_array <- c()
  for (i in 1:nrow(consensus_base_mat)) {
    
    message(paste("computing recurrence", i))
    
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
  
  count_mat <-  count_df %>%
    select(-Count, -Time) %>%
    pivot_wider(values_from = Percent,
                names_from = ObservedSample) %>%
    column_to_rownames(var = "RefSample")
  
  set.seed(1234)
  row_ha <- HeatmapAnnotation(Season = gsub(" 201\\d", "", count_mat$Season),
                              which = "row")
  
  count_mat <- count_mat %>%
    select(-Season) %>%
    as.matrix()
  
  winter <- substr(colnames(count_mat), 3, 4) %in% c("12", "01", "02")
  summer <- substr(colnames(count_mat), 3, 4) %in% c("06", "07", "08")
  spring <- substr(colnames(count_mat), 3, 4) %in% c("03", "04", "05")
  autumn <- substr(colnames(count_mat), 3, 4) %in% c("09", "10", "11")
  
  col_ha_df_bottom <- data.frame(ObservedSample = colnames(count_mat))
  
  col_ha_df_bottom$Season <- ""
  col_ha_df_bottom$Season[winter] <- "Winter"
  col_ha_df_bottom$Season[summer] <- "Summer"
  col_ha_df_bottom$Season[spring] <- "Spring"
  col_ha_df_bottom$Season[autumn] <- "Autumn"
  
  set.seed(1234)
  col_ha_bottom <- HeatmapAnnotation(Season = col_ha_df_bottom$Season,
                                     which = "col",
                                     show_legend = FALSE)
  
  col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  
  g_values <- 10 * sapply(1:ncol(consensus_base_mat), function(x) mean(snv_likelyhood[is.na(consensus_base_mat[ , x])]) / cov_df$value[x])
  
  sig_threshold1 <- g_values < 0.05
  sig_threshold2 <- g_values < 0.005
  
  pch <- rep(NA, ncol(consensus_base_mat))
  pch[sig_threshold1] <- "*"
  pch[sig_threshold2] <- "**"
  
  gvalue_col_fun <- colorRamp2(c(0, 1.5, 3), c("green", "white", "red")) 
  col_ha_top <- HeatmapAnnotation(Trust = anno_simple(-log10(g_values), col = gvalue_col_fun, pch = pch))
  lgd_gvalue <- Legend(title = "Trust", col_fun = gvalue_col_fun, at = c(0, 1, 2, 3), 
                       labels = c("1", "0.1", "0.01", "0.001"))
  lgd_sig <- Legend(pch = c("*", "**"), type = "points", labels = c("< 0.05", "< 0.005"))
  
  # report results
  rmarkdown::render("magreport.Rmd",
                    output_format = "html_document",
                    output_file = paste0("results/", cur_genome, "_report.html"))

}

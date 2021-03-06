plot_distribution <- function(cur_stat) {
  
  ggplot(df, aes(x = eval(parse(text = cur_stat)))) +
    geom_histogram(bins = 40) +
    labs(x = cur_stat,
         y = "Count") +
    theme_bw()
  
}

freq_mapper <- function(cur_pos) {
  
  freq_array <- cur_pos %>%
    select(contains(paste0("Freq", cur_pos$con_base))) %>%
    as.vector() %>%
    as.numeric()
  
  return(freq_array)
  
}

max_mapper <- function(cur_pos, pattern) {
  
  max_val <- cur_pos %>%
    select(contains(paste0(pattern, "_", cur_pos$AbundMaxIdx, "."))) %>%
    as.vector() %>%
    as.numeric()
  
  return(max_val)
  
}

sample_loader <- function(sample_name) {
  
  sample_df <- read_csv2(sample_name) %>%
    mutate(sample = gsub(".csv", "", gsub("data/", "", sample_name)),
           genome = gsub("matabat2bin.", "", gsub(".fa", "", genome)))
  
  return(sample_df)
  
}

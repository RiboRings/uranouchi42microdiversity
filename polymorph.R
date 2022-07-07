# import snv data with minimal information on
# - genome: "UUJ150418_matabat2bin.303.fa"
# - scaffold
# - position
# - con_base: "A" / "C" / "T" / "G"
# - A, C, T, G: numerical (0 - ∞)
# - sample: "UUJ150318"
snvs_df <- read_csv("data/snvs.csv")

#snvs_df <- read_csv("data/alt1.csv", col_names = FALSE)
#snvs_df <- read_csv("data/alt2.csv", col_names = FALSE)
#snvs_df <- read_csv("data/dioamico.csv", col_names = FALSE)
#snvs_df <- rbind(alt1, alt2)
#
#names(snvs_df) <- c("scaffold",
#                    "position",
#                    "position_coverage",
#                    "allele_count",
#                    "con_base",
#                    "A", "C", "T", "G",
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

  message(paste("Analysing genome", cur_genome))

  tax_df <- top_df %>%
    filter(genome == cur_genome) %>%
    transmute(Tax = paste(Phylum, Order, Family, sep = ";"))
  
  # report results
  rmarkdown::render("magreport.Rmd",
                    output_format = "html_document",
                    output_file = paste0("results/", cur_genome, "_report.html"))

}
  

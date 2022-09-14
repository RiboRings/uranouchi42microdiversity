library(tidyr)
library(tibble)
library(dplyr)
library(readr)
library(readxl)
library(stringr)

abc <- gsub("UU", "", abc)
abc <- abc[!str_detect(abc, "F")]
abc <- gsub("M", "", abc)
abc <- gsub("J", "", abc)

env_fac <- read_excel("~/Desktop/publish/supplemental_material/florian_data_supplementary_dataset_1.xlsx", sheet = 2)
ef <- env_fac %>%
  filter(station %in% c("J", "M") & year %in% c("2017", "2018") & depth_m == 5 & date %in% abc)

ef_means <- ef %>%
  select(-`Visibility [m]`, -depth_m, -category, -station) %>%
  select(`Temperature [°C]`:`Silicium [µmol/L]`) %>%
  colMeans(na.rm = TRUE)

ef_stdevs <- ef %>%
  select(-`Visibility [m]`, -depth_m, -category, -station) %>%
  select(`Temperature [°C]`:`Silicium [µmol/L]`) %>%
  as.matrix() %>%
  colSds(na.rm = TRUE)

ef_medians <- ef %>%
  select(-`Visibility [m]`, -depth_m, -category, -station) %>%
  select(`Temperature [°C]`:`Silicium [µmol/L]`) %>%
  as.matrix() %>%
  colMedians(na.rm = TRUE)

ef_quantiles <- ef %>%
  select(-`Visibility [m]`, -depth_m, -category, -station) %>%
  select(`Temperature [°C]`:`Silicium [µmol/L]`) %>%
  as.matrix() %>%
  colQuantiles(na.rm = TRUE)

ef_stats <- data.frame(Mean = ef_means,
                       StdDev = ef_stdevs,
                       Quantile = ef_quantiles) %>%
  rownames_to_column(var = "Statistics")

ef_stats$Statistics[ef_stats$Statistics == "total_depth_m"] <- "total_depth [m]"
ef_stats$Statistics[ef_stats$Statistics == "odo_conc_[percent]"] <- "odo_conc [percent]"

ef_stats <- ef_stats %>%
  separate(col = Statistics,
           into = c("Statistics", "Unit"),
           sep = "\\[") %>%
  mutate(Unit = gsub("\\]", "", Unit),
         Statistics = gsub(" ", "", Statistics))

write_csv2(ef_stats, "~/Desktop/publish/figures/annex/sample_stats.csv")

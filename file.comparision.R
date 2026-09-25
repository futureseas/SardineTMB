install.packages("diffobj")

load.lib <- c("tidyverse", "data.table", "cowplot", "sf", 
             "here", "sdmTMB", "sdmTMBextra", "rnaturalearth", "rnaturalearthhires", 
             "mregions2" # Pulls official maritime EEZ records from MarineRegions.org
) 

install.lib <- load.lib[!load.lib %in% installed.packages()]

for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

sf_use_s2(FALSE)

dir <- here(getwd())
diffobj::diffFile(here("oldIPES.and.Canadian.Data.V3.Rmd", "IPES.and.Canadian.Data.V3.Rmd"))
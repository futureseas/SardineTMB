library(diffobj)
library(here)

dir <- getwd()

diffobj::diffFile(
  here(dir, "SardineTMB", "IPES.and.Canadian.Data.V3.Rmd"),
  here("Sardine_SptlOM", "IPES.and.Canadian.Data.V3-Original.Rmd")
  )

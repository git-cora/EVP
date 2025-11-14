install.packages(c("duckdb","DBI","httr","readr","dplyr","stringr"))

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(httr)
  library(readr)
  library(dplyr)
  library(stringr)
})

# --------Block 1:configure a single plant selection here ---------
source("extract_plant_compounds.R")

plant_species <- "Digitalis purpurea"   # change this once to reuse everywhere

# normalizes inputs
plant_slug <- str_replace_all(str_to_lower(plant_species), "[^a-z0-9]+", "_") |> str_replace_all("_{2,}", "_") |> str_replace("^_+|_+$", "")

run_coconut_pipeline(plant_species, out = paste0("coconut_", plant_slug, ".csv"))

#Optional: Reuse an existing DB without downloading again:
run_coconut_pipeline(plant_species, refresh = FALSE)

#---------------end----------------------------------------


#---------Block 2: ignition code for running the wrapper--------
source("run_pipeline_expval.R")

coconut_file  <- "coconut_digitalis_purpurea.csv"  # update when your extract file changes

options(
  expval.plant_name = plant_species,
  expval.coconut_csv = coconut_file
)

# simplest run
run_pipeline_expval(
  plant_name  = plant_species,            # pulled from the single configuration above
  coconut_csv = coconut_file              # pulled from the single configuration above
)

# ideal: send outputs to a folder and prefix file names
run_pipeline_expval(
  plant_name  = plant_species,
  coconut_csv = coconut_file,
  out_dir     = "D. purpurea",
  prefix      = "D. purp"
)

#---------------end----------------------------------

# Build a tiny local database from the COCONUT "Full CSV" snapshot,
# then export all molecules associated with a given organism 
# Usage (in R):
#   source("coconut_pipeline.R")
#   run_coconut_pipeline("plant_name")  # creates coconut_plant_name_<date>.csv
#
# Requires: duckdb, DBI, httr, readr, dplyr, stringr
# Install once: install.packages(c("duckdb","DBI","httr","readr","dplyr","stringr"))

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(httr)
  library(readr)
  library(dplyr)
  library(stringr)
})

# ---- Helpers ---------------------------------------------------------------

# Try several likely URLs for the monthly "Full CSV" ZIP. Falls back to file picker.
download_full_csv_zip <- function(dest_dir = tempdir(), months_back = 2) {
  # Build candidate URLs for current month and previous 'months_back' months.
  today <- Sys.Date()
  yms <- seq(from = as.Date(format(today, "%Y-%m-01")),
             by   = "-1 month",
             length.out = months_back + 1)
  fmt <- function(d, pat) sprintf(pat, format(d, "%Y"), format(d, "%m"))
  monthly <- as.vector(sapply(yms, function(d) c(
    fmt(d, "https://coconut.s3.uni-jena.de/prod/downloads/%1$s-%2$s/coconut_csv-%2$s-%1$s.zip"),
    fmt(d, "https://coconut.s3.uni-jena.de/coconut_csv-%2$s-%1$s.zip")
  )))
  candidates <- unique(c(
    monthly,
    "https://coconut.s3.uni-jena.de/prod/downloads/coconut_csv-latest.zip"
  ))

  message("Attempting automatic download of COCONUT Full CSV…")
  for (u in candidates) {
    message("  Trying: ", u)
    tf <- file.path(dest_dir, paste0("coconut_full_", basename(u)))
    resp <- try(RETRY("GET", u, timeout(120)), silent = TRUE)
    if (!inherits(resp, "try-error") && inherits(resp, "response") && status_code(resp) == 200) {
      writeBin(content(resp, "raw"), tf)
      message("  Downloaded: ", tf)
      return(tf)
    } else {
      message("  …nope (status ", tryCatch(status_code(resp), error=function(e) NA), ")")
    }
  }
  NA_character_
}

# Inside the ZIP, pick the largest CSV (usually the main table)
pick_main_csv <- function(zip_path) {
  zinfo <- unzip(zip_path, list = TRUE)
  csvs <- subset(zinfo, grepl("\\.csv$", Name, ignore.case = TRUE))
  if (nrow(csvs) == 0) stop("No .csv found in zip: ", zip_path)
  csv_name <- csvs$Name[which.max(csvs$Length)]
  csv_name
}

# Create/refresh a DuckDB from the Full CSV
build_duckdb_from_zip <- function(zip_path, db_path = "coconut.duckdb", table = "coconut") {
  stopifnot(file.exists(zip_path))
  csv_inside <- pick_main_csv(zip_path)
  csv_file <- unzip(zip_path, files = csv_inside, exdir = tempdir(), overwrite = TRUE)
  message("Using CSV: ", csv_file)

  con <- dbConnect(duckdb(db_path))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Import via DuckDB's fast CSV reader
  csv_norm <- normalizePath(csv_file, winslash = "/", mustWork = TRUE)
  sql <- paste0(
    "CREATE OR REPLACE TABLE ", DBI::dbQuoteIdentifier(con, table),
    " AS SELECT * FROM read_csv_auto(", DBI::dbQuoteString(con, csv_norm), ");"
  )
  dbExecute(con, sql)
  message("DuckDB table '", table, "' ready in: ", normalizePath(db_path, winslash="/"))
  invisible(db_path)
}

# Pick which column contains organism names
detect_organisms_col <- function(con, table = "coconut") {
  cols <- dbGetQuery(con, paste0("PRAGMA table_info(", dbQuoteIdentifier(con, table), ")"))$name
  candidates <- c("organisms", "organism_names", "organism", "species", "taxa")
  hit <- intersect(candidates, cols)
  if (length(hit) == 0) stop("Could not find an organisms column in table '", table, "'.")
  hit[1]
}

# Export all rows matching organism name to CSV
export_organism_csv <- function(con, organism, table = "coconut",
                                out_csv = NULL) {
  org_col <- detect_organisms_col(con, table)
  org_pat <- paste0("%", organism, "%")
  if (is.null(out_csv) || !nzchar(out_csv)) {
    safe <- gsub("[^A-Za-z0-9]+", "_", organism)
    out_csv <- sprintf("coconut_%s_%s.csv", safe, format(Sys.Date(), "%Y%m%d"))
  }
  qry <- paste0(
    "SELECT * FROM ", dbQuoteIdentifier(con, table),
    " WHERE ", dbQuoteIdentifier(con, org_col), " ILIKE ", dbQuoteString(con, org_pat), ";"
  )
  df <- dbGetQuery(con, qry)
  readr::write_csv(df, out_csv)
  message("Saved ", nrow(df), " rows to ", normalizePath(out_csv, winslash="/"))
  invisible(out_csv)
}

# ---- One-call pipeline ------------------------------------------------------
# 1) Download (or pick) Full CSV ZIP
# 2) Build/refresh DuckDB table
# 3) Export organism-specific CSV
run_coconut_pipeline <- function(organism = plant_species,
                                 db_path = "coconut.duckdb",
                                 table = "coconut",
                                 refresh = TRUE,
                                 zip_path = NULL,
                                 months_back = 2,
                                 out_csv = NULL) {
  if (refresh) {
    if (is.null(zip_path) || !file.exists(zip_path)) {
      zip_path <- download_full_csv_zip(months_back = months_back)
      if (is.na(zip_path)) {
        message("\nAutomatic download failed. Please select the COCONUT Full CSV ZIP manually.")
        if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
          zip_path <- rstudioapi::selectFile(caption = "Select COCONUT Full CSV .zip", filter = "zip")
        } else {
          zip_path <- file.choose()
        }
      }
    }
    build_duckdb_from_zip(zip_path, db_path = db_path, table = table)
  } else if (!file.exists(db_path)) {
    stop("DuckDB file '", db_path, "' not found. Run with refresh=TRUE first.")
  }
  con <- dbConnect(duckdb(db_path))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  export_organism_csv(con, organism, table = table, out_csv = out_csv)
}

# ---- End --------------------------------------------------------------------

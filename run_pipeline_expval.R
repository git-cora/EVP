# run_pipeline_expval.R
# Tiny wrapper for ExpValPipeline.R using the new naming scheme:
# products.csv, products_binned.csv, targets.csv, products_omitted.csv, metadata_*.md/json


if (!exists("default_plant_name", mode = "function")) {
  default_plant_name <- function(){
    getOption("expval.plant_name", Sys.getenv("EXPVAL_PLANT_NAME", "Coffea arabica"))
  }
}

if (!exists(".expval_slug", mode = "function")) {
  .expval_slug <- function(x){
    slug <- tolower(trimws(x))
    slug <- gsub("[^a-z0-9]+", "_", slug)
    slug <- gsub("_+", "_", slug)
    gsub("^_+|_+$", "", slug)
  }
}

if (!exists("default_coconut_csv", mode = "function")) {
  default_coconut_csv <- function(){
    default_file <- paste0("coconut_", .expval_slug(default_plant_name()), ".csv")
    getOption("expval.coconut_csv", Sys.getenv("EXPVAL_COCONUT_CSV", default_file))
  }
}

run_pipeline_expval <- function(
    plant_name = default_plant_name(),
    coconut_csv = default_coconut_csv(),
    out_dir = ".",
    prefix = NULL,
    enrich_mechanism = TRUE
){
  source("ExpValPipeline.R")
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  name <- function(stem){
    if (is.null(prefix) || !nzchar(prefix)) file.path(out_dir, stem) else file.path(out_dir, paste0(prefix, "_", stem))
  }
  
  products_path   <- name("products.csv")
  products_binned <- name("products_binned.csv")
  targets_path    <- name("targets.csv")
  omitted_path    <- name("products_omitted.csv")
  metadata_md     <- name(paste0("metadata_", format(Sys.Date(), "%Y%m%d"), ".md"))
  metadata_json   <- name(paste0("metadata_", format(Sys.Date(), "%Y%m%d"), ".json"))
  
  # Stage 1: per-compound targets
  targets_v2_for_csv(
    in_csv        = coconut_csv,
    out_csv       = products_path,
    omitted_csv   = omitted_path,
    collapsed_csv = NULL,
    enrich_mechanism = enrich_mechanism
  )
  
  # Stage 2: bin + unique targets + report
  bin_targets_v2(
    in_csv            = products_path,
    out_binned        = products_binned,
    out_unique_binned = targets_path,
    report            = TRUE,
    plant_name        = plant_name,
    coconut_csv       = coconut_csv,
    omitted_csv       = omitted_path,
    report_md         = metadata_md,
    report_json       = metadata_json
  )
  
  message("Done. Outputs in: ", normalizePath(out_dir))
  invisible(list(
    products = products_path,
    products_binned = products_binned,
    targets = targets_path,
    omitted = omitted_path,
    metadata_md = metadata_md,
    metadata_json = metadata_json
  ))
}

# Auto-run if executed directly
if (sys.nframe() == 0) {
  run_pipeline_expval()
}

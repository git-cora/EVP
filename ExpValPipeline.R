# ExpValPipeline.R, now present in GitHub
# COCONUT -> ChEMBL with MoA/site + UniProt enrichment + auto metadata report
# Naming:
#   products.csv          — per-compound targets
#   products_binned.csv   — per-compound (binned)
#   targets.csv           — unique targets (collapsed by UniProt)
#   products_omitted.csv  — audit of excluded compounds
#   metadata_<plant>_<date>.md/json — run report
#
# Public entry points:
#   - targets_v2_for_csv()
#   - bin_targets_v2()
#   - write_run_report()

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr); library(purrr)
  library(tibble); library(stringr); library(readr); library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y
UA <- httr::user_agent("ExpValPipeline/1.0 (R httr)")
keep_types <- c("Ki","Kd","IC50","EC50","AC50")

.expval_slug <- function(x){
  x |>
    stringr::str_trim() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("_{2,}", "_") |>
    stringr::str_replace("^_+|_+$", "")
}

default_plant_name <- function(){
  getOption("expval.plant_name", Sys.getenv("EXPVAL_PLANT_NAME", "Coffea arabica"))
}

default_coconut_csv <- function(){
  default_file <- paste0("coconut_", .expval_slug(default_plant_name()), ".csv")
  getOption("expval.coconut_csv", Sys.getenv("EXPVAL_COCONUT_CSV", default_file))
}

set_expval_defaults <- function(plant_name = NULL, coconut_csv = NULL){
  if (!is.null(plant_name)) {
    options(expval.plant_name = plant_name)
  }
  if (is.null(coconut_csv) && !is.null(plant_name)) {
    coconut_csv <- paste0("coconut_", .expval_slug(plant_name), ".csv")
  }
  if (!is.null(coconut_csv)) {
    options(expval.coconut_csv = coconut_csv)
  }
  invisible(list(
    plant_name = default_plant_name(),
    coconut_csv = default_coconut_csv()
  ))
}

chembl_retry_get <- function(url, query = list(), timeout_s = 60, tries = 3){
  httr::RETRY("GET", url, UA, query = query,
              httr::timeout(timeout_s), httr::config(connecttimeout = 20), httr::config(low_speed_time = 0),
              times = tries, pause_base = 1, pause_min = 1, pause_cap = 8)
}

chembl_mols_from_inchikey <- function(inchikey, limit = 3){
  r <- chembl_retry_get("https://www.ebi.ac.uk/chembl/api/data/molecule.json",
                        list(molecule_structures__standard_inchi_key = inchikey, limit = limit))
  if (http_error(r)) return(character(0))
  j <- content(r, as = "parsed", type = "application/json")
  if (is.null(j$molecules)) return(character(0))
  unique(vapply(j$molecules, function(m) m$molecule_chembl_id %||% NA_character_, "")) |>
    discard(~is.na(.x) || !nzchar(.x))
}

to_nM <- function(val, units){
  v <- suppressWarnings(as.numeric(val)); u <- stringr::str_to_lower(units %||% "")
  dplyr::case_when(
    is.na(v) ~ NA_real_,
    u %in% c("nm") ~ v,
    u %in% c("um","µm","μm","umol/l","uM") ~ v*1e3,
    u %in% c("mm") ~ v*1e6,
    TRUE ~ NA_real_
  )
}

chembl_targets_for_mol <- function(mol_id, page_limit = 200){
  base <- "https://www.ebi.ac.uk/chembl/api/data/activity.json"
  out <- list(); off <- 0; total <- Inf
  repeat {
    r <- chembl_retry_get(base, list(molecule_chembl_id = mol_id, limit = page_limit, offset = off), timeout_s = 90)
    if (http_error(r)) break
    j <- content(r, as = "parsed", type = "application/json")
    acts <- j$activities
    if (length(acts) == 0) break
    out <- c(out, acts); off <- off + page_limit
    total <- j$page_meta$total_count %||% total
    if (off >= total) break
    Sys.sleep(0.12)
  }
  if (!length(out)) return(tibble())
  tibble(
    target_id        = purrr::map_chr(out, ~ .x$target_chembl_id %||% NA_character_),
    target_name      = purrr::map_chr(out, ~ .x$target_pref_name %||% NA_character_),
    target_species   = purrr::map_chr(out, ~ .x$target_species %||% NA_character_),
    std_type         = purrr::map_chr(out, ~ .x$standard_type %||% NA_character_),
    std_value        = purrr::map_chr(out, ~ .x$standard_value %||% NA_character_),
    std_units        = purrr::map_chr(out, ~ .x$standard_units %||% NA_character_),
    ref_id           = purrr::map_chr(out, ~ .x$document_chembl_id %||% NA_character_),
    assay_chembl_id  = purrr::map_chr(out, ~ .x$assay_chembl_id %||% NA_character_)
  ) |>
    mutate(std_value_num = suppressWarnings(as.numeric(std_value)),
           potency_nM    = to_nM(std_value_num, std_units))
}

chembl_target_meta <- function(target_id){
  r <- chembl_retry_get("https://www.ebi.ac.uk/chembl/api/data/target.json",
                        list(target_chembl_id = target_id, limit = 1))
  if (http_error(r)) return(tibble(target_id = target_id, target_type = NA, pref_name = NA, uniprot = NA))
  j <- content(r, as = "parsed", type = "application/json")
  if (is.null(j$targets) || !length(j$targets)) return(tibble(target_id = target_id, target_type = NA, pref_name = NA, uniprot = NA))
  t <- j$targets[[1]]
  comps <- t$target_components %||% list()
  accs  <- unique(vapply(comps, function(tc) as.character(tc$accession %||% NA_character_), ""))
  accs  <- accs[nzchar(accs)]
  tibble(target_id = target_id,
         target_type = t$target_type %||% NA_character_,
         pref_name   = t$pref_name %||% NA_character_,
         uniprot     = if (length(accs)) paste(accs, collapse = ";") else NA_character_)
}

.mechanism_cache <- new.env(parent = emptyenv())
chembl_mechanisms_for_mol <- function(mol_id){
  key <- paste0("mech:", mol_id)
  if (exists(key, envir = .mechanism_cache)) return(get(key, envir = .mechanism_cache))
  r <- chembl_retry_get("https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
                        list(molecule_chembl_id = mol_id, limit = 200))
  if (http_error(r)) { assign(key, tibble(), envir = .mechanism_cache); return(tibble()) }
  j <- content(r, as = "parsed", type = "application/json")
  mechs <- j$mechanisms %||% list()
  df <- tibble::tibble()
  if (length(mechs)) {
    df <- tibble(
      mechanism_of_action = purrr::map_chr(mechs, ~ .x$mechanism_of_action %||% NA_character_),
      action_type         = purrr::map_chr(mechs, ~ .x$action_type %||% NA_character_),
      target_chembl_id    = purrr::map_chr(mechs, ~ .x$target_chembl_id %||% NA_character_),
      target_pref_name    = purrr::map_chr(mechs, ~ .x$target_pref_name %||% NA_character_)
    )
  }
  assign(key, df, envir = .mechanism_cache)
  df
}

.assay_cache <- new.env(parent = emptyenv())
chembl_assay_desc <- function(assay_chembl_id){
  if (!nzchar(assay_chembl_id)) return(NA_character_)
  key <- paste0("assay:", assay_chembl_id)
  if (exists(key, envir = .assay_cache)) return(get(key, envir = .assay_cache))
  r <- chembl_retry_get("https://www.ebi.ac.uk/chembl/api/data/assay.json",
                        list(assay_chembl_id = assay_chembl_id, limit = 1))
  if (http_error(r)) { assign(key, NA_character_, envir = .assay_cache); return(NA_character_) }
  j <- content(r, as = "parsed", type = "application/json")
  if (is.null(j$assays) || !length(j$assays)) { assign(key, NA_character_, envir = .assay_cache); return(NA_character_) }
  desc <- j$assays[[1]]$description %||% NA_character_
  assign(key, desc, envir = .assay_cache)
  desc
}

derive_moa_site <- function(mech_rows, assay_descs){
  moa <- NA_character_; site <- NA_character_
  if (nrow(mech_rows)) {
    at <- unique(na.omit(mech_rows$action_type))
    moakeys <- c("AGONIST","ANTAGONIST","INVERSE AGONIST","PARTIAL AGONIST","INHIBITOR","MODULATOR","POSITIVE ALLOSTERIC MODULATOR","NEGATIVE ALLOSTERIC MODULATOR")
    hit <- at[at %in% moakeys]
    if (length(hit)) moa <- hit[1]
    moa <- stringr::str_to_title(tolower(moa))
    mo_text <- paste(na.omit(mech_rows$mechanism_of_action), collapse = " ")
    if (nzchar(mo_text)) {
      lo <- tolower(mo_text)
      if (grepl("allosteric", lo)) site <- "Allosteric"
      if (grepl("orthosteric", lo)) site <- "Orthosteric"
    }
  }
  if (is.na(site) || !nzchar(site)) {
    d <- paste(na.omit(assay_descs), collapse = " ")
    lo <- tolower(d)
    if (grepl("allosteric", lo)) site <- "Allosteric"
    if (grepl("orthosteric|competitive", lo)) site <- "Orthosteric"
    if (is.na(site) || !nzchar(site)) site <- "Unknown"
  }
  if (is.na(moa) || !nzchar(moa)) {
    d <- paste(na.omit(assay_descs), collapse = " ")
    lo <- tolower(d)
    if (grepl("agonist|partial agonist|inverse agonist", lo)) moa <- "Agonist/variant"
    else if (grepl("antagonist|blocker", lo)) moa <- "Antagonist"
    else if (grepl("inhibitor", lo)) moa <- "Inhibitor"
    else if (grepl("modulator|pam|nam", lo)) moa <- "Modulator"
    else moa <- "Unknown"
  }
  tibble(moa_action = moa, site_annotation = site)
}

# ----------------- UniProt enrichment helpers -----------------
.uniprot_cache <- new.env(parent = emptyenv())

.ensure_uniprot_cols <- function(df) {
  must <- c("accession","protein_name","gene_names","keyword","go_molecular_function","ec_number")
  for (m in must) if (!m %in% names(df)) df[[m]] <- NA_character_
  df
}

.uniprot_fetch_batch <- function(accs){
  accs <- unique(accs[!is.na(accs) & nzchar(accs)])
  if (!length(accs)) return(tibble())
  q <- paste0("(", paste(sprintf("accession:%s", accs), collapse = " OR "), ")")
  url <- "https://rest.uniprot.org/uniprotkb/search"
  r <- httr::RETRY(
    "GET", url,
    query = list(
      query  = q,
      fields = paste(c("accession","protein_name","gene_names","keyword","go_f","ec"), collapse=","),
      format = "tsv",
      size   = 500
    ),
    times = 3, pause_base = 1, pause_cap = 8
  )
  if (httr::http_error(r)) return(tibble())
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  if (!nzchar(txt)) return(tibble())
  con <- textConnection(txt); on.exit(close(con), add = TRUE)
  df <- tryCatch(readr::read_tsv(con, show_col_types = FALSE), error = function(e) tibble())
  nms <- tolower(names(df)); if (length(nms)) names(df) <- nms
  names(df) <- gsub(" ", "_", names(df), fixed = TRUE)
  names(df) <- gsub("\\(|\\)", "", names(df))
  df <- df %>% dplyr::rename(
    accession = dplyr::any_of("accession"),
    protein_name = dplyr::any_of("protein_name"),
    gene_names = dplyr::any_of("gene_names"),
    keyword = dplyr::any_of("keyword"),
    go_molecular_function = dplyr::any_of("go_molecular_function"),
    ec_number = dplyr::any_of("ec_number")
  )
  .ensure_uniprot_cols(df)
}

# PATCH 1: make UniProt fetch less brittle (smaller batches = fewer timeouts)
uniprot_fetch_meta <- function(accs, batch = 60){
  accs <- unique(accs[!is.na(accs) & nzchar(accs)])
  if (!length(accs)) return(tibble())
  to_get <- accs[!vapply(accs, function(a) exists(a, envir = .uniprot_cache), logical(1))]
  if (length(to_get)) {
    for (chunk in split(to_get, ceiling(seq_along(to_get)/batch))) {
      df <- .uniprot_fetch_batch(chunk)
      if (nrow(df)) {
        for (i in seq_len(nrow(df))) {
          a <- as.character(df$accession[i])
          assign(a, df[i, , drop = FALSE], envir = .uniprot_cache)
        }
      }
      missing <- setdiff(chunk, vapply(ls(.uniprot_cache), identity, character(1)))
      if (length(missing)) for (m in missing) assign(m, tibble(accession = m), envir = .uniprot_cache)
      Sys.sleep(0.2)
    }
  }
  rows <- lapply(accs, function(a) get(a, envir = .uniprot_cache))
  dplyr::bind_rows(rows)
}

class_from_uniprot_vec <- function(keyword, go_molecular_function, ec_number){
  kw <- tolower(as.character(keyword %||% ""))
  go <- tolower(as.character(go_molecular_function %||% ""))
  ec <- tolower(as.character(ec_number %||% ""))
  
  n <- max(length(kw), length(go), length(ec))
  kw <- rep_len(kw, n); go <- rep_len(go, n); ec <- rep_len(ec, n)
  
  out <- rep("Other/unknown", n)
  
  gpcr_kw <- "g[ -]?protein[ -]?coupled receptor"
  gpcr_go <- paste0(gpcr_kw, " activity")
  gpcr <- str_detect(kw, gpcr_kw) |
    str_detect(go, gpcr_go) |
    str_detect(kw, "seven transmembrane receptor") |
    str_detect(kw, "class [ai]/?[0-9]? g protein coupled receptor")
  
  out[gpcr] <- "GPCR"
  
  kin <- str_detect(kw, "\\bkinase\\b") | str_detect(go, "protein kinase activity")
  out[kin & out == "Other/unknown"] <- "Kinase"
  
  ich <- str_detect(kw, "\\bion channel\\b") | str_detect(go, "ion channel activity") |
    str_detect(kw, "voltage-gated")
  out[ich & out == "Other/unknown"] <- "Ion channel"
  
  nrx <- str_detect(kw, "\\bnuclear receptor\\b")
  out[nrx & out == "Other/unknown"] <- "Nuclear receptor"
  
  trn <- str_detect(kw, "transporter") | str_detect(go, "transporter activity")
  out[trn & out == "Other/unknown"] <- "Transporter"
  
  enz <- nzchar(ec) | str_detect(kw, "phosphodiesterase|esterase|oxidase|transferase|carboxylase|synthetase|hydrolase|dehydrogenase|cyclooxygenase|lipoxygenase|protease|peptidase")
  out[enz & out == "Other/unknown"] <- "Enzyme"
  
  out
}

# -------------------- Stage 1: build per-compound products --------------------
targets_v2_for_csv <- function(in_csv,
                               out_csv = "products.csv",
                               omitted_csv = "products_omitted.csv",
                               collapsed_csv = NULL,
                               enrich_mechanism = TRUE) {
  comp <- readr::read_csv(in_csv, show_col_types = FALSE)
  if (!"inchikey" %in% names(comp)) {
    cand <- intersect(names(comp), c("standard_inchi_key","InChIKey","INCHIKEY"))
    if (length(cand)) comp <- dplyr::rename(comp, inchikey = !!cand[1]) else stop("No inchikey column found.")
  }
  comp <- comp |> dplyr::distinct(inchikey, .keep_all = TRUE) |> dplyr::filter(!is.na(inchikey) & nzchar(inchikey))
  if (!"name" %in% names(comp)) comp$name <- NA_character_
  
  audit <- vector("list", nrow(comp))
  rows  <- vector("list", nrow(comp))
  
  for (i in seq_len(nrow(comp))) {
    ik <- comp$inchikey[i]
    nm <- comp$name[i] %||% NA_character_
    message(sprintf("[%d/%d] %s", i, nrow(comp), ik))
    
    a <- list(compound_key = ik, compound_name = nm,
              had_chembl_ids = FALSE,
              pre_rows = 0L,
              post_binding = NA_integer_,
              post_human = NA_integer_,
              post_single_protein = NA_integer_,
              omission_reason = NA_character_)
    
    mols <- tryCatch(chembl_mols_from_inchikey(ik), error = function(e) character(0))
    a$had_chembl_ids <- length(mols) > 0
    if (!length(mols)) { a$omission_reason <- "no_chembl_match"; audit[[i]] <- a; rows[[i]] <- tibble(); next }
    
    acts_full <- dplyr::bind_rows(lapply(mols, chembl_targets_for_mol))
    a$pre_rows <- nrow(acts_full)
    if (!nrow(acts_full)) { a$omission_reason <- "no_chembl_activities"; audit[[i]] <- a; rows[[i]] <- tibble(); next }
    
    acts_bind <- acts_full |> dplyr::filter(std_type %in% keep_types)
    a$post_binding <- nrow(acts_bind)
    if (!nrow(acts_bind)) { a$omission_reason <- "no_binding_assays"; audit[[i]] <- a; rows[[i]] <- tibble(); next }
    
    acts_human <- acts_bind |> dplyr::filter(is.na(target_species) | target_species == "Homo sapiens")
    a$post_human <- nrow(acts_human)
    if (!nrow(acts_human)) { a$omission_reason <- "non_human_only"; audit[[i]] <- a; rows[[i]] <- tibble(); next }
    
    tids <- unique(acts_human$target_id)
    meta <- dplyr::bind_rows(lapply(tids, chembl_target_meta))
    keep <- meta |>
      dplyr::filter(!is.na(uniprot) & nzchar(uniprot)) |>
      dplyr::filter(!is.na(target_type) & grepl("SINGLE PROTEIN", target_type, ignore.case = TRUE))
    
    acts_keep <- dplyr::inner_join(acts_human, keep, by = "target_id")
    a$post_single_protein <- nrow(acts_keep)
    if (!nrow(acts_keep)) { a$omission_reason <- "non_single_protein_or_no_uniprot"; audit[[i]] <- a; rows[[i]] <- tibble(); next }
    
    out_i <- dplyr::bind_rows(lapply(mols, function(mid){
      df <- chembl_targets_for_mol(mid)
      if (!nrow(df)) return(tibble())
      df
    }))
    out_i <- out_i |>
      dplyr::filter(std_type %in% keep_types) |>
      dplyr::filter(is.na(target_species) | target_species == "Homo sapiens") |>
      dplyr::mutate(compound_key = ik, compound_name = nm)
    
    # condense with meta + mechanism
    keep <- dplyr::inner_join(out_i, keep, by = "target_id")
    if (!nrow(keep)) { audit[[i]] <- a; rows[[i]] <- tibble(); next }
    
    mech_all <- dplyr::bind_rows(lapply(mols, chembl_mechanisms_for_mol))
    assay_descs <- unique(keep$assay_chembl_id)
    assay_map <- tibble(assay_chembl_id = assay_descs,
                        description = vapply(assay_descs, chembl_assay_desc, ""))
    
    out_i <- keep |>
      dplyr::group_by(target_id, pref_name, uniprot, compound_key, compound_name) |>
      dplyr::summarise(
        best_type = std_type[which.min(potency_nM %||% Inf)][1] %||% NA_character_,
        best_nM   = suppressWarnings(min(potency_nM, na.rm = TRUE)),
        ref_id    = ref_id[which.min(potency_nM %||% Inf)][1] %||% NA_character_,
        assays    = list(unique(assay_chembl_id)),
        .groups = "drop"
      ) |>
      dplyr::mutate(assay_descs = list(assay_map$description[match(unlist(assays), assay_map$assay_chembl_id)])) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        moa_site = list(derive_moa_site(
          mech_rows = if (nrow(mech_all)) mech_all %>% dplyr::filter(target_chembl_id == target_id) else tibble(),
          assay_descs = unlist(assay_descs)
        ))
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(moa_action = vapply(moa_site, function(x) x$moa_action, ""),
                    site_annotation = vapply(moa_site, function(x) x$site_annotation, "")) |>
      dplyr::select(-assays, -assay_descs, -moa_site) |>
      dplyr::arrange(best_nM)
    
    audit[[i]] <- a
    rows[[i]]  <- out_i
  }
  
  out <- dplyr::bind_rows(rows)
  readr::write_csv(out, out_csv)
  
  audit_df <- dplyr::bind_rows(lapply(audit, tibble::as_tibble))
  included <- out %>% dplyr::distinct(compound_key) %>% dplyr::mutate(included = TRUE)
  audit_df <- audit_df %>%
    dplyr::left_join(included, by = "compound_key") %>%
    dplyr::mutate(included = dplyr::coalesce(included, FALSE)) %>%
    dplyr::select(compound_key, compound_name, had_chembl_ids, pre_rows, post_binding, post_human, post_single_protein, included, omission_reason)
  
  readr::write_csv(audit_df %>% dplyr::filter(!included), omitted_csv)
  
  message("Wrote: ", out_csv, " and ", omitted_csv)
  
  if (!is.null(collapsed_csv)) {
    collapse <- out %>%
      dplyr::filter(!is.na(uniprot) & nzchar(uniprot)) %>%
      dplyr::group_by(uniprot) %>%
      dplyr::summarise(
        target_name  = dplyr::first(na.omit(pref_name)),
        n_compounds  = dplyr::n_distinct(compound_key),
        best_nM      = suppressWarnings(min(best_nM, na.rm = TRUE)),
        best_measure = best_type[which.min(best_nM)][1],
        example_ref  = ref_id[which.min(best_nM)][1],
        .groups = "drop"
      ) %>%
      dplyr::arrange(best_nM)
    readr::write_csv(collapse, collapsed_csv)
    message("Wrote: ", collapsed_csv)
  }
  
  invisible(out)
}

# -------------------- Stage 2: bin products + emit targets + metadata ---------
bin_targets_v2 <- function(in_csv = "products.csv",
                           out_binned = "products_binned.csv",
                           out_unique_binned = "targets.csv",
                           report = TRUE,
                           plant_name = default_plant_name(),
                           coconut_csv = default_coconut_csv(),
                           omitted_csv = "products_omitted.csv",
                           report_md = NULL,
                           report_json = NULL){
  x <- readr::read_csv(in_csv, show_col_types = FALSE)
  
  if (!"moa_action" %in% names(x)) x$moa_action <- NA_character_
  if (!"site_annotation" %in% names(x)) x$site_annotation <- NA_character_
  
  x <- x %>% mutate(uniprot = as.character(uniprot))
  x_exp <- x %>%
    mutate(uniprot = strsplit(uniprot, ";")) %>%
    tidyr::unnest(uniprot) %>%
    mutate(uniprot = str_trim(uniprot))
  
  meta <- uniprot_fetch_meta(unique(x_exp$uniprot))
  
  # PATCH 2: softened fallback to avoid "everything is Enzyme"
  if (nrow(meta) == 0) {
    warning("UniProt enrichment failed; falling back to heuristic classes (soft).")
    nm <- tolower(x$pref_name %||% "")
    class_tag <- rep("Other/unknown", nrow(x))
    
    # high-confidence matches first
    class_tag[str_detect(nm, "g[ -]?protein[ -]?coupled receptor|\bgpcr\b|adenosine receptor|adrenergic receptor|dopamine receptor|serotonin receptor")] <- "GPCR"
    class_tag[str_detect(nm, " kinase")] <- "Kinase"
    class_tag[str_detect(nm, " channel")] <- "Ion channel"
    class_tag[str_detect(nm, " nuclear receptor| nr[0-9]")] <- "Nuclear receptor"
    class_tag[str_detect(nm, " transporter")] <- "Transporter"
    
    # only call Enzyme if nothing else matched
    enzyme_like <- str_detect(nm, " phosphodiesterase|esterase|oxidase|transferase|carboxylase|synthetase|hydrolase|dehydrogenase|cyclooxygenase|lipoxygenase|protease|peptidase")
    class_tag[enzyme_like & class_tag == "Other/unknown"] <- "Enzyme"
    
    x$target_class <- class_tag
    x$class_source <- "Heuristic/Unknown"
  } else {
    meta <- .ensure_uniprot_cols(meta)
    
    meta_small <- meta %>%
      dplyr::select(dplyr::any_of(c(
        "accession","protein_name","gene_names","keyword","go_molecular_function","ec_number"
      ))) %>%
      dplyr::rename(
        uniprot               = accession,
        protein_name_uniprot  = protein_name,
        gene_names_uniprot    = gene_names
      ) %>%
      mutate(across(c(keyword, go_molecular_function, ec_number, protein_name_uniprot, gene_names_uniprot), ~ as.character(.)))
    
    xj <- x_exp %>% left_join(meta_small, by = "uniprot") %>%
      mutate(
        target_class_up = class_from_uniprot_vec(keyword, go_molecular_function, ec_number),
        class_source_up = dplyr::if_else(target_class_up != "Other/unknown", "UniProt", "Heuristic/Unknown")
      )
    
    pick_class <- function(v){
      v <- unique(na.omit(v))
      if (!length(v)) return(NA_character_)
      good <- setdiff(v, c("Other/unknown", NA_character_))
      if (length(good)) return(good[1])
      v[1]
    }
    pick_src <- function(v){
      v <- unique(na.omit(v))
      if (!length(v)) return(NA_character_)
      if ("UniProt" %in% v) return("UniProt")
      v[1]
    }
    
    x <- xj %>%
      group_by(compound_key, compound_name, target_id, pref_name, uniprot) %>%
      summarise(
        best_type = dplyr::first(best_type),
        best_nM   = dplyr::first(best_nM),
        ref_id    = dplyr::first(ref_id),
        moa_action = dplyr::first(moa_action),
        site_annotation = dplyr::first(site_annotation),
        target_class = pick_class(target_class_up),
        class_source = pick_src(class_source_up),
        .groups = "drop"
      )
  }
  
  x <- x %>%
    mutate(potency_bin = case_when(
      !is.na(best_nM) & best_nM <= 10 ~ "Ultra-high (≤10 nM)",
      !is.na(best_nM) & best_nM <= 100 ~ "High (10–100 nM)",
      !is.na(best_nM) & best_nM <= 1000 ~ "Moderate (100–1000 nM)",
      !is.na(best_nM) & best_nM <= 10000 ~ "Weak (1–10 µM)",
      !is.na(best_nM) & best_nM > 10000 ~ "Very weak (>10 µM)",
      TRUE ~ "Unquantified"
    )) %>%
    mutate(assay_bin = case_when(
      best_type %in% c("Ki","Kd") ~ "Binding (Ki/Kd)",
      best_type %in% c("IC50","EC50","AC50") ~ "Functional (IC50/EC50/AC50)",
      TRUE ~ "Other/NA"
    ))
  
  readr::write_csv(x, out_binned)
  
  uniq <- x %>%
    filter(!is.na(uniprot) & nzchar(uniprot)) %>%
    group_by(uniprot) %>%
    summarise(
      target_name  = dplyr::first(na.omit(pref_name)),
      target_class = dplyr::first(na.omit(target_class)),
      n_compounds  = dplyr::n_distinct(compound_key),
      best_nM      = suppressWarnings(min(best_nM, na.rm = TRUE)),
      best_measure = best_type[which.min(best_nM)][1],
      .groups = "drop"
    ) %>%
    mutate(coverage_bin = case_when(
      n_compounds >= 5 ~ "Broadly hit (≥5)",
      n_compounds >= 3 ~ "Often hit (3–4)",
      n_compounds == 2 ~ "Occasionally hit (2)",
      n_compounds == 1 ~ "Single-hit (1)",
      TRUE ~ "NA"
    )) %>%
    arrange(best_nM)
  
  readr::write_csv(uniq, out_unique_binned)
  message("Wrote: ", out_binned, " and ", out_unique_binned)
  
  if (isTRUE(report)) {
    write_run_report(
      plant_name = plant_name,
      coconut_csv = coconut_csv,
      targets_csv = in_csv,
      binned_csv = out_binned,
      unique_binned_csv = out_unique_binned,
      omitted_csv = omitted_csv,
      out_md = report_md,
      out_json = report_json
    )
  }
  
  invisible(list(per_compound = x, unique_targets = uniq))
}

# -------------------- Metadata report -----------------------------------------
write_run_report <- function(
    plant_name =  default_plant_name(),
    coconut_csv = default_coconut_csv(),
    targets_csv = "products.csv",
    binned_csv = "products_binned.csv",
    unique_binned_csv = "targets.csv",
    omitted_csv = "products_omitted.csv",
    out_md = NULL,
    out_json = NULL
){
  suppressPackageStartupMessages({ library(readr); library(dplyr); library(stringr); library(tidyr); library(jsonlite) })
  
  read_safe <- function(path){
    if (is.null(path) || !file.exists(path)) return(NULL)
    tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  }
  
  coco   <- read_safe(coconut_csv)
  targ   <- read_safe(targets_csv)
  bin    <- read_safe(binned_csv)
  uniq   <- read_safe(unique_binned_csv)
  omit   <- read_safe(omitted_csv)
  
  total_compounds <- NA_integer_
  if (!is.null(coco)) {
    nm <- tolower(names(coco))
    ik_idx <- match("inchikey", nm, nomatch = NA)
    if (is.na(ik_idx)) ik_idx <- match("standard_inchi_key", nm, nomatch = NA)
    if (is.na(ik_idx)) ik_idx <- match("inchikey", names(coco), nomatch = NA)
    if (!is.na(ik_idx)) {
      ik_col <- names(coco)[ik_idx]
      total_compounds <- coco %>%
        dplyr::filter(!is.na(.data[[ik_col]]) & nzchar(.data[[ik_col]])) %>%
        dplyr::pull(dplyr::all_of(ik_col)) %>%
        dplyr::n_distinct()
    } else {
      total_compounds <- nrow(coco)
    }
  }
  
  included_compounds <- if (!is.null(targ)) targ %>% distinct(compound_key) %>% nrow() else NA_integer_
  omitted_compounds  <- if (!is.null(omit)) omit %>% distinct(compound_key) %>% nrow() else NA_integer_
  unique_targets     <- if (!is.null(uniq)) nrow(uniq) else if (!is.null(targ)) targ %>% filter(!is.na(uniprot) & nzchar(uniprot)) %>% distinct(uniprot) %>% nrow() else NA_integer_
  
  class_dist   <- if (!is.null(bin) && "target_class" %in% names(bin)) bin %>% count(target_class, sort = TRUE) else tibble(target_class = character(), n = integer())
  potency_dist <- if (!is.null(bin) && "potency_bin" %in% names(bin)) bin %>% count(potency_bin, sort = TRUE) else tibble(potency_bin = character(), n = integer())
  assay_dist   <- if (!is.null(bin) && "assay_bin" %in% names(bin)) bin %>% count(assay_bin, sort = TRUE) else tibble(assay_bin = character(), n = integer())
  moa_dist     <- if (!is.null(bin) && "moa_action" %in% names(bin)) bin %>% mutate(moa_action = ifelse(is.na(moa_action) | !nzchar(moa_action), "Unknown", moa_action)) %>% count(moa_action, sort = TRUE) else tibble(moa_action = character(), n = integer())
  
  site_dist <- if (!is.null(bin) && "site_annotation" %in% names(bin)) {
    bin %>%
      mutate(site = case_when(
        str_detect(tolower(site_annotation %||% ""), "allosteric") ~ "Allosteric",
        str_detect(tolower(site_annotation %||% ""), "orthosteric|competitive") ~ "Orthosteric",
        TRUE ~ "Unknown"
      )) %>%
      count(site, sort = TRUE)
  } else tibble(site = character(), n = integer())
  
  omit_reasons <- if (!is.null(omit) && "omission_reason" %in% names(omit)) omit %>% mutate(reason = ifelse(is.na(omission_reason) | !nzchar(omission_reason), "unspecified", omission_reason)) %>% count(reason, sort = TRUE) else tibble(reason = character(), n = integer())
  
  top_targets <- if (!is.null(uniq)) {
    uniq %>% arrange(best_nM) %>% select(uniprot, target_name, target_class, best_nM, best_measure) %>% head(10)
  } else if (!is.null(targ)) {
    targ %>%
      filter(!is.na(uniprot) & nzchar(uniprot)) %>%
      group_by(uniprot) %>%
      summarise(target_name = first(na.omit(pref_name)), best_nM = suppressWarnings(min(best_nM, na.rm = TRUE)), best_measure = first(best_type[which.min(best_nM)]), .groups="drop") %>%
      arrange(best_nM) %>% head(10) %>%
      mutate(target_class = NA_character_) %>% select(uniprot, target_name, target_class, best_nM, best_measure)
  } else tibble()
  
  top_compounds <- if (!is.null(targ)) {
    targ %>%
      filter(!is.na(uniprot) & nzchar(uniprot)) %>%
      count(compound_key, compound_name, sort = TRUE) %>%
      rename(n_targets = n) %>% head(10)
  } else tibble()
  
  platforms <- c("COCONUT (source plant compounds)", "ChEMBL (activities/targets)", "UniProt (target classification)")
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  
  slug <- tolower(gsub("[^a-z0-9]+", "_", plant_name))
  if (is.null(out_md))   out_md   <- sprintf("metadata_%s_%s.md", slug, format(Sys.Date(), "%Y%m%d"))
  if (is.null(out_json)) out_json <- sprintf("metadata_%s_%s.json", slug, format(Sys.Date(), "%Y%m%d"))
  
  cat(
    "# Metadata\n",
    sprintf("**Plant:** %s\n\n", plant_name),
    sprintf("**Timestamp:** %s\n\n", timestamp),
    "## Sources\n",
    paste0("- ", platforms, collapse = "\n"), "\n\n",
    "## Totals\n",
    sprintf("- Total compounds in COCONUT: **%s**\n", ifelse(is.na(total_compounds), "NA", total_compounds)),
    sprintf("- Compounds with ≥1 human protein target: **%s**\n", ifelse(is.na(included_compounds), "NA", included_compounds)),
    sprintf("- Omitted compounds: **%s**\n", ifelse(is.na(omitted_compounds), "NA", omitted_compounds)),
    sprintf("- Unique protein targets (UniProt): **%s**\n\n", ifelse(is.na(unique_targets), "NA", unique_targets)),
    "## Target classes\n",
    if (nrow(class_dist)) paste0(paste0("- ", class_dist$target_class, ": ", class_dist$n), collapse = "\n") else "_NA_", "\n\n",
    "## Potency bins\n",
    if (nrow(potency_dist)) paste0(paste0("- ", potency_dist$potency_bin, ": ", potency_dist$n), collapse = "\n") else "_NA_", "\n\n",
    "## Assay types\n",
    if (nrow(assay_dist)) paste0(paste0("- ", assay_dist$assay_bin, ": ", assay_dist$n), collapse = "\n") else "_NA_", "\n\n",
    "## Mechanism of action (where available)\n",
    if (nrow(moa_dist)) paste0(paste0("- ", moa_dist$moa_action, ": ", moa_dist$n), collapse = "\n") else "_NA_", "\n\n",
    "## Binding site annotation\n",
    if (nrow(site_dist)) paste0(paste0("- ", site_dist$site, ": ", site_dist$n), collapse = "\n") else "_NA_", "\n\n",
    "## Omission reasons\n",
    if (nrow(omit_reasons)) paste0(paste0("- ", omit_reasons$reason, ": ", omit_reasons$n), collapse = "\n") else "_NA_", "\n\n",
    "## Top 10 targets by potency (best_nM)\n",
    if (nrow(top_targets)) paste0(
      paste0("- ", top_targets$target_name, " (", top_targets$uniprot, ", ", coalesce(top_targets$target_class, "NA"),
             ") — best_nM=", signif(top_targets$best_nM, 4), " [", coalesce(top_targets$best_measure, "NA"), "]"),
      collapse = "\n") else "_NA_", "\n\n",
    "## Top 10 compounds by number of targets\n",
    if (nrow(top_compounds)) paste0(
      paste0("- ", coalesce(top_compounds$compound_name, top_compounds$compound_key), " — ", top_compounds$n_targets, " targets"),
      collapse = "\n") else "_NA_", "\n",
    file = out_md
  )
  
  json <- list(
    plant_name = plant_name,
    timestamp = timestamp,
    files = list(
      coconut_csv = coconut_csv,
      products_csv = targets_csv,
      products_binned_csv  = binned_csv,
      targets_csv = unique_binned_csv,
      products_omitted_csv = omitted_csv
    ),
    totals = list(
      total_compounds = total_compounds,
      included_compounds = included_compounds,
      omitted_compounds = omitted_compounds,
      unique_targets = unique_targets
    ),
    distributions = list(
      target_class = if (nrow(class_dist)) setNames(as.list(class_dist$n), class_dist$target_class) else list(),
      potency_bin  = if (nrow(potency_dist)) setNames(as.list(potency_dist$n), potency_dist$potency_bin) else list(),
      assay_bin    = if (nrow(assay_dist)) setNames(as.list(assay_dist$n), assay_dist$assay_bin) else list(),
      moa_action   = if (nrow(moa_dist)) setNames(as.list(moa_dist$n), moa_dist$moa_action) else list(),
      site_annotation = if (nrow(site_dist)) setNames(as.list(site_dist$n), site_dist$site) else list(),
      omission_reasons = if (nrow(omit_reasons)) setNames(as.list(omit_reasons$n), omit_reasons$reason) else list()
    ),
    platforms = platforms
  )
  writeLines(jsonlite::toJSON(json, auto_unbox = TRUE, pretty = TRUE), out_json)
  
  message("Wrote metadata: ", out_md, " and ", out_json)
  invisible(list(markdown = out_md, json = out_json))
}

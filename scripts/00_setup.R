\
    # scripts/00_setup.R
    # Project bootstrap: directories, packages, config, options

    message(">>> Step 0: Project setup starting...")

    # ---------- 1) Robust project root handling ----------
    .pkgs_min <- c("pacman","here","yaml")
    to_install <- setdiff(.pkgs_min, rownames(installed.packages()))
    if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

    suppressPackageStartupMessages({
      library(pacman)
      pacman::p_load(here, yaml)
    })

    # Confirm project root (the folder that contains this script)
    root <- here::here()
    message("Project root: ", root)

    # ---------- 2) Create folders idempotently ----------
    .dirs <- c(
      "data/raw", "data/processed", "data/external",
      "outputs/figures", "outputs/tables", "outputs/models", "outputs/logs",
      "scripts", "config", "reports"
    )
    for (d in .dirs) if (!dir.exists(here::here(d))) dir.create(here::here(d), recursive = TRUE, showWarnings = FALSE)

    # ---------- 3) OS checks: Rtools/Mac tools ----------
    os <- Sys.info()[["sysname"]]
    if (os == "Windows") {
      if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild", repos = "https://cloud.r-project.org")
      has_bt <- tryCatch(pkgbuild::has_build_tools(debug = TRUE), error = function(e) FALSE)
      if (!isTRUE(has_bt)) {
        message("\n*** Windows build tools (Rtools) not detected. Some packages may fail to build.\n",
                "Install Rtools: https://cran.r-project.org/bin/windows/Rtools/\n")
      }
    }

    # ---------- 4) Package set ----------
    pkgs <- c(
      # core genetics & stats
      "adegenet","poppr","pegas","ade4","hierfstat","vegan","ape",
      # spatial & distances
      "spdep","geosphere","sf",
      # modeling & dissimilarity
      "gdm",
      # data handling & viz
      "data.table","dplyr","tidyr","readr","janitor","stringr","purrr",
      "ggplot2","RColorBrewer","reshape2",
      # utilities
      "here","yaml","magrittr"
    )

    message("Installing/loading packages (this may take a bit the first time)...")
    pacman::p_load(char = pkgs, install = TRUE, update = FALSE)

    # ---------- 5) Optional: reproducible env via renv ----------
    # if (interactive() && !requireNamespace("renv", quietly = TRUE)) install.packages("renv")
    # if (interactive()) {
    #   renv::init(bare = TRUE)
    #   # After successful installs, run renv::snapshot()
    # }

    # ---------- 6) Read config ----------
    cfg_path <- here::here("config","config.yml")
    if (!file.exists(cfg_path)) stop("Missing config/config.yml. Please create it (a template is included).")
    cfg <- yaml::read_yaml(cfg_path)

    # Expose cfg to global env for convenience
    assign("cfg", cfg, envir = .GlobalEnv)

    # ---------- 7) Global options & theme ----------
    set.seed(cfg$seed %||% 1234)
    options(stringsAsFactors = FALSE, scipen = 999, repos = c(CRAN = "https://cloud.r-project.org"))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      ggplot2::theme_set(ggplot2::theme_minimal(base_size = 12))
    }

    # ---------- 8) Parallel setup (optional) ----------
    if (isTRUE(cfg$parallel$enable)) {
      n_workers <- cfg$parallel$workers
      if (is.character(n_workers) && n_workers == "auto") {
        n_workers <- max(1L, parallel::detectCores() - 1L)
      }
      options(mc.cores = n_workers)
      message("Parallel enabled with workers: ", n_workers)
    }

    # ---------- 9) Helper functions ----------
    source_if_exists <- function(path) if (file.exists(path)) source(path, local = FALSE)

    source_if_exists(here::here("scripts","utils_io.R"))
    source_if_exists(here::here("scripts","utils_genetics.R"))

    # ---------- 10) Sanity checks for expected input files ----------
    exp_files <- c(
      genotypes = here::here("data","raw","final_individual_genotypes.csv"),
      coords    = here::here(cfg$coords$file),
      malaria   = here::here(cfg$malaria$file)
    )
    exists_vec <- file.exists(unname(exp_files))
    if (!all(exists_vec)) {
      message("\nSome expected files are missing:\n",
              paste(names(exp_files)[!exists_vec], "->", unname(exp_files)[!exists_vec], collapse = "\n"),
              "\n\nAdd these before running analysis scripts.\n",
              "Templates for coords and malaria burden are included in data/raw/.\n")
    } else {
      message("All expected input files found.")
    }

    message(">>> Step 0 complete. Proceed to data preprocessing in scripts/01_preprocess.R.")

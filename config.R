## =============================================================================
## config.R -- paths, gene-set labels and figure output shared by every notebook
##
## Each notebook sources this file first. Point the three directories at your
## copies of the files, either by editing the defaults below or by setting the
## environment variables before starting R or Jupyter:
##
##   RR_DATA_DIR      the files downloaded from EGA and the paper's supplementary
##                    data (default: data/)
##   RR_RESOURCE_DIR  the resources included with the code (default: resources/)
##   RR_FIG_DIR       where figures are written (created if missing)
##
## data/README.md and resources/README.md describe every file.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!exists("REPO_DIR")) REPO_DIR <- Sys.getenv("RR_REPO", unset = ".")
REPO_DIR <- normalizePath(REPO_DIR, mustWork = TRUE)

DATA_DIR     <- Sys.getenv("RR_DATA_DIR",     file.path(REPO_DIR, "data"))
RESOURCE_DIR <- Sys.getenv("RR_RESOURCE_DIR", file.path(REPO_DIR, "resources"))
FIG_DIR      <- Sys.getenv("RR_FIG_DIR",      file.path(REPO_DIR, "figures"))

## ---- Gene sets ---------------------------------------------------------------
## Every variant record and resources/gene_annotation.tsv.gz carry each gene's
## set (Methods): DD-associated, the other genes with pLI >= 0.9, and the rest.
## A gene absent from gene_annotation.tsv.gz is in the pLI < 0.9 set.
GENE_SETS <- c("DD-associated", "pLI \u2265 0.9", "pLI < 0.9")
GS_DD     <- GENE_SETS[1]
GS_HIGH   <- GENE_SETS[2]
GS_LOW    <- GENE_SETS[3]
CONSTRAINT_LEVELS  <- c("Constrained", "Unconstrained")
GENE_SET_COLOURS   <- setNames(c("#1B7837", "#C1852B", "#5B9BD5"), GENE_SETS)
CONSTRAINT_COLOURS <- c("Constrained" = "#1B7837", "Unconstrained" = "#5B9BD5")

## ---- Figures -----------------------------------------------------------------
## save_figure() looks a figure's label and file name up in figures.tsv by key,
## and writes it as PNG (300 dpi) and SVG to FIG_DIR. Inside Jupyter the PNG is
## also displayed in the notebook.
FIG_DPI <- 300
FIGURES <- fread(file.path(REPO_DIR, "figures.tsv"))

fig_item <- function(fig_key) {
  ## base indexing: inside `[.data.table` a bare `key` would be the column
  it <- FIGURES[which(FIGURES[["key"]] == fig_key)]
  if (nrow(it) != 1L) stop("No figures.tsv entry for '", fig_key, "'.", call. = FALSE)
  it
}

#' "Figure 2", "Fig. S4", "Figure N3": for messages
fig_label <- function(key) fig_item(key)$label

save_figure <- function(key, plot, width, height, devices = c("png", "svg")) {
  if (is.null(plot)) {
    message("!! ", key, ": no plot; nothing written")
    return(invisible(NULL))
  }
  it <- fig_item(key)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(FIG_DIR, paste0(it$file, ".", devices))
  for (i in seq_along(devices)) {
    ok <- tryCatch({
      suppressWarnings(ggsave(files[i], plot, device = devices[i], dpi = FIG_DPI,
                              width = width, height = height))
      TRUE
    }, error = function(e) {
      message("!! ", basename(files[i]), " failed: ", conditionMessage(e))
      FALSE
    })
    if (ok) message("  wrote ", it$label, "  ->  ", basename(files[i]))
  }
  png <- files[devices == "png"]
  if (isTRUE(getOption("jupyter.in_kernel")) && length(png) && file.exists(png) &&
      requireNamespace("IRdisplay", quietly = TRUE))
    IRdisplay::display_png(file = png)
  invisible(plot)
}

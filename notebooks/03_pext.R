# %% [markdown]
# # 03 · Inherited PTVs are characterised by partial, rather than complete, transcriptional loss from the affected allele
#
# Transmission of unrescued inherited protein-truncating variants (PTVs) as a
# function of PEXT, the proportion of a gene's expression from the transcripts
# a variant affects. PTVs are placed in the per-transcript NMD regions of Blakes
# et al.: distal PTVs, which escape NMD, are removed, and every other variant's
# PEXT is multiplied by the NMD efficiency of its region (NMD-scaled PEXT).
# Every analysis below uses NMD-scaled PEXT, except the rescue regression
# (Fig. S13A), which uses PEXT.
#
# Figures drawn:
#
# - **Figure 3**: A odds of transmission per unit NMD-scaled PEXT; B observed /
#   expected share of PTVs by NMD-scaled PEXT bin; C continuous dose-response;
#   D within-carrier *de novo* burden
# - **Fig. S11**: PTV origin (A) and transmission (B) by NMD region
# - **Fig. S12**: variance in transmission explained by NMD efficiency, PEXT and
#   NMD-scaled PEXT
# - **Fig. S13**: odds of rescue (A), maternal origin (B) and a diagnosed
#   proband (C) per unit PEXT
# - **Fig. S14**: expected distribution across NMD-scaled PEXT bins
# - **Fig. S15**: transmission excess by NMD-scaled PEXT bin
# - **Fig. S16**: Figure 3D by gene set, with three index rules
# - **Fig. S18**: PEXT of clinically interpreted PTVs against a permuted null
# - **Figure N3**: NMD-scaled PEXT against gene-level constraint (Supplementary Note 3)

# %%
REPO_DIR <- Sys.getenv("RR_REPO", unset = if (file.exists("config.R")) "." else "..")
source(file.path(REPO_DIR, "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(splines)
})
## MASS::mvrnorm draws the coefficient vectors behind the dose-response
## intervals; without it the draws are independent normals
HAVE_MASS <- requireNamespace("MASS", quietly = TRUE)

# %% [markdown]
# ## Inputs
#
# The released files (`data/README.md`), in `DATA_DIR`:
#
# - `parental_variants.tsv.gz`: rare PTVs and synonymous variants in the
#   parents, one row per variant, transcript and parent, with whether each was
#   transmitted, its parent of origin, VEP / LOFTEE annotation, the RescueRanger
#   prediction, the NMD region and NMD-scaled PEXT on that transcript, the gene
#   set, the clinical interpretation of transmitted PTVs in DD-associated genes
#   and the proband's diagnostic status.
# - `proband_variants.tsv.gz`: the probands' calls, for the *de novo* PTVs and
#   their clinical interpretation.
# - `within_carrier_summary.tsv`: one row per carrier and gene set, for the
#   within-carrier fits (Figure 3D, Fig. S16).
#
# and in `RESOURCE_DIR` (`resources/README.md`), the gene annotation (sHet, for
# Figure N3) and the Roulette mutation rate of every possible PTV by gene and
# NMD-scaled PEXT (the *de novo* baseline). Call quality filters, the
# heterozygous-in-proband check on transmitted variants, the removal of
# hypermutated trios, NMD region placement, NMD-scaled PEXT, gene sets and the
# clinical annotation were applied before release.
#
# `RR_BOOT` sets the number of bootstrap replicates for Fig. S12 (default 500).

# %%
PATHS <- list(
  parents        = file.path(DATA_DIR, "parental_variants.tsv.gz"),
  probands       = file.path(DATA_DIR, "proband_variants.tsv.gz"),
  within_carrier = file.path(DATA_DIR, "within_carrier_summary.tsv"),
  genes          = file.path(RESOURCE_DIR, "gene_annotation.tsv.gz"),
  roulette       = file.path(RESOURCE_DIR, "roulette_ptv_grid.tsv.gz"))

# %% [markdown]
# ## Parameters
#
# NMD efficiency per region comes from the allele-specific expression of rare
# nonsense variants in Blakes et al., relative to NMD-target regions (long exon
# about 0.52, start-proximal about 0.29; Methods). A region with half the NMD
# efficiency gets half the PEXT score. Synonymous variants are scaled by the
# region they sit in, so they remain the matched neutral baseline. The released
# tables carry NMD-scaled PEXT; dividing by the region's efficiency gives PEXT.

# %%
PARAMS <- list(
  ## NMD efficiency (Methods): mean alternate allele fractions of rare nonsense
  ## variants in blood RNA-seq by region (Blakes et al. 2026, Fig. 2), with
  ## synonymous variants as the no-NMD reference
  ase_mean_af = c(Synonymous = 0.487, `NMD-target` = 0.383, `Start-proximal` = 0.461,
                  `Long exon` = 0.438, Distal = 0.486),
  ## regions removed outright: distal PTVs escape NMD completely and leave a
  ## near-full-length protein at normal levels, so they can act through altered
  ## function rather than reduced dosage
  scaled_drop = "Distal",
  ## Fig. S11: the regions shown
  region_levels = c("Start-proximal", "Long exon", "Distal", "NMD-target"),
  ## the rescue regression (Fig. S13A) uses raw PEXT and every placed region;
  ## distal unrescued PTVs, which escape NMD, count with the rescued PTVs
  rescue_distal_as_rescued = TRUE,
  ## Fig. S12: bootstrap replicates over genes, and the seed
  nmd_pext_boot = as.integer(Sys.getenv("RR_BOOT", "500")),
  nmd_pext_seed = 20261005L,
  ## within-carrier de novo burden (Figure 3D, Fig. S16)
  within_index_rules = c("sum", "max"),
  ## Fig. S16's open circles: the total (sum) rule refitted without carriers
  ## whose HIGHEST transmitted PEXT is below this (the Low bin's edge)
  within_floor = 0.1,
  ## Figure 3D: the strata drawn. The pooled stratum is the union of
  ## DD-associated and pLI >= 0.9 genes, i.e. every gene with pLI >= 0.9, so it
  ## takes the pLI >= 0.9 colour.
  within_main_sets = c("Constrained (pooled)", GS_LOW),
  ## Fig. S16: the three gene sets in one panel
  within_supp_sets = GENE_SETS,
  within_pooled_colour = unname(GENE_SET_COLOURS[[GS_HIGH]]),
  ## Low [0,0.1), Intermediate [0.1,0.9), High [0.9,1]
  pext_breaks = c(0, 0.1, 0.9, 1),
  pext_levels = c("Low", "Intermediate", "High"),
  constrained_sets = c(GS_DD, GS_HIGH),
  constrained_label = "Constrained (pooled)",
  ptv_class = "Unrescued PTVs", syn_class = "Synonymous", any_ptv_class = "PTVs",
  spline_df = 3, n_sim = 2000L, grid_n = 101L,
  n_perm = 10000L, seed = 1,
  ## Supplementary Note 3: PEXT floor before the log, and the beta grid of the
  ## exchange-rate profile (coarse, then refined near the minimum)
  pext_floor = 1e-3, exchange_grid_coarse = seq(0, 4, by = 0.25),
  exchange_refine_step = 0.05, fig_n3_size = c(width = 14, height = 10),
  ## Figure 3: two rows of two, 17 in wide like Figure 2
  ##   A transmission  | B enrichment
  ##   C dose-response | D within-carrier
  fig3_size = c(width = 17, height = 11.5),
  ## widths of the plotting areas within each row. Bottom row: D's outcome
  ## labels are wrapped onto two lines so C gets the room.
  fig3_top_widths = c(1, 1.3),
  fig3_bottom_widths = c(1, 1.45),
  fig3_row_heights = c(1, 1),
  ## approximate height of each plotting area in inches, used only to give the
  ## error bars Figure 2's bar-end length (cap_in); update if the figure sizes
  ## change
  cap_in = 0.06,
  panel_in = c(fig3 = 3.9, f48 = 3.6, within = 8.9),
  ## and the width of each plotting area (one facet), for vertical error bars
  panel_w_in = c(enrich = 7.2, nmd_tx = 4.4, nmd_r2 = 3.5),
  supp_regression_size = c(width = 18, height = 5.6),
  supp_expected_size = c(width = 10, height = 5.2),
  supp_excess_size = c(width = 11, height = 10.5),
  supp_within_size = c(width = 12, height = 11),
  fig_perm_size = c(width = 12, height = 5.6),
  nmd_supp_fig_size = c(width = 16, height = 10.5),
  nmd_pext_fig_size = c(width = 13, height = 5.4),
  dose_response_point_dodge = 0.022, dose_response_point_size = 3.6)
set.seed(PARAMS$seed)

#' NMD efficiency of each region relative to NMD-target: allele fraction ->
#' alt:ref expression ratio relative to synonymous variants; efficiency = 1 -
#' that ratio; scaled so NMD-target = 1.
nmd_efficiency <- function(af) {
  odds <- function(f) f / (1 - f)
  eff <- 1 - odds(af[setdiff(names(af), "Synonymous")]) / odds(af[["Synonymous"]])
  eff / eff[["NMD-target"]]
}
NMD_ETA  <- nmd_efficiency(PARAMS$ase_mean_af)
NMD_KEPT <- setdiff(names(NMD_ETA), PARAMS$scaled_drop)
## what the analyses call PEXT, for axis titles
PEXT_LAB <- "NMD-scaled PEXT"
## the rescue regression is on raw PEXT
RESCUE_PEXT_LAB <- "PEXT"

# %% [markdown]
# ## Theme and shared helpers
#
# The figures use Figure 2's scale (base font 15, tags 20, axis lines 1.2,
# points 5, error bars 0.8 with about 0.06 in bar-ends), so the main figures and
# the supplement read as one set at the same print width.

# %%
STYLE <- list(base_size = 15, annot_size = 4.2, annot_small = 3.8, tag_size = 20,
              star_size = 5.5,
              point_size = 5, point_size_dense = 3.6, line_width = 0.8,
              axis_width = 1.2, grid_colour = "grey85", grid_width = 0.5,
              axis_colour = "black", text_colour = "black")
PALETTES <- list(
  series = c("All PTVs" = "black", "Unrescued PTVs" = "tomato4",
             "Synonymous" = "grey50"),
  gene_set = GENE_SET_COLOURS,
  pext_bin = setNames(c("#EDF8B1", "#7FCDBB", "#2C7FB8"), PARAMS$pext_levels),
  pext_bin_text = setNames(c("grey15", "grey15", "white"), PARAMS$pext_levels),
  permutation = c(`FALSE` = "grey80", `TRUE` = "tomato4"),
  origin = c("De novo" = "grey30", "Inherited" = "grey80"),
  nmd_model = c(`NMD efficiency` = "tomato4", PEXT = "grey55",
                `NMD-scaled PEXT` = "#1B5E7A"))
SERIES_ORDER <- names(PALETTES$series)
FIG3_TAGS <- c(transmission = "A", enrichment = "B", dose_response = "C",
               within = "D")

msg <- function(...) message(sprintf(...))

## Region names as printed: "NMD target" without the hyphen, as in the other
## figures. The released tables keep "NMD-target"; only labels change.
region_display <- function(x) sub("^NMD-target$", "NMD target", as.character(x))

theme_pub <- function(base_size = STYLE$base_size, grid = c("x", "y", "none")) {
  grid <- match.arg(grid)
  th <- theme_minimal(base_size = base_size) +
    theme(text = element_text(family = "sans", colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
          axis.line.x = element_line(colour = STYLE$axis_colour, linewidth = STYLE$axis_width),
          axis.line.y = element_line(colour = STYLE$axis_colour, linewidth = STYLE$axis_width),
          axis.ticks = element_line(colour = STYLE$axis_colour, linewidth = 0.6),
          axis.ticks.length = unit(3, "pt"),
          axis.text = element_text(size = base_size, colour = STYLE$text_colour),
          axis.text.y = element_text(size = base_size, hjust = 1,
                                     margin = margin(r = 2)),
          axis.title = element_text(size = base_size),
          axis.title.x = element_text(size = base_size, margin = margin(t = 6)),
          strip.text = element_text(face = "bold", size = base_size),
          plot.title = element_text(face = "bold", size = base_size, hjust = 0.5,
                                    margin = margin(b = 6)),
          plot.tag = element_text(face = "bold", size = STYLE$tag_size),
          plot.tag.position = c(0.002, 1),
          legend.position = "bottom",
          legend.title = element_text(size = base_size - 1),
          legend.text = element_text(size = base_size - 1),
          legend.key = element_blank(),
          plot.margin = margin(6, 20, 6, 6))
  gl <- element_line(colour = STYLE$grid_colour, linewidth = STYLE$grid_width)
  if (grid == "x") th <- th + theme(panel.grid.major.x = gl)
  if (grid == "y") th <- th + theme(panel.grid.major.y = gl)
  th
}

#' Error-bar width, in data units, that draws Figure 2's bar-ends
#' (PARAMS$cap_in long) on an axis spanning `span` data units over `panel_in`
#' inches. Under position_dodge() ggplot divides the width by the number of
#' dodged groups, so pass that as `n_dodge`.
cap_width <- function(span, panel_in, n_dodge = 1L, cap_in = PARAMS$cap_in)
  cap_in / panel_in * span * n_dodge
## Forces one text size across a whole figure, applied with `&` at assembly, so
## a panel's own theme() call cannot make it drift.
harmonise_text <- function(base_size = STYLE$base_size)
  theme(text = element_text(size = base_size), axis.text = element_text(size = base_size),
        axis.title = element_text(size = base_size),
        legend.text = element_text(size = base_size - 1),
        legend.title = element_text(size = base_size - 1),
        ## size only: each panel keeps its own face (Figure 3B's bins are plain)
        strip.text = element_text(size = base_size),
        plot.tag = element_text(size = STYLE$tag_size, face = "bold"))

assign_pext_bin <- function(x, breaks = PARAMS$pext_breaks, levels = PARAMS$pext_levels)
  as.character(cut(x, breaks = breaks, labels = levels, include.lowest = TRUE,
                   right = FALSE))
fmt_pct <- function(p) ifelse(abs(p) >= 1, sprintf("%.1f%%", p), sprintf("%.2f%%", p))

## "Nx" labels with a multiplication sign, marked as UTF-8: sprintf() and
## paste() return strings of unknown encoding, which a non-UTF-8 locale draws
## as two stray characters
times <- function(x) { x <- paste0(x, "\u00d7"); Encoding(x) <- "UTF-8"; x }
fold_label <- function(x) times(format(signif(x, 3), drop0trailing = TRUE, trim = TRUE))

## log2 axis labelled as multiples; at most about six breaks over the data
fold_axis <- function(lo, hi) {
  rng <- log2(range(c(lo, hi, 1), na.rm = TRUE, finite = TRUE))
  lo2 <- floor(rng[1]); hi2 <- ceiling(rng[2])
  by <- max(1, ceiling((hi2 - lo2) / 6))
  br <- 2^sort(unique(c(seq(0, lo2, by = -by), seq(0, hi2, by = by))))
  scale_y_continuous(trans = "log2", breaks = br, labels = fold_label,
                     expand = expansion(mult = c(0.05, 0.12)))
}

## Strings carrying "de novo" as plotmath, so it is set in italics without an
## optional package (axis and legend labels accept expressions)
dn_expr <- function(x) {
  out <- lapply(as.character(x), function(s) {
    m <- regmatches(s, regexec("^(.*?)([Dd]e novo)(.*)$", s))[[1]]
    if (!length(m)) return(s)
    q <- function(t) gsub('"', '\\\\"', t)
    parse(text = sprintf('"%s"*italic("%s")*"%s"', q(m[2]), m[3], q(m[4])))[[1]]
  })
  do.call(expression, out)
}

## PEXT bins with their ranges, for axes where the bin is the category
BIN_AXIS_LABELS <- setNames(c("Low\n(< 0.1)", "Intermediate\n(0.1\u20130.9)",
                              "High\n(\u2265 0.9)"), PARAMS$pext_levels)

# %% [markdown]
# ## The released variant tables
#
# Both tables have one row per transcript, and the parental table one per
# parent. `obs_id` labels the rows of one variant in one family (its transcript
# rows and, for parents, the mother's and the father's records): the analyses
# collapse rows to one record per variant with it, and find variants carried by
# both parents. It links nothing else, so no record can be tied to another
# variant of the same child.

# %%
## read as UTF-8, so the gene set names match GENE_SETS in any locale
read_release <- function(f)
  fread(f, colClasses = list(character = "CHROM"), encoding = "UTF-8")
parents_all  <- read_release(PATHS$parents)
probands_all <- read_release(PATHS$probands)

# %% [markdown]
# ## One record per variant
#
# Each transcript row carries the NMD region of the variant on that transcript
# and its NMD-scaled PEXT. Two ways of reducing a variant's rows to one record:
#
# - `collapse_scaled()` (the main analyses): rows in distal regions, and rows
#   that could not be placed (`non_coding`, `unplaced`), are removed. The
#   variant takes the region of its remaining rows with the highest NMD
#   efficiency (`nmd_eta`), and `pext` is its PEXT scaled by that efficiency,
#   the analogue of keeping a variant that is NMD-target on any transcript.
#   `raw` is each row's PEXT, so rows that disagree on it stay apart.
# - `collapse_region()` (the rescue regression, Figs. S11 and S12): one row per
#   variant, NMD-target before long exon before distal before start-proximal
#   before rows that could not be placed (whose region is NA); `pext` is raw
#   PEXT.

# %%
## a row's NMD efficiency; NA for rows that could not be placed
row_eta <- function(region) unname(NMD_ETA[region])

collapse_scaled <- function(d, keep = NMD_KEPT, by = "obs_id") {
  d <- copy(d[nmd_region %in% keep])
  d[, .eta := row_eta(nmd_region)]
  d[, `:=`(.region = nmd_region[which.max(.eta)], nmd_eta = max(.eta)), by = by]
  d[, raw := signif(nmd_scaled_pext / .eta, 10)]
  ## the variant's NMD-scaled PEXT: the value on its row in the chosen region
  d[, pext := if (any(.eta == nmd_eta)) nmd_scaled_pext[.eta == nmd_eta][1L]
              else nmd_scaled_pext[1L] / .eta[1L] * nmd_eta[1L], by = c(by, "raw")]
  d[, nmd_region := .region][, c(".eta", ".region") := NULL][]
}

REGION_ORDER <- c("NMD-target", "Long exon", "Distal", "Start-proximal")
collapse_region <- function(d, by) {
  d <- copy(d)
  d[, .rk := match(nmd_region, REGION_ORDER, nomatch = length(REGION_ORDER) + 1L)]
  setorderv(d, c(by, ".rk"))
  d <- unique(d, by = by)
  d[!nmd_region %in% REGION_ORDER, nmd_region := NA_character_]
  d[, pext := nmd_scaled_pext / row_eta(nmd_region)]
  d[, .rk := NULL][]
}

# %% [markdown]
# ## Parental variants: the transmission table
#
# Three classes of parental variant are kept: unrescued PTVs (RescueRanger
# `not_rescued`, LOFTEE HC), rescued PTVs and synonymous variants. Variants
# carried by both parents are removed because their parent of origin is
# undefined; so is a variant whose kept transcript rows disagree on class,
# consequence, rescue call or PEXT.
#
# `nmd_mode = "filter"` (the main table): `collapse_scaled()` over the
# unrescued PTVs and synonymous variants, and over the rescued PTVs keeping
# distal rows (no figure uses the rescued PTVs of this table; they are carried
# along). `nmd_mode = "annotate"` (the rescue regression, Figs. S11 and S12):
# every PTV by consequence and every synonymous variant, one row per variant
# and class by `collapse_region()`, with `unrescued_hc` (not rescued, LOFTEE
# HC) and `rescued`.

# %%
## columns of the variant records; the child's clinical status and the
## clinical interpretation are the same on every row of a variant
TX_COLS <- c("obs_id", "CHROM", "POS", "REF", "ALT", "study", "vepSYMBOL",
             "Relationship", "Transmitted", "class", "gene_set", "rescue",
             "vepConsequence", "clinical_interpretation", "proband_diagnosis_status",
             "inherited_ptv_diagnosis")

load_parents <- function(nmd_mode = c("filter", "annotate")) {
  nmd_mode <- match.arg(nmd_mode)
  p <- parents_all
  if (nmd_mode == "filter") {
    classes <- rbindlist(list(
      p[is_ptv == TRUE & rescue %in% "rescued"][, class := "Rescued PTVs"],
      p[is_ptv == TRUE & rescue %in% "not_rescued" & vepLoF %in% "HC"][
        , class := "Unrescued PTVs"],
      p[is_ptv == FALSE][, class := "Synonymous"]))
    classes <- rbind(collapse_scaled(classes[class != "Rescued PTVs"]),
                     collapse_scaled(classes[class == "Rescued PTVs"],
                                     keep = names(NMD_ETA)))
    keep <- c(TX_COLS, "raw", "pext", "nmd_region", "nmd_eta")
  } else {
    classes <- rbindlist(list(
      ## every PTV by consequence, no LOFTEE or rescue filter; unrescued_hc marks
      ## the unrescued PTVs of the main table
      p[is_ptv == TRUE][, `:=`(class = PARAMS$any_ptv_class,
                               unrescued_hc = rescue %in% "not_rescued" & vepLoF %in% "HC",
                               rescued = rescue %in% "rescued")],
      p[is_ptv == FALSE][, `:=`(class = PARAMS$syn_class, unrescued_hc = FALSE,
                                rescued = FALSE)]))
    classes <- collapse_region(classes, c("obs_id", "class"))
    keep <- c(TX_COLS, "pext", "nmd_region", "unrescued_hc", "rescued")
  }
  unique(classes[, ..keep])
}

build_overtransmission_table <- function(nmd_mode = "filter") {
  tbl <- load_parents(nmd_mode)
  ## variants present in both parents have undefined parent of origin
  tbl[, parent_count := .N, by = obs_id]
  tbl <- tbl[parent_count == 1][, parent_count := NULL]
  tbl[, gene_set := relevel(factor(gene_set, levels = GENE_SETS), ref = GS_LOW)]
  tbl[, label := relevel(factor(paste0(gene_set, "\n", class)),
                         ref = paste0(GS_LOW, "\nSynonymous"))]
  tbl[]
}

# %% [markdown]
# ## *De novo* PTVs
#
# *De novo* PTVs in protein-coding transcripts with depth > 10 in both parents
# (`parental_depth_ok`) and no RescueRanger UTR flag. In the main set
# (`collapse_scaled()`) distal PTVs are removed and a PEXT value is required,
# as for the inherited variants; rows of one variant on different genes stay
# apart.

# %%
load_observed_dnms <- function(nmd_mode = c("filter", "annotate")) {
  nmd_mode <- match.arg(nmd_mode)
  dnm <- probands_all[is_denovo == TRUE & is_ptv == TRUE & vepBIOTYPE == "protein_coding" &
                        parental_depth_ok == TRUE & utr_flag %in% FALSE]
  cols <- c("obs_id", "CHROM", "POS", "REF", "ALT", "study", "vepSYMBOL", "gene_set",
            "clinical_interpretation")
  out <- if (nmd_mode == "filter") {
    x <- unique(collapse_scaled(dnm)[, c(cols, "raw", "pext", "nmd_region", "nmd_eta"),
                                     with = FALSE])
    x[!is.na(pext)]
  } else collapse_region(dnm, "obs_id")[, c(cols, "pext", "nmd_region"), with = FALSE]
  out[]
}

# %% [markdown]
# ## Roulette expected mutational target
#
# The *de novo* baseline: Roulette mutation rates summed over all possible
# single-nucleotide PTVs, by gene and NMD-scaled PEXT (in percent), built with
# the observed variants' rules (`resources/README.md`): distal transcript rows
# dropped, and each site's PEXT the maximum over transcripts after scaling by
# the NMD efficiency of its region. Genes take the gene set of the observed
# variants; genes with no observed variant are left out.

# %%
get_roulette_grid <- function(gene_set_map) {
  g <- fread(PATHS$roulette)
  g[, PEXT := PEXT_pct / 100]
  g <- merge(g, gene_set_map, by = "GENE")   # inner join: unmapped genes dropped
  stopifnot(all(g$PEXT >= 0 & g$PEXT <= 1), all(is.finite(g$MR)))
  g[]
}
build_gene_set_map <- function(...) {
  parts <- lapply(list(...), function(d)
    as.data.table(d)[, .(GENE = vepSYMBOL, gene_set = as.character(gene_set))])
  unique(rbindlist(parts))[!is.na(GENE) & !is.na(gene_set)][
    , .(gene_set = gene_set[1L]), by = GENE]
}

# %% [markdown]
# ## Per-unit PEXT regressions (Figure 3A, Fig. S13)
#
# Logistic regressions with an interaction between PEXT, gene set and variant
# class (`label`: gene set × unrescued PTV or synonymous). The odds ratio per
# unit PEXT within each class is a linear contrast on the model's coefficients;
# for transmission it is averaged over parent of origin. P values are
# Bonferroni-corrected within each outcome (the five or six strata of one
# outcome are the family), and the asterisks follow the corrected value.
#
# - Transmission (Figure 3A): `Transmitted ~ pext * label * Relationship` over
#   untransmitted and transmitted parental variants.
# - Maternal origin (Fig. S13B): `Relationship ~ pext * label`; the response is
#   ordered Mother < Father, so the sign is flipped to report maternal origin.
# - Rescue (Fig. S13A): from the annotated table at raw PEXT, rescued PTVs in
#   any region plus unrescued LOFTEE HC PTVs in distal regions (which escape
#   NMD) against unrescued LOFTEE HC PTVs in every other region, by gene set.
#   Unplaced variants and splice consequences are excluded.
# - Diagnosis (Fig. S13C): `diagnosed ~ pext * label` over transmitted
#   variants, where a proband is diagnosed with a confident diagnosis or a VUS
#   (`proband_diagnosis_status`). Probands with no clinical record are left
#   out, and so are probands with an inherited PTV diagnosis
#   (`inherited_ptv_diagnosis`: a transmitted unrescued PTV in a DD-associated
#   gene interpreted pathogenic, likely pathogenic or uncertain).

# %%
fit_models <- function(tbl) {
  transmission_dat <- tbl[class != "Rescued PTVs"]
  list(transmission = glm(Transmitted ~ pext * label * Relationship,
                          binomial("logit"), transmission_dat),
       relationship = glm(factor(Relationship, levels = c("Mother", "Father")) ~
                            pext * label, binomial("logit"), transmission_dat))
}

fit_rescue_model <- function(tbl_all, distal_as_rescued = PARAMS$rescue_distal_as_rescued) {
  d <- tbl_all[class == PARAMS$any_ptv_class & (rescued | unrescued_hc) &
                 !grepl("splice", vepConsequence)]
  d <- copy(d[!is.na(nmd_region) & !is.na(pext)])
  d[, rescued_or_escape := as.integer(rescued |
                                        (distal_as_rescued & unrescued_hc &
                                           nmd_region == "Distal"))]
  glm(rescued_or_escape ~ pext * gene_set, binomial("logit"), d)
}

#' Odds ratio for a one-unit increase in PEXT within each level of `class_var`,
#' as a linear contrast L %*% beta on the model's design matrix.
pext_or <- function(model, class_var, classes, avg_over = NULL, flip = FALSE,
                    pext = "pext") {
  tt <- delete.response(terms(model)); vars <- all.vars(tt)
  xl <- model$xlevels; beta <- coef(model); V <- vcov(model)
  classes <- classes[classes %in% xl[[class_var]]]
  contrast_row <- function(cl, avg_val = NULL) {
    a <- list()
    for (v in vars) {
      if (v == pext) next
      else if (v == class_var) a[[v]] <- cl
      else if (!is.null(avg_over) && v == avg_over) a[[v]] <- avg_val
      else if (v %in% names(xl)) a[[v]] <- xl[[v]][1]
      else a[[v]] <- 0
    }
    mm <- function(nd) model.matrix(tt, nd, xlev = xl, contrasts.arg = model$contrasts)
    mm(data.frame(c(a, setNames(list(1), pext)), stringsAsFactors = FALSE)) -
      mm(data.frame(c(a, setNames(list(0), pext)), stringsAsFactors = FALSE))
  }
  rbindlist(lapply(classes, function(cl) {
    avg <- if (is.null(avg_over)) list(NULL) else as.list(xl[[avg_over]])
    L <- Reduce(`+`, lapply(avg, function(v) contrast_row(cl, v))) / length(avg)
    est <- as.numeric(L %*% beta); se <- sqrt(as.numeric(L %*% V %*% t(L)))
    if (flip) est <- -est
    data.table(class = cl, or = exp(est), lo = exp(est - 1.96 * se),
               hi = exp(est + 1.96 * se), p = 2 * pnorm(-abs(est / se)))
  }))
}

#' Odds that a transmitted variant landed in a diagnosed proband (a confident
#' diagnosis or a VUS), per unit PEXT. Probands with an inherited PTV diagnosis
#' are dropped: asking whether high-PEXT inherited PTVs land in diagnosed
#' probands is circular when the diagnosis is an inherited PTV. Probands with no
#' clinical record are dropped.
fit_diagnosis_model <- function(tbl) {
  d <- copy(tbl[class != "Rescued PTVs" & Transmitted == TRUE &
                  proband_diagnosis_status != "not_available"])
  if (!nrow(d)) return(NULL)
  d[, diagnosed := proband_diagnosis_status %in% c("diagnosed", "VUS")]
  d <- d[inherited_ptv_diagnosis == FALSE]
  if (uniqueN(d$diagnosed) < 2L) return(NULL)
  glm(diagnosed ~ pext * label, binomial("logit"), d)
}

regression_config <- function(models) {
  cfg <- list(
    Rescue = list(model = models$rescue, class_var = "gene_set",
                  classes = models$rescue$xlevels$gene_set, avg_over = NULL,
                  flip = FALSE),
    Transmission = list(model = models$transmission, class_var = "label",
                        classes = models$transmission$xlevels$label,
                        avg_over = if ("Relationship" %in%
                                       names(models$transmission$xlevels))
                          "Relationship" else NULL, flip = FALSE),
    Parent = list(model = models$relationship, class_var = "label",
                  classes = models$relationship$xlevels$label, avg_over = NULL,
                  flip = TRUE))
  if (!is.null(models$diagnosis))
    cfg$Diagnosis <- list(model = models$diagnosis, class_var = "label",
                          classes = models$diagnosis$xlevels$label,
                          avg_over = NULL, flip = FALSE)
  cfg
}

collect_estimates <- function(cfg) {
  d <- rbindlist(lapply(names(cfg), function(nm) {
    x <- cfg[[nm]]
    pext_or(x$model, x$class_var, x$classes, x$avg_over, x$flip)[, outcome := nm][]
  }))
  ## Model levels are either "<gene set>" or "<gene set>\n<class>". Matching is
  ## pattern-based so whitespace or encoding differences cannot drop rows.
  cls <- gsub("\\\\n", "\n", as.character(d$class))
  has_nl <- grepl("\n", cls)
  gene_raw <- trimws(ifelse(has_nl, sub("\n.*$", "", cls), cls))
  cons_raw <- trimws(ifelse(has_nl, sub("^.*\n", "", cls), ""))
  d[, gene_set := factor(
    fifelse(grepl("DD|develop", gene_raw, ignore.case = TRUE), GS_DD,
            fifelse(grepl("\u2265|>=", gene_raw), GS_HIGH,
                    fifelse(grepl("<", gene_raw), GS_LOW, NA_character_))),
    levels = rev(GENE_SETS))]                    # y = 1..3, DD-associated on top
  d[, gene_num := as.integer(gene_set)]
  d[, consequence := factor(fifelse(grepl("syn", cons_raw, ignore.case = TRUE),
                                    "Synonymous", "Unrescued PTVs"),
                            levels = c("Unrescued PTVs", "Synonymous"))]
  ## the rescue model (Fig. S13A) models the odds of RESCUE and so covers all PTVs
  d[, series := factor(fifelse(outcome == "Rescue", "All PTVs",
                               as.character(consequence)), levels = SERIES_ORDER)]
  d[, offset := fifelse(outcome == "Rescue", 0,
                        fifelse(consequence == "Unrescued PTVs", 0.18, -0.18))]
  d[, ypos := gene_num + offset]
  ## Bonferroni within each outcome: the five or six strata of one outcome are
  ## the family, and the four outcomes are separate questions.
  d[, p_adj := pmin(1, p * .N), by = outcome]
  d[, stars := fcase(p_adj < 1e-3, "***", p_adj < 1e-2, "**", p_adj < 5e-2, "*",
                     default = "")]
  d[, val := sprintf("%.2f [%.2f\u2013%.2f]", or, lo, hi)]
  d[!is.na(gene_set)]
}

#' Breaks for a log axis spanning lo..hi: round values first (1, 0.5, 2, ...),
#' each kept only if it is at least `min_frac` of the axis away from every
#' break already kept, so labels never run into each other.
log_breaks_within <- function(lo, hi, min_frac = 0.1) {
  priority <- c(1, 0.5, 2, 0.25, 4, 0.1, 10, 1.5, 0.7, 3, 5, 0.2, 20, 1.2, 0.8, 7, 0.9,
                15, 30, 50, 100)
  cand <- priority[priority >= lo * 0.999 & priority <= hi * 1.001]
  gap <- min_frac * log10(hi / lo)
  out <- numeric(0)
  for (x in cand) if (!length(out) || min(abs(log10(out) - log10(x))) >= gap) out <- c(out, x)
  sort(out)
}
#' Log-axis limits covering lo..hi and 1, with room on the right for the
#' estimate labels (`label_room` of the data's log span)
log_limits <- function(lo, hi, label_room = 0.25, pad = 0.06) {
  lo <- min(lo, 1); hi <- max(hi, 1)
  lr <- log(hi / lo)
  c(lo * exp(-pad * lr), hi * exp(max(label_room * lr, log(1.12))))
}

#' Ball-and-stick forest panel: one row per gene set (DD-associated on top),
#' series offset within it, odds ratio on a log axis, error bars with Figure 2's
#' bar-ends. With show_values (Fig. S13) each estimate and its Bonferroni
#' asterisks are printed in a column beside the intervals; without (Figure 3A)
#' only the asterisks, just beyond each interval at Figure 2's size.
#' @param series the legend's entries. Every panel of one figure passes the
#'   same set, so patchwork collects the panels' legends into one.
#' @param panel_in height of the plotting area in inches (bar-end length only)
forest_panel <- function(d, xlab, tag = NULL, xlim = NULL, xbreaks = NULL,
                         series = SERIES_ORDER, panel_in = PARAMS$panel_in[["f48"]],
                         show_y = TRUE, label_room = 0.25, show_values = TRUE) {
  d <- copy(d)
  d[, lab := paste0(sprintf("%.2f", or), fifelse(stars == "", "", paste0(" ", stars)))]
  if (!show_values) label_room <- min(label_room, 0.1)
  if (is.null(xlim)) xlim <- log_limits(min(d$lo), max(d$hi), label_room)
  ylim <- c(0.5, 3.5); yadd <- 0.2
  cap <- cap_width(diff(ylim) + 2 * yadd, panel_in)
  ## one invisible layer carrying every series, so each panel emits an identical
  ## legend and patchwork can collect them into a single guide
  legend_df <- data.frame(series = factor(series, levels = series))
  ggplot(d) +
    geom_vline(xintercept = 1, colour = "black", linewidth = 0.6) +
    geom_errorbar(aes(y = ypos, xmin = lo, xmax = hi, colour = series), orientation = "y",
                  width = cap, linewidth = STYLE$line_width, show.legend = FALSE) +
    geom_point(aes(x = or, y = ypos, colour = series), size = STYLE$point_size,
               show.legend = FALSE) +
    ## one aligned column right of the widest interval, so a label can never
    ## sit on the line at 1 or on another interval
    {if (show_values)
      geom_text(aes(x = max(hi) * 1.08, y = ypos, label = lab, colour = series), hjust = 0,
                vjust = 0.45, size = STYLE$annot_small, fontface = "bold",
                show.legend = FALSE)
     else
      geom_text(data = d[stars != ""], aes(x = hi * 1.03, y = ypos, label = stars,
                                            colour = series),
                hjust = 0, vjust = 0.75, size = STYLE$star_size, show.legend = FALSE)} +
    geom_point(data = legend_df, aes(x = 1, y = 2, colour = series), size = 0,
               na.rm = TRUE, show.legend = TRUE) +
    scale_colour_manual(values = PALETTES$series, limits = series, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 4.5))) +
    scale_x_log10(breaks = if (is.null(xbreaks)) log_breaks_within(xlim[1], xlim[2]) else xbreaks,
                  labels = function(x) formatC(x, format = "fg"), limits = xlim) +
    scale_y_continuous(breaks = 1:3, labels = if (show_y) rev(GENE_SETS) else NULL,
                       limits = ylim, expand = expansion(add = yadd)) +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = NULL, tag = tag) +
    theme_pub(grid = "x")
}

# %% [markdown]
# ## Observed versus expected by NMD-scaled PEXT bin (Figure 3B, Fig. S14)
#
# PTVs are binned at NMD-scaled PEXT 0.1 and 0.9 (Low, Intermediate, High). For
# *de novo* PTVs, the share in each bin is compared with the share of the
# Roulette mutational target (exact binomial test and interval). For
# transmitted inherited unrescued PTVs, the share is compared with the share of
# transmitted synonymous variants in the same gene set (risk ratio with a 0.5
# Haldane-Anscombe correction and log risk-ratio intervals). Each cell is tested
# against its own null, Bonferroni-corrected across all gene set × bin ×
# variant type cells together.

# %%
dnm_enrichment <- function(observed, expected) {
  res <- merge(CJ(gene_set = unique(expected$gene_set), BIN = PARAMS$pext_levels,
                  unique = TRUE), expected[, .(gene_set, BIN, exp_frac)],
               by = c("gene_set", "BIN"), all.x = TRUE)
  res <- merge(res, observed, by = c("gene_set", "BIN"), all.x = TRUE)
  res[is.na(obs_n), obs_n := 0L]
  res[, N_total := sum(obs_n), by = gene_set]
  res <- res[!is.na(exp_frac) & exp_frac > 0 & N_total > 0]
  bt <- Map(function(o, n, p) binom.test(as.integer(o), as.integer(n), p),
            res$obs_n, res$N_total, res$exp_frac)
  res[, obs_prop := obs_n / N_total]
  res[, lo := vapply(bt, function(b) b$conf.int[1], numeric(1))]
  res[, hi := vapply(bt, function(b) b$conf.int[2], numeric(1))]
  res[, p_value := vapply(bt, function(b) b$p.value, numeric(1))]
  res[, `:=`(log2_fold = log2(obs_prop / exp_frac), log2_lo = log2(lo / exp_frac),
             log2_hi = log2(hi / exp_frac))]
  res[, variant_type := "De novo"]
  res[, .(gene_set, BIN, variant_type, baseline = "Roulette mutational target",
          n_observed = obs_n, n_total = N_total, observed_share = obs_prop,
          expected_share = exp_frac, n_expected = exp_frac * N_total,
          fold = obs_prop / exp_frac, fold_lo = lo / exp_frac, fold_hi = hi / exp_frac,
          log2_fold, log2_lo, log2_hi, p_value)]
}

## Inherited PTVs against the synonymous baseline within each gene set
## (Haldane-Anscombe 0.5 correction, log risk-ratio CIs).
inherited_enrichment <- function(ptv, syn, conf = 0.95) {
  cells <- CJ(gene_set = unique(c(ptv$gene_set, syn$gene_set)),
              BIN = PARAMS$pext_levels, unique = TRUE)
  r <- merge(merge(cells, ptv, by = c("gene_set", "BIN"), all.x = TRUE), syn,
             by = c("gene_set", "BIN"), all.x = TRUE)
  r[is.na(a), a := 0][is.na(c), c := 0]
  r[, n1 := sum(a), by = gene_set][, n2 := sum(c), by = gene_set]
  z <- qnorm(1 - (1 - conf) / 2)
  r[, a2 := a + 0.5][, c2 := c + 0.5]
  r[, rr := (a2 / n1) / (c2 / n2)]
  r[, se := sqrt(1 / a2 - 1 / n1 + 1 / c2 - 1 / n2)]
  r[!is.finite(se) | se <= 0, se := sqrt(1 / a2 + 1 / c2)]
  r[, `:=`(log2_fold = log2(rr), log2_lo = log2(rr) - z * se / log(2),
           log2_hi = log2(rr) + z * se / log(2))]
  r[, p_value := 2 * pnorm(-abs(log(rr) / se))]
  r[, variant_type := "Inherited PTV"]
  r[, .(gene_set, BIN, variant_type, baseline = "Transmitted synonymous",
        n_observed = a, n_total = n1, n_baseline = c, n_baseline_total = n2,
        observed_share = a / n1, expected_share = c / n2, n_expected = (c / n2) * n1,
        fold = rr, fold_lo = 2^log2_lo, fold_hi = 2^log2_hi, log2_fold, log2_lo,
        log2_hi, p_value)]
}

expected_target <- function(roulette_grid) {
  e <- roulette_grid[, .(MR = sum(MR)), by = .(gene_set, BIN = assign_pext_bin(PEXT))]
  e <- e[!is.na(BIN)]
  e[, exp_frac := MR / sum(MR), by = gene_set]
  e[]
}

pext_bin_enrichment <- function(tbl, dnms_observed, roulette_grid) {
  expected <- expected_target(roulette_grid)
  obs <- copy(dnms_observed)[!is.na(pext)]
  obs[, gene_set := as.character(gene_set)][, BIN := assign_pext_bin(pext)]
  inherited <- copy(tbl[Transmitted == TRUE & !is.na(pext)])
  inherited[, gene_set := as.character(gene_set)][, BIN := assign_pext_bin(pext)]
  ## Transmitted synonymous variants ARE the baseline the inherited enrichment
  ## is measured against, so this draws that denominator explicitly (Fig. S14A).
  expected_inh <- inherited[class == PARAMS$syn_class & !is.na(BIN),
                            .(n_syn = .N), by = .(gene_set, BIN)]
  expected_inh[, exp_frac := n_syn / sum(n_syn), by = gene_set]
  both <- rbind(
    dnm_enrichment(obs[!is.na(BIN), .(obs_n = .N), by = .(gene_set, BIN)], expected),
    inherited_enrichment(inherited[class == PARAMS$ptv_class & !is.na(BIN),
                                   .(a = .N), by = .(gene_set, BIN)],
                         inherited[class == PARAMS$syn_class & !is.na(BIN),
                                   .(c = .N), by = .(gene_set, BIN)]),
    use.names = TRUE, fill = TRUE)
  ## each row tests that cell against its own null, Bonferroni-corrected across
  ## all gene set x bin x variant type cells together
  both[, p_adj := p.adjust(p_value, "bonferroni")]
  both[, sig := fcase(p_adj < 1e-3, "***", p_adj < 1e-2, "**", p_adj < 5e-2, "*",
                      default = "ns")]
  list(expected = expected, expected_inherited = expected_inh, both = both)
}

#' Figure 3B: one facet per PEXT bin, labelled along the bottom like a
#' categorical axis; within each, the gene sets side by side, de novo (circle)
#' and inherited (square) PTVs next to each other. Colour is the gene set
#' throughout. The axis is log2 but labelled as multiples (1x = as expected).
#' Error bars with Figure 2's bar-ends; asterisks: Bonferroni (ns not shown).
panel_enrichment_grouped <- function(both, tag, panel_w_in = PARAMS$panel_w_in[["enrich"]]) {
  d <- copy(both)[, `:=`(group = as.character(gene_set),
                         series = fifelse(variant_type == "De novo", "dn", "inh"))]
  groups <- c(GS_DD, GS_HIGH, GS_LOW)
  d[, `:=`(colour_set = as.character(gene_set),
           series = factor(series, levels = c("dn", "inh")),
           BIN = factor(as.character(BIN), levels = PARAMS$pext_levels,
                        labels = unname(BIN_AXIS_LABELS[PARAMS$pext_levels])))]
  d[, x := match(group, groups) + fifelse(series == "dn", -0.18, 0.18)]
  d[, `:=`(star_y = fifelse(log2_fold >= 0, log2_hi + 0.06, log2_lo - 0.06),
           star_vjust = fifelse(log2_fold >= 0, 0, 1))]
  xlim <- c(0.5, length(groups) + 0.5)
  cap <- cap_width(diff(xlim) * 1.1, panel_w_in / length(PARAMS$pext_levels))
  yr <- range(c(d$log2_lo, d$log2_hi, 0), finite = TRUE)
  ybr <- seq(floor(yr[1]), ceiling(yr[2]), by = max(1, ceiling(diff(yr) / 6)))
  labs_series <- c(dn = "Observed de novo PTVs (versus Roulette expected)",
                   inh = "Observed inherited PTVs (vs rare transmitted synonymous)")
  present <- levels(droplevels(d$series))
  ggplot(d, aes(x, log2_fold)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.6) +
    geom_errorbar(aes(ymin = log2_lo, ymax = log2_hi), width = cap, colour = "black",
                  linewidth = 0.6) +
    geom_point(aes(shape = series, fill = colour_set), colour = "black", size = 4,
               stroke = 0.7) +
    geom_text(data = d[!sig %in% c("ns", "")],
              aes(y = star_y, label = sig, colour = colour_set, vjust = star_vjust),
              size = STYLE$annot_small, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ BIN, nrow = 1, strip.position = "bottom") +
    scale_shape_manual(values = c(dn = 21, inh = 22)[present],
                       labels = dn_expr(labs_series[present]), name = NULL,
                       guide = guide_legend(order = 2, ncol = 1, override.aes = list(
                         fill = "grey45", colour = "black", size = 4.5))) +
    scale_fill_manual(values = PALETTES$gene_set, breaks = GENE_SETS, name = NULL,
                      guide = guide_legend(order = 1, override.aes = list(
                        shape = 21, colour = "black", size = 4.5))) +
    scale_colour_manual(values = PALETTES$gene_set, guide = "none") +
    scale_x_continuous(breaks = NULL, limits = xlim) +
    scale_y_continuous(breaks = ybr, labels = function(b) fold_label(2^b),
                       expand = expansion(mult = 0.08)) +
    labs(x = paste(PEXT_LAB, "bin"), y = "Observed / expected\nshare of PTVs", tag = tag) +
    theme_pub(grid = "y") +
    theme(strip.placement = "outside", strip.text = element_text(face = "plain",
                                                                 lineheight = 0.9),
          panel.spacing.x = unit(1.2, "lines"), legend.box = "vertical",
          legend.spacing.y = unit(0.05, "cm"), legend.text = element_text(hjust = 0))
}

#' Fig. S14: the share of each gene set's expected PTVs falling in each PEXT bin
panel_expected <- function(expected, tag, xlab) {
  d <- copy(expected)
  d[, gene_set := factor(gene_set, levels = rev(GENE_SETS))]
  d[, BIN := factor(BIN, levels = PARAMS$pext_levels)]
  d[, lab := fifelse(exp_frac >= 0.06, paste0(round(100 * exp_frac), "%"), "")]
  d[, lab_colour := PALETTES$pext_bin_text[as.character(BIN)]]
  ggplot(d, aes(exp_frac, gene_set, fill = BIN)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3,
             position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = lab, colour = lab_colour),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = STYLE$annot_small) +
    scale_colour_identity() +
    scale_fill_manual(values = PALETTES$pext_bin, name = paste(PEXT_LAB, "bin"),
                      guide = guide_legend(title.position = "left")) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1),
                       expand = expansion(c(0, 0.01))) +
    labs(x = xlab, y = NULL, tag = tag) + theme_pub(grid = "none") +
    theme(legend.title = element_text(size = STYLE$base_size - 1))
}

# %% [markdown]
# ## Continuous dose-response (Figure 3C)
#
# One model per gene set, `Transmitted ~ pext * class + Relationship`, over
# unrescued PTVs and synonymous variants. The curve is the fitted excess
# transmission rate (PTV minus synonymous), averaged over parent of origin, with
# a 95% band from 2,000 draws of the coefficient vector. The points are the
# observed excess in each PEXT bin (at the bin's mean PEXT, 95% CI; bins with at
# least 30 PTVs and 30 synonymous variants).

# %%
design_matrix <- function(model, newdata) {
  tt <- delete.response(terms(model))
  ## going through model.frame() is what picks up predvars, i.e. the ns() knots;
  ## building from terms() directly would refit the basis on the grid
  model.matrix(tt, model.frame(tt, newdata, xlev = model$xlevels),
               contrasts.arg = model$contrasts)
}
predict_draws <- function(model, newdata, B = PARAMS$n_sim, seed = PARAMS$seed) {
  X <- design_matrix(model, newdata); beta <- coef(model); ok <- !is.na(beta)
  X <- X[, names(beta)[ok], drop = FALSE]; V <- vcov(model)
  set.seed(seed)
  draws <- if (HAVE_MASS) MASS::mvrnorm(B, beta[ok], V)
           else matrix(rnorm(B * sum(ok), rep(beta[ok], each = B),
                             rep(sqrt(diag(V)), each = B)), nrow = B)
  plogis(draws %*% t(X))
}

observed_bins <- function(d, min_n = 30L) {
  x <- as.data.table(d)[!is.na(px) & !is.na(Transmitted)]
  brk <- PARAMS$pext_breaks
  brk[1] <- min(brk[1], min(x$px, na.rm = TRUE))
  brk[length(brk)] <- max(brk[length(brk)], max(x$px, na.rm = TRUE))
  ## right = FALSE matches assign_pext_bin(), so a variant at exactly 0.1 or 0.9
  ## falls in the same bin as everywhere else in the paper
  x[, bin := cut(px, breaks = brk, include.lowest = TRUE, right = FALSE, labels = FALSE)]
  x <- x[!is.na(bin)][, is_ptv := as.character(class) == PARAMS$ptv_class]
  agg <- x[, .(mean_pext = mean(px[is_ptv]), N_ptv = sum(is_ptv),
               T_ptv = sum(is_ptv & Transmitted), N_syn = sum(!is_ptv),
               T_syn = sum(!is_ptv & Transmitted)), by = bin]
  agg <- agg[N_ptv >= min_n & N_syn >= min_n & is.finite(mean_pext)]
  if (!nrow(agg)) return(NULL)
  agg[, bin_label := PARAMS$pext_levels[bin]]
  agg[, `:=`(p_ptv = T_ptv / N_ptv, p_syn = T_syn / N_syn)]
  agg[, e := p_ptv - p_syn]
  agg[, se_e := sqrt(p_ptv * (1 - p_ptv) / N_ptv + p_syn * (1 - p_syn) / N_syn)]
  agg[, `:=`(excess_pp = 100 * e, lo_pp = 100 * (e - 1.96 * se_e),
             hi_pp = 100 * (e + 1.96 * se_e))]
  setorder(agg, bin)
  agg[]
}

fit_dose_response <- function(tbl, gs, B = PARAMS$n_sim) {
  d <- as.data.table(tbl)[gene_set == gs &
                            as.character(class) %in% c(PARAMS$syn_class,
                                                       PARAMS$ptv_class)]
  d[, px := pext]
  d <- d[!is.na(px) & !is.na(Transmitted)]
  d[, class := factor(as.character(class), levels = c(PARAMS$syn_class,
                                                      PARAMS$ptv_class))]
  d[, Transmitted := as.logical(Transmitted)]
  n_ptv <- sum(d$class == PARAMS$ptv_class)
  if (n_ptv < 50 || uniqueN(d$class) < 2L) { msg("!! %s: too few PTVs; skipped", gs)
    return(NULL) }
  has_rel <- "Relationship" %in% names(d) && uniqueN(na.omit(d$Relationship)) > 1L
  if (has_rel) d[, Relationship := droplevels(factor(Relationship))]
  rhs_rel <- if (has_rel) " + Relationship" else ""

  m_lin <- glm(as.formula(sprintf("Transmitted ~ px * class%s", rhs_rel)),
               binomial("logit"), data = d)

  grid <- seq(0, 1, length.out = PARAMS$grid_n)
  rel_levels <- if (has_rel) levels(d$Relationship) else NA_character_
  n_rel <- if (has_rel) length(rel_levels) else 1L
  make_nd <- function(cl) {
    nd <- expand.grid(px = grid, Relationship = rel_levels,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    nd$class <- factor(cl, levels = levels(d$class))
    if (has_rel) nd$Relationship <- factor(nd$Relationship, levels = rel_levels)
    else nd$Relationship <- NULL
    nd
  }
  nd_ptv <- make_nd(PARAMS$ptv_class); nd_syn <- make_nd(PARAMS$syn_class)
  n1 <- nrow(nd_ptv)
  ## expand.grid varies its FIRST argument fastest, so the columns run px 1..n
  ## in one block per parent; each block is weighted 1/n_parents
  avg_parents <- function(M) {
    out <- matrix(0, nrow(M), length(grid))
    for (r in seq_len(n_rel))
      out <- out + (1 / n_rel) * M[, (r - 1L) * length(grid) + seq_along(grid), drop = FALSE]
    out
  }
  P <- predict_draws(m_lin, rbind(nd_ptv, nd_syn), B = B)
  E_lin <- 100 * (avg_parents(P[, seq_len(n1), drop = FALSE]) -
                    avg_parents(P[, n1 + seq_len(nrow(nd_syn)), drop = FALSE]))
  curve <- data.table(gene_set = gs, pext = grid,
                      excess_pp = apply(E_lin, 2, median),
                      lo_pp = apply(E_lin, 2, quantile, 0.025),
                      hi_pp = apply(E_lin, 2, quantile, 0.975))
  obs <- observed_bins(d)
  if (!is.null(obs)) obs[, gene_set := gs]
  list(curve = curve, observed = obs)
}

#' Figure 3C: the linear fit with its 95% band, and the observed excess in each
#' PEXT bin. The binned intervals carry no bar-ends, as in Figure 2E. The dotted
#' lines are the bin edges of Figure 3B.
panel_dose_response <- function(fits, tag) {
  curves <- rbindlist(lapply(fits, `[[`, "curve"))
  bins   <- rbindlist(Filter(Negate(is.null), lapply(fits, `[[`, "observed")),
                      use.names = TRUE)
  curves[, gene_set := factor(gene_set, levels = GENE_SETS)]
  p <- ggplot(mapping = aes(colour = gene_set, fill = gene_set)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.6) +
    geom_vline(xintercept = PARAMS$pext_breaks[2:3], linetype = "dotted",
               colour = "grey45", linewidth = 0.6) +
    geom_ribbon(data = curves, aes(pext, ymin = lo_pp, ymax = hi_pp), alpha = 0.15,
                colour = NA, show.legend = FALSE) +
    geom_line(data = curves, aes(pext, excess_pp), linewidth = 1.2)
  if (nrow(bins)) {
    bins[, gene_set := factor(gene_set, levels = GENE_SETS)]
    ## fixed offset per gene set rather than jitter: the three sets have almost
    ## the same mean PEXT within a bin and a random offset would move between runs
    dodge <- PARAMS$dose_response_point_dodge
    off <- setNames(seq(-dodge, dodge, length.out = length(GENE_SETS)), GENE_SETS)
    bins[, x_pos := pmin(pmax(mean_pext + off[as.character(gene_set)], 0), 1)]
    p <- p +
      geom_linerange(data = bins, aes(x_pos, ymin = lo_pp, ymax = hi_pp),
                     colour = "black", linewidth = STYLE$line_width, inherit.aes = FALSE) +
      geom_point(data = bins, aes(x_pos, excess_pp, fill = gene_set), shape = 21,
                 colour = "black", size = PARAMS$dose_response_point_size,
                 stroke = 0.7, inherit.aes = FALSE, show.legend = FALSE)
  }
  p + scale_colour_manual(values = PALETTES$gene_set, name = NULL) +
    scale_fill_manual(values = PALETTES$gene_set, guide = "none") +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.6))) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = PEXT_LAB, y = "Excess transmission rate (pp)", tag = tag) +
    theme_pub(grid = "y")
}

# %% [markdown]
# ## Within-carrier *de novo* burden (Figure 3D, Fig. S16)
#
# Among probands carrying at least one transmitted unrescued PTV in a stratum
# (a gene set, or the pooled constrained genes), does lower NMD-scaled PEXT of
# those inherited (index) PTVs go with more damaging *de novo* burden, i.e. did
# carriers of lower-PEXT variants need more *de novo* burden to be ascertained?
# Each carrier's index PTVs are summarised as their highest NMD-scaled PEXT
# (Max), their total (Sum; Figure 3D), or the total refitted without carriers
# whose highest value is below 0.1 (Sum, floor). Outcomes: any *de novo* PTV and
# a *de novo* PTV in a DD-associated gene (logistic), the sHet-weighted *de novo*
# PTV burden and the synonymous *de novo* count as a negative control
# (quasi-Poisson). All models adjust for study and ancestry group, and all but
# the control for log(synonymous *de novo* count + 1). Estimates are per unit
# DECREASE in PEXT, exp(-beta): above 1 means carriers of lower-PEXT variants
# carried more *de novo* burden.
#
# > **Per-individual data:** read from `within_carrier_summary.tsv`: one row per
# > carrier and stratum, with the sum and the maximum NMD-scaled PEXT of the
# > carrier's transmitted unrescued PTVs in the stratum, the number of those
# > PTVs, the same proband's *de novo* outcomes, study and ancestry group.
# > Probands with more than 20 *de novo* calls of any class were left out.

# %%
WITHIN_POOLED <- PARAMS$constrained_label
WITHIN_STRATA <- c(setNames(as.list(GENE_SETS), GENE_SETS),
                   setNames(list(PARAMS$constrained_sets), WITHIN_POOLED))
WITHIN_OUTCOMES <- c("Any de novo PTV", "De novo PTV in a DD-associated gene",
                     "sHet-weighted de novo PTV burden",
                     "Synonymous de novo count (control)")
WITHIN_RULE_CODES <- c(sum = "Sum", max = "Max", floor = "Sum, floor")
## Fig. S16: highest = triangle, total = filled circle, total without the
## carriers whose highest PEXT is below the floor = open circle
WITHIN_FIG_RULES <- c("Max", "Sum", "Sum, floor")
WITHIN_RULE_LABELS <- c(Max = "Highest PEXT", Sum = "Total PEXT",
                        `Sum, floor` = sprintf("Total PEXT, carriers with highest < %.1f removed",
                                               PARAMS$within_floor))
WITHIN_RULE_SHAPES <- c(Max = 17, Sum = 16, `Sum, floor` = 21)
WITHIN_COLOURS <- c(GENE_SET_COLOURS, setNames(PARAMS$within_pooled_colour, WITHIN_POOLED))
## stratum names as printed in Figure 3D (legend)
WITHIN_DISPLAY <- setNames(c("Constrained (pooled pLI \u2265 0.9)",
                             paste0("Unconstrained (", GS_LOW, ")")),
                           c(WITHIN_POOLED, GS_LOW))
within_label <- function(x) { x <- as.character(x)
  ifelse(x %in% names(WITHIN_DISPLAY), WITHIN_DISPLAY[x], x) }
## outcome labels on two lines, for Figure 3D (narrower, so C gets the room)
WITHIN_OUTCOME_LINES <- list(
  "Any de novo PTV" = "Any de novo PTV",
  "De novo PTV in a DD-associated gene" = c("De novo PTV in a", "DD-associated gene"),
  "sHet-weighted de novo PTV burden" = c("sHet-weighted", "de novo PTV burden"),
  "Synonymous de novo count (control)" = c("Synonymous de novo", "count (control)"))
#' plotmath labels with "de novo" in italics; two-line labels are stacked with
#' atop() (plotmath has no line break)
within_outcome_labels <- function(outcomes, wrap = FALSE) {
  one <- function(t) {
    m <- regmatches(t, regexec("^(.*?)([Dd]e novo)(.*)$", t))[[1]]
    if (!length(m)) return(sprintf('"%s"', t))
    sprintf('"%s"*italic("%s")*"%s"', m[2], m[3], m[4])
  }
  do.call(expression, lapply(outcomes, function(o) {
    lines <- if (wrap) WITHIN_OUTCOME_LINES[[o]] else o
    if (is.null(lines)) lines <- o
    txt <- if (length(lines) == 1L) one(lines) else
      sprintf("atop(displaystyle(%s), displaystyle(%s))", one(lines[1]), one(lines[2]))
    parse(text = txt)[[1]]
  }))
}

#' The within-carrier file, as one table per stratum: carriers with at least one
#' transmitted unrescued PTV with a PEXT value in the stratum.
load_within_base <- function() {
  w <- fread(PATHS$within_carrier, encoding = "UTF-8")
  split(w, by = "stratum", keep.by = FALSE)
}

#' Every within-carrier estimate: stratum x index rule x outcome x adjustment
#' for the proband's synonymous de novo count.
fit_within <- function(wb, rules = PARAMS$within_index_rules,
                       floor_pext = PARAMS$within_floor,
                       covars = " + study + ancestry_group") {
  within <- rbindlist(lapply(names(WITHIN_STRATA), function(nm) {
    dc0 <- wb[[nm]]
    if (is.null(dc0) || nrow(dc0) < 50) return(NULL)
    rbindlist(lapply(c(rules, "floor"), function(rule) {
      dc <- if (rule == "floor") copy(dc0[index_pext_max >= floor_pext]) else copy(dc0)
      if (nrow(dc) < 50) return(NULL)
      dc[, index_pext := switch(rule, max = index_pext_max, index_pext_sum)]
      if (!is.finite(sd(dc$index_pext, na.rm = TRUE)) ||
          sd(dc$index_pext, na.rm = TRUE) == 0) return(NULL)
      grab <- function(m, outcome, adj) {
        if (is.null(m)) return(NULL)
        s <- summary(m)$coefficients
        if (!"index_pext" %in% rownames(s)) return(NULL)
        e <- s["index_pext", 1]; se <- s["index_pext", 2]
        data.table(gene_set = nm, outcome = outcome,
                   index_rule = WITHIN_RULE_CODES[[rule]],
                   adjusted_for_synonymous = adj, n = nrow(dc),
                   rr_per_unit_decrease = exp(-e), lo = exp(-e - 1.96 * se),
                   hi = exp(-e + 1.96 * se), p_value = s["index_pext", 4])
      }
      fmla <- function(lhs, adj) as.formula(paste0(lhs, " ~ index_pext", covars,
                                                   if (adj) " + log(dnm_synonymous_n + 1)"
                                                   else ""))
      safe_glm <- function(f, fam) tryCatch(glm(f, fam, data = dc),
                                            error = function(e) NULL)
      out <- list()
      for (adj in c(TRUE, FALSE)) {
        if (uniqueN(dc$dnm_ptv_any) > 1L)
          out <- c(out, list(grab(safe_glm(fmla("dnm_ptv_any", adj), binomial("logit")),
                                  "Any de novo PTV", adj)))
        if (uniqueN(dc$dnm_ptv_dd_any) > 1L)
          out <- c(out, list(grab(safe_glm(fmla("dnm_ptv_dd_any", adj), binomial("logit")),
                                  "De novo PTV in a DD-associated gene", adj)))
      }
      if (sum(dc$dnm_synonymous_n, na.rm = TRUE) > 0 &&
          sd(dc$dnm_synonymous_n, na.rm = TRUE) > 0)
        out <- c(out, list(grab(safe_glm(fmla("dnm_synonymous_n", FALSE),
                                         quasipoisson(link = "log")),
                                "Synonymous de novo count (control)", FALSE)))
      if (sd(dc$dnm_shet_burden, na.rm = TRUE) > 0)
        for (adj in c(TRUE, FALSE))
          out <- c(out, list(grab(safe_glm(fmla("dnm_shet_burden", adj),
                                           quasipoisson(link = "log")),
                                  "sHet-weighted de novo PTV burden", adj)))
      rbindlist(Filter(Negate(is.null), out), use.names = TRUE)
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE)
  if (is.null(within) || !nrow(within)) return(NULL)
  within[, gene_set := factor(gene_set, levels = c(GENE_SETS, WITHIN_POOLED))]
  within[, index_rule := factor(index_rule, levels = unname(WITHIN_RULE_CODES))]
  setorder(within, gene_set, index_rule, outcome, -adjusted_for_synonymous)
  within[]
}

## the estimates drawn: the adjusted fit, and the control (unadjusted by
## construction, since the synonymous count is its own outcome)
within_drawn <- function(est, sets, rules)
  est[as.character(gene_set) %in% sets & as.character(index_rule) %in% rules &
        outcome %in% WITHIN_OUTCOMES &
        (adjusted_for_synonymous == TRUE | outcome == "Synonymous de novo count (control)")]

#' Figure 3D and Fig. S16: one panel, outcomes as rows. Within each outcome
#' the strata are offset and coloured; with more than one index rule (Fig.
#' S16) each stratum's rules are offset again, in a block, and told apart by
#' shape. Ball-and-stick with Figure 2's bar-ends; asterisks at nominal P.
panel_within <- function(est, sets = PARAMS$within_main_sets, rules = "Sum", tag = NULL,
                         xlab = NULL, panel_in = PARAMS$panel_in[["fig3"]],
                         x_range = NULL, relabel = TRUE, wrap = FALSE) {
  d <- copy(within_drawn(est, sets, rules))
  if (!nrow(d)) return(NULL)
  d[, `:=`(stratum = factor(as.character(gene_set), levels = sets),
           rule = factor(as.character(index_rule), levels = rules))]
  d[, y0 := match(outcome, rev(WITHIN_OUTCOMES))]   # the control at the bottom
  multi <- length(rules) > 1L
  ## slot of each estimate within its outcome, first stratum (and rule) on top;
  ## a gap of half a slot between strata when the rules are blocked
  ks <- length(sets); kr <- length(rules)
  gap <- if (multi) 0.5 else 0
  d[, slot := (as.integer(stratum) - 1L) * (kr + gap) + (as.integer(rule) - 1L)]
  n_slot <- (ks - 1L) * (kr + gap) + kr
  step <- min(0.3, 0.84 / max(1, n_slot - 1))
  d[, ypos := y0 + ((n_slot - 1) / 2 - slot) * step]
  d[, stars := fcase(p_value < 1e-3, "***", p_value < 1e-2, "**", p_value < 5e-2, "*",
                     default = "")]
  d[, lab := stars]
  ylim <- c(0.5, length(WITHIN_OUTCOMES) + 0.5); yadd <- 0.1
  cap <- cap_width(diff(ylim) + 2 * yadd, panel_in)
  ## aligned column right of the widest interval
  d[, lab_x := hi * 1.03]
  rng <- if (is.null(x_range)) range(c(d$lo, d$hi), finite = TRUE) else x_range
  xlim <- log_limits(rng[1], rng[2], label_room = 0.1)
  ggplot(d, aes(y = ypos, colour = stratum)) +
    geom_vline(xintercept = 1, colour = "black", linewidth = 0.6) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = cap,
                  linewidth = STYLE$line_width, show.legend = FALSE) +
    geom_point(aes(x = rr_per_unit_decrease, shape = rule),
               size = if (multi) STYLE$point_size_dense else STYLE$point_size,
               fill = "white", stroke = 1.2) +
    geom_text(data = d[lab != ""], aes(x = lab_x, label = lab), hjust = 0, vjust = 0.75,
              size = STYLE$star_size, show.legend = FALSE) +
    scale_colour_manual(values = WITHIN_COLOURS, breaks = sets,
                        labels = if (relabel) within_label else waiver(), name = NULL,
                        guide = guide_legend(order = 1, override.aes = list(shape = 16,
                                                                            size = 4.5))) +
    scale_shape_manual(values = WITHIN_RULE_SHAPES[rules],
                       labels = WITHIN_RULE_LABELS[rules], name = NULL,
                       guide = if (multi) guide_legend(order = 2, override.aes = list(
                         colour = "black", size = 4.5)) else "none") +
    scale_x_log10(breaks = log_breaks_within(xlim[1], xlim[2]),
                  labels = function(x) formatC(x, format = "fg"), limits = xlim) +
    scale_y_continuous(breaks = seq_along(WITHIN_OUTCOMES),
                       labels = within_outcome_labels(rev(WITHIN_OUTCOMES), wrap),
                       limits = ylim, expand = expansion(add = yadd)) +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = NULL, tag = tag) +
    theme_pub(grid = "x") +
    theme(legend.box = "vertical", legend.spacing.y = unit(0.1, "cm"))
}

#' Figure 3: two rows, each its own grid, so a row's panels share a baseline
#' while the long outcome labels of D do not push C away from the left edge.
figure3 <- function(pA, pB, pC, pD) {
  row <- function(widths, ...) wrap_plots(Filter(Negate(is.null), list(...)), nrow = 1,
                                          widths = widths)
  wrap_plots(row(PARAMS$fig3_top_widths, pA, pB), row(PARAMS$fig3_bottom_widths, pC, pD),
             ncol = 1, heights = PARAMS$fig3_row_heights) & harmonise_text()
}

# %% [markdown]
# ## Building the analysis tables
#
# `tbl` is the main transmission table (NMD-scaled PEXT, distal PTVs removed);
# `tbl_regions` labels every variant with its NMD region at raw PEXT, for the
# rescue regression and Figs. S11 and S12.

# %%
tbl <- build_overtransmission_table()
tbl_regions <- build_overtransmission_table(nmd_mode = "annotate")

models <- fit_models(tbl)
models$rescue <- fit_rescue_model(tbl_regions)
models$diagnosis <- fit_diagnosis_model(tbl)
## one set of estimates covering all four outcomes: Figure 3A draws the
## transmission rows and Fig. S13 the other three
estimates <- collect_estimates(regression_config(models))

dnms_observed <- load_observed_dnms()
gene_set_map  <- build_gene_set_map(dnms_observed, tbl[Transmitted == TRUE])
roulette_grid <- get_roulette_grid(gene_set_map)
enrichment    <- pext_bin_enrichment(tbl, dnms_observed, roulette_grid)
## the denominators the Figure 3B legend quotes: n_total repeats on every bin
## row of a gene set, so it is taken once per gene set
den <- enrichment$both[, .(n_total = n_total[1L], n_baseline_total = n_baseline_total[1L]),
                       by = .(variant_type, gene_set)]
msg("Figure 3B denominators: %s transmitted unrescued PTVs, %s de novo PTVs, %s transmitted synonymous variants",
    format(den[variant_type != "De novo", sum(n_total)], big.mark = ","),
    format(den[variant_type == "De novo", sum(n_total)], big.mark = ","),
    format(den[variant_type != "De novo", sum(n_baseline_total)], big.mark = ","))

dose_fits <- Filter(Negate(is.null),
                    setNames(lapply(GENE_SETS, function(gs) fit_dose_response(tbl, gs)),
                             GENE_SETS))

within_base <- load_within_base()
within_est  <- fit_within(within_base)
if (!is.null(within_est)) {
  shown <- within_drawn(within_est, PARAMS$within_main_sets, "Sum")
  for (i in seq_len(nrow(shown)))
    msg("Figure 3D  %-22s %-38s %.2f (%.2f-%.2f), P = %.2g, %s carriers",
        shown$gene_set[i], shown$outcome[i], shown$rr_per_unit_decrease[i], shown$lo[i],
        shown$hi[i], shown$p_value[i], format(shown$n[i], big.mark = ","))
} else msg("!! within-carrier fits: no stratum had 50 carriers; Figure 3 has no panel D")

# %% [markdown]
# ## Figure 3
#
# A: odds of transmission per unit NMD-scaled PEXT for unrescued PTVs and
# synonymous variants in each gene set (only the Bonferroni asterisks are
# printed). B: observed / expected share of *de novo* and inherited PTVs in each
# NMD-scaled PEXT bin. C: fitted excess transmission rate against NMD-scaled
# PEXT. D: within-carrier *de novo* burden per unit decrease in the total
# NMD-scaled PEXT of the carrier's inherited PTVs, constrained genes pooled
# against pLI < 0.9 genes.
#
# > **Per-individual data:** panel D only, from `within_carrier_summary.tsv`
# > (see the within-carrier section above).

# %%
p3A <- forest_panel(estimates[outcome == "Transmission"],
                    paste("Odds of transmission per unit", PEXT_LAB, "(95% CI)"),
                    FIG3_TAGS[["transmission"]],
                    series = c(PARAMS$ptv_class, PARAMS$syn_class),
                    panel_in = PARAMS$panel_in[["fig3"]], show_values = FALSE)
p3B <- panel_enrichment_grouped(enrichment$both, FIG3_TAGS[["enrichment"]])
p3C <- panel_dose_response(dose_fits, FIG3_TAGS[["dose_response"]])
p3D <- if (!is.null(within_est))
  panel_within(within_est, tag = FIG3_TAGS[["within"]], wrap = TRUE,
               xlab = paste0("Per unit decrease in the total ", PEXT_LAB,
                             "\nof the carrier's inherited PTVs (95% CI)"))
save_figure("fig3", figure3(p3A, p3B, p3C, p3D),
            PARAMS$fig3_size[["width"]], PARAMS$fig3_size[["height"]])

# %% [markdown]
# ## Fig. S11: PTV origin and transmission by NMD region
#
# Every variant labelled with its Blakes et al. region, split by gene set. Regions are ordered by the share of *de novo*
# PTVs in DD-associated genes. A: of the PTVs carried by probands in each region
# (every PTV by consequence, no LOFTEE or rescue filter), the share that arose
# *de novo* rather than being transmitted. B: one logistic regression of
# transmission over parental variants, `Transmitted ~ label`, with one label per
# gene set × {synonymous (pooled over regions), unrescued PTVs in each region}
# and synonymous variants in pLI < 0.9 genes as the reference (open circle).
# Wald 95% CIs; Bonferroni across every contrast with the reference.

# %%
## exact binomial 95% interval for k of n
binom_interval <- function(k, n) {
  ci <- Map(function(a, b) if (b > 0) binom.test(a, b)$conf.int else c(NA_real_, NA_real_),
            k, n)
  list(lo = vapply(ci, `[`, numeric(1), 1L), hi = vapply(ci, `[`, numeric(1), 2L))
}

denovo_share_by_region <- function(tbl_all, dnms_all, levels = PARAMS$region_levels) {
  dn <- dnms_all[!is.na(nmd_region), .(n_de_novo = .N),
                 by = .(gene_set = as.character(gene_set), nmd_region)]
  inh <- tbl_all[class == PARAMS$any_ptv_class & as.logical(Transmitted) %in% TRUE &
                   !is.na(nmd_region), .(n_inherited = .N),
                 by = .(gene_set = as.character(gene_set), nmd_region)]
  out <- merge(CJ(gene_set = GENE_SETS, nmd_region = levels), dn,
               by = c("gene_set", "nmd_region"), all.x = TRUE)
  out <- merge(out, inh, by = c("gene_set", "nmd_region"), all.x = TRUE)
  out[is.na(n_de_novo), n_de_novo := 0L][is.na(n_inherited), n_inherited := 0L]
  out[, n_total := n_de_novo + n_inherited]
  ci <- binom_interval(out$n_de_novo, out$n_total)
  out[, `:=`(share_de_novo = fifelse(n_total > 0, n_de_novo / n_total, NA_real_),
             lo = ci$lo, hi = ci$hi)]
  out[, nmd_region := factor(nmd_region, levels = levels)]
  out[]
}

nmd_star <- function(p) fcase(p < 1e-3, "***", p < 1e-2, "**", p < 5e-2, "*", default = "")

transmission_model_by_region <- function(tbl_all, levels = PARAMS$region_levels) {
  d <- tbl_all[!is.na(nmd_region) & !is.na(Transmitted) &
                 ((class == PARAMS$any_ptv_class & unrescued_hc) |
                    class == PARAMS$syn_class)]
  d <- copy(d)[, `:=`(gene_set = as.character(gene_set),
                      y = as.logical(Transmitted),
                      category = fifelse(class == PARAMS$syn_class, PARAMS$syn_class,
                                         as.character(nmd_region)))]
  ref <- paste(GS_LOW, PARAMS$syn_class, sep = " | ")
  d[, label := relevel(factor(paste(gene_set, category, sep = " | ")), ref = ref)]
  fit <- glm(y ~ label, binomial("logit"), d)
  s <- summary(fit)$coefficients
  ci <- suppressMessages(confint.default(fit))
  ## rows matched to labels by position, not by name: in a non-UTF-8 locale
  ## R mangles the ">=" sign of the pLI >= 0.9 label in coefficient names
  stopifnot(nrow(s) == nlevels(d$label), nrow(ci) == nlevels(d$label))
  est <- data.table(label = levels(d$label), or = exp(s[, "Estimate"]),
                    or_lo = exp(ci[, 1]), or_hi = exp(ci[, 2]), p_value = s[, "Pr(>|z|)"])
  ## exp(intercept) is the reference group's transmission odds, not a ratio
  est[label == ref, `:=`( or = 1, or_lo = NA_real_,
                                   or_hi = NA_real_, p_value = NA_real_)]
  out <- merge(d[, .(n = .N, n_transmitted = sum(y)), by = .(label, gene_set, category)],
               est, by = "label", all.x = TRUE)
  out[, is_reference := label == ref]
  out[, p_adj := p.adjust(p_value, "bonferroni")]   # NA reference not counted
  out[, sig := fifelse(is_reference, "", nmd_star(p_adj))]
  out[, `:=`(rate = n_transmitted / n,
             category = factor(category, levels = c(PARAMS$syn_class, levels)),
             variant_class = factor(fifelse(category == PARAMS$syn_class,
                                            PARAMS$syn_class, PARAMS$ptv_class),
                                    levels = c(PARAMS$ptv_class, PARAMS$syn_class)))]
  out[, label := NULL][]
}

#' Order of the regions on both panels and in every facet: highest share of de
#' novo PTVs first, in DD-associated genes, where the share varies (it is near
#' zero everywhere in pLI < 0.9 genes). Ties and regions without PTVs keep the
#' order of PARAMS$region_levels.
nmd_region_order <- function(denovo, by_gene_set = GS_DD,
                             levels = PARAMS$region_levels) {
  d <- denovo[gene_set == by_gene_set]
  share <- d$share_de_novo[match(levels, as.character(d$nmd_region))]
  levels[order(-share, seq_along(levels), na.last = TRUE)]
}

## Panel A keeps the % de novo in each bar and the number of PTVs above it.
panel_nmd_denovo <- function(b, region_order, tag = "A") {
  d <- copy(b)[, `:=`(gene_set = factor(gene_set, levels = GENE_SETS),
                      nmd_region = factor(as.character(nmd_region), levels = region_order))]
  ## shares computed here rather than by position_fill(): scale limits act on
  ## the raw values, so counts would fall outside a 0-1 axis and be dropped
  long <- rbind(d[n_total > 0, .(gene_set, nmd_region, origin = "De novo",
                                 share = n_de_novo / n_total)],
                d[n_total > 0, .(gene_set, nmd_region, origin = "Inherited",
                                 share = n_inherited / n_total)])
  long[, origin := factor(origin, levels = c("De novo", "Inherited"))]
  lab <- d[n_total > 0]
  ggplot(long, aes(nmd_region, share, fill = origin)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7, colour = "white",
             linewidth = 0.3) +
    geom_text(data = lab, aes(nmd_region, pmin(share_de_novo + 0.05, 0.93),
                              label = sprintf("%.0f%%", 100 * share_de_novo)),
              inherit.aes = FALSE, size = STYLE$annot_small, colour = "black") +
    geom_text(data = lab, aes(nmd_region, 1.05,
                              label = format(n_total, big.mark = ",", trim = TRUE)),
              inherit.aes = FALSE, size = STYLE$annot_small - 0.7, colour = "grey30") +
    facet_wrap(~ gene_set, nrow = 1) +
    scale_x_discrete(labels = region_display) +
    scale_fill_manual(values = PALETTES$origin, labels = dn_expr(c("De novo", "Inherited")),
                      name = NULL) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25),
                       limits = c(0, 1.09), expand = expansion(0)) +
    labs(x = NULL, y = "PTVs in probands", tag = tag) +
    theme_pub(grid = "none") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right",
          legend.text = element_text(hjust = 0))
}

## regions in the order of panel A, then synonymous variants
panel_nmd_transmission <- function(tx, region_order, tag = "B",
                                   facet_w_in = PARAMS$panel_w_in[["nmd_tx"]]) {
  d <- copy(tx)[, `:=`(gene_set = factor(gene_set, levels = GENE_SETS),
                       category = factor(as.character(category),
                                         levels = c(region_order, PARAMS$syn_class)))]
  d[!is.na(or_hi) & !is.finite(or_hi), or_hi := NA_real_]
  ref <- d[is_reference == TRUE]
  ## discrete axis: n categories plus ggplot's 0.6 on either side
  cap <- cap_width(nlevels(d$category) + 1.2, facet_w_in)
  ggplot(d, aes(category, or, colour = variant_class)) +
    scale_x_discrete(drop = FALSE, labels = region_display) +
    geom_hline(yintercept = 1, colour = "black", linewidth = 0.6) +
    geom_errorbar(aes(ymin = or_lo, ymax = or_hi), width = cap,
                  linewidth = STYLE$line_width, na.rm = TRUE) +
    geom_point(size = STYLE$point_size_dense) +
    geom_point(data = ref, size = STYLE$point_size_dense, shape = 21, fill = "white",
               stroke = 1.1, show.legend = FALSE) +
    geom_text(aes(y = or_hi, label = sig), vjust = -0.25, size = STYLE$annot_small,
              na.rm = TRUE, show.legend = FALSE) +
    facet_wrap(~ gene_set, nrow = 1, drop = FALSE) +
    fold_axis(d$or_lo, d$or_hi) +
    scale_colour_manual(values = PALETTES$series[c(PARAMS$ptv_class, PARAMS$syn_class)],
                        name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    labs(x = NULL, y = sprintf("Transmission odds relative to\nsynonymous in %s genes", GS_LOW),
         tag = tag) +
    theme_pub(grid = "y") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
}

nmd_supp_figure <- function(denovo, tx, region_order = nmd_region_order(denovo))
  wrap_plots(panel_nmd_denovo(denovo, region_order),
             panel_nmd_transmission(tx, region_order), ncol = 1) & harmonise_text()

# %%
nmd_denovo <- denovo_share_by_region(tbl_regions, load_observed_dnms(nmd_mode = "annotate"))
nmd_tx     <- transmission_model_by_region(tbl_regions)
nmd_order  <- nmd_region_order(nmd_denovo)
save_figure("figS11", nmd_supp_figure(nmd_denovo, nmd_tx, nmd_order),
            PARAMS$nmd_supp_fig_size[["width"]], PARAMS$nmd_supp_fig_size[["height"]])

# %% [markdown]
# ## Fig. S12: variance in transmission explained
#
# Does NMD-scaled PEXT explain more of the variation in the transmission of
# unrescued PTVs than PEXT alone or NMD efficiency alone? Unrescued PTVs (RESCUE
# 0, LOFTEE HC) in NMD-target, long-exon and start-proximal regions (distal and
# unplaced PTVs excluded). One predictor per model, each with parent of origin:
# NMD efficiency, PEXT, or PEXT × NMD efficiency. Variance explained is
# McFadden's pseudo-R² against the parent-only model. 95% intervals come from a
# bootstrap over genes (PTVs in one gene are not independent).

# %%
NMD_R2_MODELS <- c("NMD efficiency", "PEXT", "NMD-scaled PEXT", "PEXT + NMD efficiency")
NMD_R2_PLOTTED <- NMD_R2_MODELS[1:3]

nmd_pext_r2 <- function(tbl_all, n_boot = PARAMS$nmd_pext_boot, seed = PARAMS$nmd_pext_seed) {
  regions <- NMD_KEPT
  ptv <- tbl_all[class == PARAMS$any_ptv_class & unrescued_hc &
                   !is.na(Transmitted) & !is.na(Relationship)]
  ptv <- copy(ptv[nmd_region %in% regions & !is.na(pext)])
  ptv[, `:=`(gene_set = as.character(gene_set), gene = as.character(vepSYMBOL),
             y = as.numeric(as.logical(Transmitted)), eta = unname(NMD_ETA[nmd_region]))]
  ptv[, pext_eta := pext * eta]
  ## glm.fit on prebuilt design matrices: the same fits as glm(), fast enough to
  ## bootstrap. Returns the log-likelihood of each model.
  loglik <- function(X, y) {
    f <- suppressWarnings(glm.fit(X, y, family = binomial()))
    p <- pmin(pmax(f$fitted.values, 1e-12), 1 - 1e-12)
    sum(y * log(p) + (1 - y) * log(1 - p))
  }
  rel_matrix <- function(d) if (uniqueN(d$Relationship) > 1L)
    model.matrix(~ factor(Relationship), d) else matrix(1, nrow(d), 1L)
  fit_all <- function(d) {
    rel <- rel_matrix(d)
    ll0 <- loglik(rel, d$y)
    ll <- c(loglik(cbind(rel, d$eta), d$y), loglik(cbind(rel, d$pext), d$y),
            loglik(cbind(rel, d$pext_eta), d$y), loglik(cbind(rel, d$pext, d$eta), d$y))
    names(ll) <- NMD_R2_MODELS
    list(ll0 = ll0, ll = ll, r2 = 1 - ll / ll0, k0 = ncol(rel))
  }
  set.seed(seed)
  per_set <- lapply(setNames(GENE_SETS, GENE_SETS), function(gs) {
    d <- ptv[gene_set == gs]
    if (nrow(d) < 50L || uniqueN(d$y) < 2L) return(NULL)
    est <- fit_all(d)
    rows <- split(seq_len(nrow(d)), d$gene)
    boot <- t(vapply(seq_len(n_boot), function(b) {
      i <- unlist(rows[sample.int(length(rows), replace = TRUE)], use.names = FALSE)
      fit_all(d[i])$r2
    }, numeric(length(NMD_R2_MODELS))))
    colnames(boot) <- NMD_R2_MODELS
    extra <- c(1, 1, 1, 2)
    chi2 <- 2 * (est$ll - est$ll0)
    fit <- data.table(
      gene_set = gs, model = NMD_R2_MODELS, n_ptvs = nrow(d), n_genes = length(rows),
      parameters = extra, r2_pct = 100 * est$r2,
      r2_lo = 100 * apply(boot, 2, quantile, 0.025, na.rm = TRUE),
      r2_hi = 100 * apply(boot, 2, quantile, 0.975, na.rm = TRUE),
      lr_chi2 = chi2, lr_p = pchisq(chi2, extra, lower.tail = FALSE),
      AIC = -2 * est$ll + 2 * (est$k0 + extra))
    fit[, delta_AIC := AIC - min(AIC)]
    ## direction: R2 counts an association either way, so the odds ratio per
    ## unit of each single predictor (> 1 = more transmission with more
    ## predicted dosage loss)
    rel <- rel_matrix(d)
    slope <- function(x) { b <- suppressWarnings(glm.fit(cbind(rel, x), d$y,
                                                         family = binomial()))$coefficients
      exp(b[length(b)]) }
    fit[, or_per_unit := c(slope(d$eta), slope(d$pext), slope(d$pext_eta), NA_real_)]
    fit
  })
  list(fit = rbindlist(Filter(Negate(is.null), per_set)))
}

## Free y scales: the three gene sets differ by orders of magnitude.
panel_nmd_pext <- function(r2, facet_w_in = PARAMS$panel_w_in[["nmd_r2"]]) {
  d <- r2$fit[model %in% NMD_R2_PLOTTED]
  if (!nrow(d)) return(NULL)
  d[, `:=`(gene_set = factor(gene_set, levels = GENE_SETS),
           model = factor(model, levels = NMD_R2_PLOTTED))]
  cap <- cap_width(length(NMD_R2_PLOTTED) + 1.2, facet_w_in)
  ggplot(d, aes(model, r2_pct, fill = model)) +
    geom_col(width = 0.65) +
    geom_errorbar(aes(ymin = r2_lo, ymax = r2_hi), width = cap,
                  linewidth = STYLE$line_width) +
    facet_wrap(~ gene_set, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = PALETTES$nmd_model, guide = "none") +
    scale_x_discrete(labels = c(`NMD efficiency` = "NMD\nefficiency", PEXT = "PEXT",
                                `NMD-scaled PEXT` = "NMD-scaled\nPEXT")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL,
         y = expression(atop(displaystyle("Variance in transmission"),
                             displaystyle("explained (pseudo-" * R^2 * ", %)")))) +
    theme_pub(grid = "y") +
    theme(panel.spacing.x = unit(2, "lines"), axis.text.x = element_text(lineheight = 0.9))
}

# %%
nmd_r2 <- nmd_pext_r2(tbl_regions)
p_np <- panel_nmd_pext(nmd_r2)
if (!is.null(p_np))
  save_figure("figS12", p_np & harmonise_text(),
              PARAMS$nmd_pext_fig_size[["width"]], PARAMS$nmd_pext_fig_size[["height"]])

# %% [markdown]
# ## Fig. S13: rescue, maternal origin and diagnosis per unit PEXT
#
# The three outcomes of the per-unit PEXT regressions not drawn in Figure 3,
# from the same estimates (see "Per-unit PEXT regressions" above). A: odds of
# rescue per unit raw PEXT, over all PTVs. B: odds of maternal origin per unit
# NMD-scaled PEXT. C: odds that a transmitted variant is carried by a diagnosed
# proband, per unit NMD-scaled PEXT. Each estimate is printed with its
# Bonferroni asterisks. Panel C reads the proband's diagnostic status and
# inherited PTV diagnosis carried on every parental record.

# %%
s13_spec <- list(Rescue = paste0("Odds of rescue\nper unit ", RESCUE_PEXT_LAB, " (95% CI)"),
                 Parent = paste0("Odds of maternal origin\nper unit ", PEXT_LAB, " (95% CI)"),
                 Diagnosis = paste0("Odds of a diagnosed proband\nper unit ", PEXT_LAB,
                                    " (95% CI)"))
s13_spec <- s13_spec[names(s13_spec) %in% unique(estimates$outcome)]
s13_panels <- lapply(names(s13_spec), function(oc)
  forest_panel(estimates[outcome == oc], s13_spec[[oc]], label_room = 0.32))
if (length(s13_panels))
  save_figure("figS13", wrap_plots(s13_panels, nrow = 1, guides = "collect") +
                plot_annotation(tag_levels = "A") & theme(legend.position = "bottom") &
                harmonise_text(),
              PARAMS$supp_regression_size[["width"]], PARAMS$supp_regression_size[["height"]])

# %% [markdown]
# ## Fig. S14: expected distribution across NMD-scaled PEXT bins
#
# How each baseline of Figure 3B is split across the Low, Intermediate and High
# NMD-scaled PEXT bins, by gene set. A: transmitted synonymous variants, the
# baseline for inherited PTVs, each binned by the NMD-scaled PEXT at its
# position. B: the Roulette mutational target, the baseline for *de novo* PTVs.

# %%
save_figure("figS14", wrap_plots(
  panel_expected(enrichment$expected_inherited, "A",
                 "Transmitted synonymous variants (inherited baseline)"),
  ## plotmath rather than markdown, so the italic survives without an optional
  ## package: the label is one expression, not a string with markup in it
  panel_expected(enrichment$expected, "B",
                 expression("Roulette expected mutational target (" * italic("de novo") *
                              " baseline)")),
  ncol = 1, guides = "collect") & theme(legend.position = "bottom") & harmonise_text(),
  PARAMS$supp_expected_size[["width"]], PARAMS$supp_expected_size[["height"]])

# %% [markdown]
# ## Fig. S15: transmission excess by NMD-scaled PEXT bin
#
# Within each gene set and bin, the excess transmission rate is the
# transmission rate of unrescued PTVs minus that of synonymous variants in the
# same bin, and the number of excess transmissions is that rate times the
# number of PTVs. A: total excess transmissions per gene set, stacked by bin;
# bold labels give the total and the excess as a percentage of PTVs. B: each
# bin's share of the excess minus its share of PTVs, in percentage points,
# against the null that every PTV contributes equally.

# %%
excess_by_bin <- function(tbl, ptv_class = PARAMS$ptv_class,
                          syn_class = PARAMS$syn_class) {
  dt <- tbl[class %in% c(ptv_class, syn_class) & !is.na(Transmitted) & !is.na(pext)]
  dt[, gene_set := as.character(gene_set)]
  dt[, BIN := factor(assign_pext_bin(pext), levels = PARAMS$pext_levels)]
  dt <- dt[!is.na(BIN)]
  dt[, `:=`(is_ptv = class == ptv_class, tr = as.logical(Transmitted))]
  agg <- dt[, .(T_ptv = sum(is_ptv & tr), U_ptv = sum(is_ptv & !tr),
                T_syn = sum(!is_ptv & tr), U_syn = sum(!is_ptv & !tr)),
            by = .(gene_set, BIN)]
  agg[, `:=`(N_ptv = T_ptv + U_ptv, N_syn = T_syn + U_syn)]
  agg <- agg[N_ptv > 0 & N_syn > 0]
  agg[, `:=`(p_ptv = T_ptv / N_ptv, p_syn = T_syn / N_syn)]
  agg[, e := p_ptv - p_syn]
  agg[, se_e := sqrt(p_ptv * (1 - p_ptv) / N_ptv + p_syn * (1 - p_syn) / N_syn)]
  agg[, excess := N_ptv * e]
  agg[, `:=`(obs_share = excess / sum(excess), exp_share = N_ptv / sum(N_ptv)),
      by = gene_set]
  agg[, `:=`(delta_pp = 100 * (obs_share - exp_share), ratio = obs_share / exp_share)]
  setorderv(agg, c("gene_set", "BIN"))
  agg[]
}

panel_excess_total <- function(agg, tag = NULL) {
  d <- copy(agg)
  d[, gene_set := factor(gene_set, levels = GENE_SETS)]
  d[, pct := 100 * excess / N_ptv]
  setorder(d, gene_set, BIN)
  d[, sgn := sign(excess)]
  d[, cum_end := cumsum(excess), by = .(gene_set, sgn)]
  d[, cum_start := cum_end - excess][, mid := (cum_start + cum_end) / 2]
  d[, bar_span := sum(abs(excess)), by = gene_set]
  d[, lab := sprintf("%.0f (%s)", excess, fmt_pct(pct))]
  d[, inside := abs(excess) >= 0.10 * bar_span]
  d[, txt_colour := PALETTES$pext_bin_text[as.character(BIN)]]
  net <- d[, .(net = sum(excess), N_ptv = sum(N_ptv), top = sum(pmax(excess, 0))),
           by = gene_set]
  net[, lab := sprintf("%.0f (%s)", net, fmt_pct(100 * net / N_ptv))]
  ymax <- max(net$top, 0); if (ymax <= 0) ymax <- max(abs(d$excess))
  PAD <- 0.13 * ymax
  ymin <- min(0, d$cum_start, d$cum_end) - PAD
  d[, y_below := -PAD / 2]
  ybrk <- scales::breaks_pretty(4)(c(0, ymax)); ybrk <- ybrk[ybrk >= 0]
  ggplot(d, aes(gene_set, excess, fill = BIN)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3,
             position = position_stack(reverse = TRUE)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_text(data = d[inside == TRUE], aes(y = mid, label = lab, colour = txt_colour),
              size = STYLE$annot_small) +
    geom_text(data = d[inside == FALSE], aes(y = y_below, label = lab),
              colour = "grey25", size = STYLE$annot_small) +
    geom_text(data = net, aes(x = gene_set, y = top, label = lab), inherit.aes = FALSE,
              vjust = -0.6, size = STYLE$annot_small, fontface = "bold",
              colour = "grey20") +
    scale_fill_manual(values = PALETTES$pext_bin, name = paste(PEXT_LAB, "bin"),
                      breaks = PARAMS$pext_levels) +
    scale_colour_identity() + scale_x_discrete(limits = GENE_SETS) +
    scale_y_continuous(limits = c(ymin, ymax * 1.16), breaks = ybrk,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Total excess transmissions", tag = tag) + theme_pub(grid = "y")
}

panel_excess_departure <- function(agg, tag) {
  d <- copy(agg)
  d[, gene_set := factor(gene_set, levels = rev(GENE_SETS))]
  setorder(d, gene_set, BIN)
  d[, sgn := sign(delta_pp)]
  d[, cum_end := cumsum(delta_pp), by = .(gene_set, sgn)]
  d[, cum_start := cum_end - delta_pp][, mid := (cum_start + cum_end) / 2]
  d[, lab := sprintf("%+.0f", delta_pp)]
  d[, txt_colour := PALETTES$pext_bin_text[as.character(BIN)]]
  ggplot(d, aes(delta_pp, gene_set, fill = BIN)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3,
             position = position_stack(reverse = TRUE)) +
    geom_text(aes(x = mid, label = lab, colour = txt_colour), size = STYLE$annot_small) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
    scale_fill_manual(values = PALETTES$pext_bin, name = paste(PEXT_LAB, "bin"),
                      breaks = PARAMS$pext_levels) +
    scale_colour_identity() +
    scale_x_continuous(labels = function(x) paste0(x, " pp"),
                       expand = expansion(mult = 0.05)) +
    labs(x = "Departure from an equal contribution to the excess",
         y = NULL, tag = tag) + theme_pub(grid = "none")
}

# %%
excess <- excess_by_bin(tbl)
## A above B, the legend below both
save_figure("figS15", wrap_plots(panel_excess_total(excess, "A"),
                                 panel_excess_departure(excess, "B"),
                                 ncol = 1, heights = c(1.2, 1), guides = "collect") &
              theme(legend.position = "bottom", legend.box = "horizontal",
                    legend.justification = "center") & harmonise_text(),
            PARAMS$supp_excess_size[["width"]], PARAMS$supp_excess_size[["height"]])

# %% [markdown]
# ## Fig. S16: within-carrier *de novo* burden, three index rules
#
# Figure 3D for DD-associated, pLI ≥ 0.9 and pLI < 0.9 genes separately (Figure
# 3D pools the first two), with each carrier's inherited PTVs summarised as the
# highest NMD-scaled PEXT (triangles), the total (filled circles, as in Figure
# 3D), or the total after removing carriers whose highest value is below 0.1
# (open circles).
#
# > **Per-individual data:** as Figure 3D, from `within_carrier_summary.tsv`.

# %%
if (!is.null(within_est))
  save_figure("figS16", panel_within(
    within_est, sets = PARAMS$within_supp_sets, rules = WITHIN_FIG_RULES, relabel = FALSE,
    xlab = paste("Per unit decrease in", PEXT_LAB, "of the carrier's inherited PTVs (95% CI)"),
    panel_in = PARAMS$panel_in[["within"]]) & harmonise_text(),
    PARAMS$supp_within_size[["width"]], PARAMS$supp_within_size[["height"]])

# %% [markdown]
# ## Fig. S18: PEXT of clinically interpreted PTVs against a permuted null
#
# Stratified label-permutation tests (10,000 permutations): the null is that the
# clinical label is independent of PEXT. Each permutation holds the per-study
# number of focal variants fixed and redraws that many variants from the same
# study's population. A: are pathogenic or likely pathogenic transmitted
# unrescued PTVs in DD-associated genes at higher NMD-scaled PEXT than chance?
# The population is every transmitted unrescued PTV in a DD-associated gene;
# unreported variants stay in as background. B: are clinically uncertain
# *de novo* PTVs in DD-associated genes at lower NMD-scaled PEXT? The population
# is the *de novo* PTVs labelled pathogenic/likely pathogenic, uncertain or
# benign/likely benign; unreported ones are dropped. P values are one-sided.
# Labels are the released `clinical_interpretation`: the most severe
# interpretation reported at the variant's position in the proband's clinical
# record.

# %%
perm_test <- function(dt, direction = c("greater", "less"), n_perm = PARAMS$n_perm,
                      stat_fun = mean, seed = PARAMS$seed) {
  direction <- match.arg(direction)
  dt <- as.data.table(dt)[!is.na(pext) & !is.na(focal) & !is.na(study)]
  if (!any(dt$focal)) stop("perm_test: no focal variants left after filtering.",
                           call. = FALSE)
  set.seed(seed)
  pext <- dt$pext; obs <- stat_fun(pext[dt$focal])
  idx_by_study <- split(seq_len(nrow(dt)), dt$study)
  k_by_study <- tapply(dt$focal, dt$study, sum)
  studies <- names(idx_by_study)
  null <- vapply(seq_len(n_perm), function(b) {
    sel <- unlist(lapply(studies, function(s) {
      ix <- idx_by_study[[s]]; k <- k_by_study[[s]]
      if (k == 0L) integer(0) else ix[sample.int(length(ix), k)]
    }), use.names = FALSE)
    stat_fun(pext[sel])
  }, numeric(1))
  p <- if (direction == "greater") (1 + sum(null >= obs)) / (n_perm + 1)
       else (1 + sum(null <= obs)) / (n_perm + 1)
  list(observed = obs, null = null, null_mean = mean(null), pop_mean = stat_fun(pext),
       p_value = p, direction = direction, n_focal = sum(dt$focal), n_pop = nrow(dt),
       n_perm = n_perm,
       by_study = data.table(study = studies,
                             n_focal = as.integer(k_by_study[studies]),
                             n_pop = as.integer(table(dt$study)[studies])))
}
format_perm_p <- function(res) {
  floor_p <- 1 / (res$n_perm + 1)
  if (res$p_value <= floor_p) paste0("< ", signif(floor_p, 2))
  else paste0("= ", signif(res$p_value, 3))
}
permutation_panel <- function(res, title_plain, tag) {
  d <- data.table(null = res$null)
  d[, extreme := if (res$direction == "greater") null >= res$observed
                 else null <= res$observed]
  rng <- range(c(res$null, res$observed)); span <- diff(rng)
  ## annotate on the side away from the observed line, so it never collides
  ## with the shaded tail
  x_ann <- if (res$direction == "greater") rng[1] + 0.02 * span else rng[2] - 0.02 * span
  h_ann <- if (res$direction == "greater") 0 else 1
  ggplot(d, aes(null, fill = extreme)) +
    geom_histogram(bins = 60, colour = NA) +
    geom_vline(xintercept = res$observed, colour = "black", linewidth = 1) +
    annotate("text", x = x_ann, y = Inf,
             label = sprintf("observed = %.3f\np %s\nn = %d", res$observed,
                             format_perm_p(res), res$n_focal),
             hjust = h_ann, vjust = 1.2, size = STYLE$annot_size, colour = "black",
             lineheight = 1.05) +
    scale_fill_manual(values = PALETTES$permutation, breaks = c("FALSE", "TRUE"),
                      name = NULL,
                      labels = c("Permuted null", "As or more extreme than observed")) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01),
                       expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = paste("Mean", PEXT_LAB, "of resampled variant set"), y = "Permutations",
         title = title_plain, tag = tag) +
    theme_pub(grid = "y") + theme(legend.position = "bottom")
}

# %%
## A: pathogenic/likely pathogenic INHERITED PTVs at higher PEXT than chance?
q1 <- copy(tbl[Transmitted == TRUE & class == PARAMS$ptv_class & gene_set == GS_DD])
q1[, focal := clinical_interpretation %in% "Pathogenic/likely pathogenic"]

## B: clinically 'Uncertain' DE NOVO PTVs at lower PEXT? Variants with no
## reported interpretation are dropped rather than treated as background. A
## variant on two DD-associated genes takes the mean of their PEXT.
q2 <- dnms_observed[gene_set == GS_DD,
                    .(pext = mean(pext, na.rm = TRUE),
                      clinical_interpretation = clinical_interpretation[1L]),
                    by = .(study, obs_id)]
q2 <- q2[clinical_interpretation %in% c("Pathogenic/likely pathogenic", "Uncertain",
                                        "Benign/likely benign")]
q2[, focal := clinical_interpretation == "Uncertain"]

if (any(q1$focal) && any(q2$focal)) {
  res_q1 <- perm_test(q1[, .(pext, focal, study)], direction = "greater")
  res_q2 <- perm_test(q2[, .(pext, focal, study)], direction = "less")
  save_figure("figS18", wrap_plots(
    permutation_panel(res_q1,
                      "Clinically interpreted\npathogenic/likely pathogenic inherited PTVs",
                      "A"),
    permutation_panel(res_q2, "Clinically interpreted\nuncertain de novo PTVs", "B"),
    nrow = 1, guides = "collect") & theme(legend.position = "bottom"),
    PARAMS$fig_perm_size[["width"]], PARAMS$fig_perm_size[["height"]])
} else msg("!! Fig. S18: no clinically interpreted focal variants; nothing written")

# %% [markdown]
# ## Figure N3: NMD-scaled PEXT against gene-level constraint (Supplementary Note 3)
#
# Unrescued PTVs and synonymous variants of the main table, with NMD-scaled
# PEXT and the sHet of their gene, on log scales (PEXT floored at 0.001). A:
# the exchange rate β by profile likelihood over a single-index model,
# logit P(transmitted) ~ f(log10 sHet + β · log10 NMD-scaled PEXT) × class, with
# f a natural spline; the sHet coefficient is fixed at 1, so β = 1 is the
# product axis. B: fitted excess transmission rate against NMD-scaled PEXT at
# four fixed sHet levels, from the full sHet × PEXT interaction model. C: the
# same split by DD association at two matched sHet values taken from the
# DD-associated genes.

# %%
N3_ACCENT <- "tomato4"

#' Unrescued PTVs and synonymous variants with NMD-scaled PEXT and the sHet of
#' their gene, on the log scales the models use.
note3_input <- function(tbl, shet) {
  dd <- merge(tbl[as.character(class) %in% c(PARAMS$syn_class, PARAMS$ptv_class),
                  .(vepSYMBOL, class, gene_set, pext, Transmitted, Relationship)],
              shet[, .(vepSYMBOL = GENE, shet)], by = "vepSYMBOL", all.x = FALSE)
  dd <- dd[!is.na(pext) & !is.na(shet) & shet > 0 & !is.na(Transmitted)]
  dd[, class := factor(as.character(class), levels = c(PARAMS$syn_class,
                                                       PARAMS$ptv_class))]
  dd[, Transmitted := as.logical(Transmitted)]
  dd[, ls := log10(shet)][, lp := log10(pext + PARAMS$pext_floor)]
  has_rel <- uniqueN(na.omit(dd$Relationship)) > 1L
  if (has_rel) dd[, Relationship := droplevels(factor(Relationship))]
  list(dd = dd, rel = if (has_rel) " + Relationship" else "")
}

#' The curves of Figure N3B and N3C.
note3_constraint <- function(dd, rel, df = PARAMS$spline_df) {
  has_rel <- nzchar(rel)
  fit <- function(rhs) glm(as.formula(paste0("Transmitted ~ ", rhs, rel)),
                           binomial("logit"), data = dd)
  m_int  <- fit(sprintf("ns(ls, %d) * ns(lp, %d) * class", df, df))

  ## B: fitted excess against PEXT at fixed constraint. Parallel curves mean the
  ##    PEXT effect does not depend on constraint; fanning means it does. The
  ##    sHet levels are quantiles over the genes carrying an unrescued PTV.
  ## C: the same, split by DD association at MATCHED constraint, so a gap
  ##    between the two is not explained by constraint.
  grid_p <- seq(0, 1, length.out = PARAMS$grid_n)
  by_gene <- dd[as.character(class) == PARAMS$ptv_class, .(shet = shet[1]), by = vepSYMBOL]
  q_ls <- quantile(by_gene$shet, c(0.25, 0.5, 0.9, 0.99), na.rm = TRUE)
  newdata_at <- function(shet_value, extra = NULL) {
    nd <- as.data.frame(CJ(pext = grid_p, class = levels(dd$class), sorted = FALSE))
    nd$class <- factor(nd$class, levels = levels(dd$class))
    nd$lp <- log10(nd$pext + PARAMS$pext_floor)
    nd$ls <- log10(unname(shet_value))
    if (has_rel)
      nd$Relationship <- factor(levels(dd$Relationship)[1], levels = levels(dd$Relationship))
    if (!is.null(extra)) for (k in names(extra)) nd[[k]] <- extra[[k]]
    nd
  }
  excess_at <- function(model, nd) {
    P <- predict_draws(model, nd, B = PARAMS$n_sim)
    ip <- which(nd$class == PARAMS$ptv_class); is_ <- which(nd$class == PARAMS$syn_class)
    100 * (P[, ip[order(nd$pext[ip])], drop = FALSE] -
             P[, is_[order(nd$pext[is_])], drop = FALSE])
  }
  curve_row <- function(E, nd) data.table(pext = sort(unique(nd$pext)),
                                          excess_pp = apply(E, 2, median),
                                          lo_pp = apply(E, 2, quantile, 0.025),
                                          hi_pp = apply(E, 2, quantile, 0.975))
  shet_curves <- rbindlist(lapply(seq_along(q_ls), function(i) {
    nd <- newdata_at(q_ls[i])
    cbind(data.table(shet_level = sprintf("sHet = %.3g (%s)", q_ls[i], names(q_ls)[i]),
                     shet = unname(q_ls[i])), curve_row(excess_at(m_int, nd), nd))
  }))
  shet_curves[, shet_level := factor(shet_level, levels = unique(shet_level))]

  ## DD-stratified. Rank deficiency is the expected failure: DD-associated genes
  ## occupy a narrow, high band of sHet, so interacting the sHet spline with the
  ## DD indicator can leave basis columns aliased. The fallback drops that
  ## interaction, which is not the term of interest -- what matters is whether
  ## the level and the PEXT slope differ at FIXED sHet.
  dd[, is_dd := factor(fifelse(as.character(gene_set) == GS_DD,
                               "DD-associated", "Not DD-associated"),
                       levels = c("Not DD-associated", "DD-associated"))]
  dd_curves <- NULL
  if (dd[is_dd == "DD-associated", .N] > 200L && uniqueN(dd$is_dd) > 1L) {
    m_dd <- tryCatch(fit(sprintf("(ns(ls, %d) + ns(lp, %d)) * class * is_dd", df, df)),
                     error = function(e) NULL)
    if (!is.null(m_dd) && anyNA(coef(m_dd))) {
      msg("!! the DD-stratified model is rank deficient (%d aliased term(s)); refitting without the sHet-by-DD interaction",
          sum(is.na(coef(m_dd))))
      m_dd <- tryCatch(fit(sprintf("ns(ls, %d) * class + ns(lp, %d) * class * is_dd",
                                   df, df)), error = function(e) NULL)
      if (!is.null(m_dd) && anyNA(coef(m_dd))) m_dd <- NULL
    }
    if (!is.null(m_dd)) {
      shet_dd  <- 10^dd[is_dd == "DD-associated", ls]
      q_dd <- quantile(shet_dd, c(0.5, 0.9), na.rm = TRUE)
      dd_curves <- rbindlist(lapply(seq_along(q_dd), function(i) {
        rbindlist(lapply(levels(dd$is_dd), function(g) {
          nd <- newdata_at(q_dd[i], list(is_dd = factor(g, levels = levels(dd$is_dd))))
          cbind(data.table(group = g,
                           shet_level = sprintf("sHet = %.3g (%s of DD genes)",
                                                q_dd[i], names(q_dd)[i])),
                curve_row(excess_at(m_dd, nd), nd))
        }), use.names = TRUE)
      }), use.names = TRUE)
      dd_curves[, group := factor(group, levels = levels(dd$is_dd))]
      dd_curves[, shet_level := factor(shet_level, levels = unique(shet_level))]
    }
  }
  dd[, is_dd := NULL]
  list(shet_curves = shet_curves, dd_curves = dd_curves)
}

#' Figure N3A: the exchange rate by profile likelihood over a single-index
#' model, logit P(transmitted) ~ f(log10 sHet + beta * log10 PEXT) * class.
#' The index is identified up to scale, so the sHet coefficient is fixed at 1
#' and beta is the exchange rate. beta = 1 is the product axis, tested rather
#' than assumed. beta is ATTENUATED by measurement error in whichever variable
#' is noisier -- PEXT almost certainly -- so it is a LOWER bound on the weight
#' of PEXT.
note3_exchange <- function(dd, rel, df = PARAMS$spline_df,
                           grid = PARAMS$exchange_grid_coarse,
                           step_fine = PARAMS$exchange_refine_step) {
  prof_dt <- dd[, intersect(c("Transmitted", "ls", "lp", "class", "Relationship"),
                            names(dd)), with = FALSE]
  f_idx <- as.formula(sprintf("Transmitted ~ ns(idx, %d) * class%s", df, rel))
  profile_over <- function(betas) {
    dev <- rep(NA_real_, length(betas))
    for (i in seq_along(betas)) {
      ## idx is overwritten in place: copying the table per beta ran out of memory
      set(prof_dt, j = "idx", value = prof_dt$ls + betas[i] * prof_dt$lp)
      m <- tryCatch(glm(f_idx, binomial("logit"), data = prof_dt, model = FALSE,
                        y = FALSE), error = function(e) NULL)
      if (!is.null(m)) dev[i] <- deviance(m)
      rm(m); if (i %% 5L == 0L) gc(verbose = FALSE)
    }
    data.table(beta = betas, deviance = dev)[is.finite(deviance)]
  }
  coarse <- profile_over(grid)
  prof <- coarse
  if (nrow(coarse)) {
    coarse[, delta_deviance := deviance - min(deviance)]
    step <- diff(grid)[1]
    near <- coarse[delta_deviance <= 6, range(beta)]
    fine <- seq(max(min(grid), near[1] - step), min(max(grid), near[2] + step),
                by = step_fine)
    fine <- fine[!round(fine, 6) %in% round(coarse$beta, 6)]
    if (length(fine)) prof <- rbindlist(list(coarse[, .(beta, deviance)],
                                             profile_over(fine)))
  }
  if (!nrow(prof)) return(NULL)
  setorder(prof, beta)
  prof[, delta_deviance := deviance - min(deviance)]
  cut <- qchisq(0.95, 1)
  inside <- prof[delta_deviance <= cut, beta]
  ## whether the interval ran into the edge of the grid, in which case it is a
  ## bound rather than an interval
  if (min(inside) <= min(prof$beta) || max(inside) >= max(prof$beta))
    msg("!! the profile interval reaches the edge of the beta grid; widen PARAMS$exchange_grid_coarse")
  list(profile = prof, interval = c(min(inside), max(inside)))
}

panel_n3_profile <- function(ex) {
  if (is.null(ex)) return(NULL)
  prof <- ex$profile
  lab <- sprintf("\u03b2 = %.2f [%.2f\u2013%.2f]", prof[which.min(deviance), beta],
                 ex$interval[1], ex$interval[2])
  ggplot(prof, aes(beta, delta_deviance)) +
    geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", colour = "grey50",
               linewidth = 0.5) +
    ## the product axis, tested rather than assumed
    geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50", linewidth = 0.5) +
    geom_line(linewidth = STYLE$line_width, colour = "grey15") +
    geom_point(data = prof[which.min(deviance)], colour = N3_ACCENT, size = 3) +
    annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.05, vjust = 1.4,
             size = STYLE$annot_small, colour = "grey20") +
    coord_cartesian(ylim = c(0, 20)) +
    labs(x = expression(beta * ": log units of NMD-scaled PEXT per log unit of sHet"),
         y = "Profile deviance\nabove minimum") +
    theme_pub(grid = "y") +
    theme(panel.grid.major.x = element_line(colour = STYLE$grid_colour,
                                            linewidth = STYLE$grid_width))
}

panel_n3_shet <- function(curves) {
  if (is.null(curves)) return(NULL)
  ggplot(curves, aes(pext, excess_pp, colour = shet_level, fill = shet_level)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lo_pp, ymax = hi_pp), alpha = 0.13, colour = NA) +
    geom_line(linewidth = STYLE$line_width) +
    scale_colour_viridis_d(option = "C", end = 0.85, name = NULL) +
    scale_fill_viridis_d(option = "C", end = 0.85, guide = "none") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    coord_cartesian(ylim = c(-3, 12)) +
    labs(x = "NMD-scaled PEXT", y = "Excess transmission\nrate (pp)") +
    theme_pub(grid = "y")
}

panel_n3_dd <- function(curves) {
  if (is.null(curves)) return(NULL)
  cols <- c(`Not DD-associated` = "grey45", `DD-associated` = GENE_SET_COLOURS[[GS_DD]])
  ggplot(curves, aes(pext, excess_pp, colour = group, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lo_pp, ymax = hi_pp), alpha = 0.15, colour = NA) +
    geom_line(linewidth = STYLE$line_width) +
    facet_wrap(~ shet_level, nrow = 1) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = "NMD-scaled PEXT", y = "Excess transmission\nrate (pp)") +
    theme_pub(grid = "y")
}

figure_n3 <- function(pa, pb, pc) {
  if (is.null(pa) && is.null(pb) && is.null(pc)) return(NULL)
  blank <- patchwork::plot_spacer()
  or_blank <- function(p) if (is.null(p)) blank else p
  wrap_plots(A = or_blank(pa), B = or_blank(pb), C = or_blank(pc), design = "AB\nCC") +
    plot_annotation(tag_levels = "A") & harmonise_text()
}

# %%
## sHet of each gene (Weghorn et al. 2019), from the gene annotation
shet_scores <- unique(fread(PATHS$genes, encoding = "UTF-8")[
  !is.na(sHet), .(GENE = gene, shet = sHet)])
n3_in <- note3_input(tbl, shet_scores)
if (nrow(n3_in$dd) >= 500) {
  n3 <- note3_constraint(n3_in$dd, n3_in$rel)
  n3_ex <- note3_exchange(n3_in$dd, n3_in$rel)
  fig <- figure_n3(panel_n3_profile(n3_ex), panel_n3_shet(n3$shet_curves),
                   panel_n3_dd(n3$dd_curves))
  if (!is.null(fig)) save_figure("figN3", fig, PARAMS$fig_n3_size[["width"]],
                                 PARAMS$fig_n3_size[["height"]])
} else msg("!! Supplementary Note 3: fewer than 500 variants with PEXT and sHet; nothing written")

# %%
sessionInfo()

# %% [markdown]
# # 04 · Incompletely penetrant PTVs occur in most haploinsufficient DD-associated genes
#
# This notebook tests the transmission of unrescued protein-truncating variants (PTVs)
# gene by gene against each gene set's synonymous transmission rate, and asks how the
# transmission excess is spread across genes. Every quantity is a per-gene count of
# parental variant records; nothing is pooled within individuals.
#
# - **Figure 4A**: the haploinsufficient DDG2P genes tested, and how many carry at least
#   one unrescued parental PTV in the cohort
# - **Figure 4B**: per-gene transmission excess (one-sided binomial P) against a simulated
#   null, with a matched set of unconstrained genes
# - **Figure 4C**: cumulative share of the transmission excess, genes ranked by their
#   leave-one-out influence
# - **Fig. S17**: DD-associated genes with an excess of at least two transmitted PTVs,
#   coloured by published evidence of incomplete penetrance

# %%
REPO_DIR <- Sys.getenv("RR_REPO", unset = if (file.exists("config.R")) "." else "..")
source(file.path(REPO_DIR, "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

# %% [markdown]
# ## Inputs
#
# - `parental_variants.tsv.gz` (EGA release): the parents' PTV and synonymous calls, one
#   row per variant, transcript and parent, with transmission to the proband, the
#   RescueRanger prediction (`rescue`) and the gene set. Hypermutated trios, call quality
#   filters and the heterozygous-in-proband check on transmitted variants were applied
#   before release (`data/README.md`).
# - `gene_annotation.tsv.gz`: the DDG2P mechanism class of every DDG2P gene.
# - The curated table of published penetrance evidence that colours Fig. S17.

# %%
PATHS <- list(
  parents         = file.path(DATA_DIR, "parental_variants.tsv.gz"),
  gene_annotation = file.path(RESOURCE_DIR, "gene_annotation.tsv.gz"),
  ## curated by hand from the literature (Methods: Concentration of the excess)
  penetrance      = file.path(RESOURCE_DIR, "gene_penetrance_consolidated.tsv"))

# %% [markdown]
# ## Parameters, helpers and plotting theme

# %%
PARAMS <- list(
  ptv_class = "Unrescued PTVs", syn_class = "Synonymous",
  ## RECURRENCE WEIGHTING for the per-gene scan. A rare PTV seen in k unrelated
  ## trios is usually one ancestral mutation on a shared haplotype, not k
  ## independent inheritances. A variant seen k times contributes k^alpha.
  ## BOTH sides are weighted: transmission is a ratio, so downweighting the
  ## denominator alone would change the rate rather than the evidence.
  jackknife_recurrence_alpha = 0.5,
  ## Minimum PTV CALLS for a gene to enter the scan. 1 means every gene in which
  ## the statistic is defined; the scan also sets the exome-wide BH correction.
  scan_min_ptv = 1L,
  ## Display-only threshold for Figure 4B: a one-sided binomial on n = 1 can
  ## return only two distinct p-values, so those genes draw a step function at
  ## the bottom-left.
  qq_min_ptv = 3L,
  qq_n_sim = 2000L, qq_control_draws = 200L,
  ## Fig. S17 draws every DD-associated gene whose own excess (transmitted PTVs
  ## minus those expected at the synonymous rate) is at least this many PTVs.
  ## A threshold rather than a count, so genes tied at the edge are never split
  ## by name; ties are common because the excess is a function of two integers.
  s17_min_excess = 2,
  seed = 1)
set.seed(PARAMS$seed)

msg <- function(...) message(sprintf(...))

FONT_MS <- 13; FONT_TAG <- 15; FONT_ANNOT_SMALL <- 3.0
ACCENT <- "tomato4"
HAVE_GGTEXT <- requireNamespace("ggtext", quietly = TRUE)
element_md <- function(...) {
  if (HAVE_GGTEXT) ggtext::element_markdown(...) else element_text(...)
}

theme_pub <- function(base = FONT_MS, grid = c("y", "x", "both", "none")) {
  grid <- match.arg(grid)
  th <- theme(text = element_text(family = "sans", size = base, colour = "black"),
              panel.background = element_blank(), plot.background = element_blank(),
              panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linewidth = 0.9),
              axis.line.x.top = element_blank(), axis.line.y.right = element_blank(),
              axis.ticks = element_line(colour = "black", linewidth = 0.5),
              axis.ticks.length = unit(3, "pt"),
              axis.text = element_text(size = base, colour = "grey10"),
              axis.title = element_text(size = base),
              legend.position = "bottom", legend.box = "horizontal",
              legend.title = element_blank(), legend.text = element_text(size = base),
              legend.key = element_blank(), strip.background = element_blank(),
              strip.text = element_text(size = base, face = "bold"),
              plot.title = element_md(size = base, face = "bold"),
              plot.subtitle = element_md(size = base - 2, colour = "grey30"),
              plot.tag = element_text(face = "bold", size = FONT_TAG),
              plot.tag.position = c(0.002, 1), plot.margin = margin(6, 20, 6, 6))
  gl <- element_line(colour = "grey88", linewidth = 0.4)
  if (grid %in% c("y", "both")) th <- th + theme(panel.grid.major.y = gl)
  if (grid %in% c("x", "both")) th <- th + theme(panel.grid.major.x = gl)
  th
}

## The mechanisms TESTED: a transmission excess bears on penetrance under
## dominant loss of function, and under dual-mode, where the same excess could
## equally come from compound heterozygous probands. Both are tested; the
## distinction is carried in the mechanism column.
QQ_KEEP_MECHANISMS <- c("Dominant loss of function", "Dominant LoF and recessive")

QQ_MECH_COLS <- c("Dominant loss of function" = "tomato4",
                  ## dark orange: distinct from the dominant genes' red beside it
                  ## (Figure 4A and 4B); Figure 4C's gold is a separate panel and key
                  "Dominant LoF and recessive" = "darkorange3",
                  "Dominant LoF and other mechanism" = "#8C6BB1",
                  "Dominant LoF, unconstrained" = "#D6604D",
                  "Dominant, other mechanism" = "#5B9BD5",
                  "Recessive" = "#7B1E12", "Other or undetermined" = "grey60")
QQ_CONTROL_LABEL <- "Unconstrained, non-G2P (matched)"
QQ_CONTROL_COL   <- "skyblue3"
QQ_SERIES_COLS   <- c(QQ_MECH_COLS, setNames(QQ_CONTROL_COL, QQ_CONTROL_LABEL))
TIER_LABELS <- c(documented = "Documented reduced penetrance", uncertain = "Uncertain",
                 none = "No verified evidence identified")
TIER_COLS   <- c(documented = ACCENT, uncertain = "#C1852B", none = "grey55")

# %% [markdown]
# ## Analysis table
#
# Rare parental variants in three classes: unrescued PTVs (LOFTEE high confidence, not
# predicted to be rescued), rescued PTVs and synonymous variants. The transcript rows of a
# variant in one parent are one record when they agree on gene, class, consequence,
# LOFTEE call and PEXT (on rows placed in an NMD region, NMD-scaled PEXT divided by the
# region's NMD efficiency, as in notebook 03). A variant with more than one record in a
# family (`obs_id`) is removed: either both parents carry it, so its parent of origin is
# undefined, or its transcript rows disagree.

# %%
g2p_mechanism <- unique(fread(PATHS$gene_annotation, encoding = "UTF-8")[
  !is.na(ddg2p_mechanism), .(gene, mechanism = ddg2p_mechanism)])

## NMD efficiency of each region relative to NMD-target (Methods; as in
## notebook 03), to recover PEXT from NMD-scaled PEXT
ASE_MEAN_AF <- c(Synonymous = 0.487, `NMD-target` = 0.383, `Start-proximal` = 0.461,
                 `Long exon` = 0.438, Distal = 0.486)
nmd_efficiency <- function(af) {
  odds <- function(f) f / (1 - f)
  eff <- 1 - odds(af[setdiff(names(af), "Synonymous")]) / odds(af[["Synonymous"]])
  eff / eff[["NMD-target"]]
}
NMD_ETA <- nmd_efficiency(ASE_MEAN_AF)

build_table <- function() {
  parents <- fread(PATHS$parents, colClasses = list(character = "CHROM"),
                   encoding = "UTF-8")
  parents[, pext := signif(nmd_scaled_pext / unname(NMD_ETA[nmd_region]), 10)]
  classes <- rbindlist(list(
    parents[is_ptv == TRUE & rescue %in% "rescued"][, class := "Rescued PTVs"],
    parents[is_ptv == TRUE & rescue %in% "not_rescued" & vepLoF == "HC"][
      , class := PARAMS$ptv_class],
    parents[is_ptv == FALSE][, class := PARAMS$syn_class]),
    use.names = TRUE, fill = TRUE)
  tbl <- unique(classes[, .(study, obs_id, CHROM, POS, REF, ALT, vepSYMBOL, Relationship,
                            Transmitted, class, gene_set, rescue, vepConsequence, vepLoF,
                            pext)])
  tbl[, n_records := .N, by = obs_id]
  tbl <- tbl[n_records == 1]
  tbl[, gene_set := factor(gene_set, levels = GENE_SETS)]
  tbl[, class := factor(class, levels = c(PARAMS$syn_class, PARAMS$ptv_class,
                                          "Rescued PTVs"))]
  tbl[]
}

tbl <- build_table()

# ## Per-gene transmission test
#
# For every gene with at least one unrescued PTV call, transmitted PTVs are compared with
# the gene set's synonymous transmission rate by a one-sided binomial test (Methods:
# Quantifying per-gene transmission excess). A variant observed in *k* unrelated trios
# contributes *√k* in total, each call carrying 1/*k* of it, on both the transmitted and
# the total count. P values are Benjamini-Hochberg corrected over the DDG2P genes tested
# (dominant loss of function, and dominant loss of function with a recessive entry) and,
# separately, exome-wide.
#
# *k* is the number of records of the same variant (chromosome, position, alleles) among
# unrescued PTVs: in the analysis table each carrier of a variant is one record.

# %%
d <- tbl[as.character(class) %in% c(PARAMS$syn_class, PARAMS$ptv_class) &
           !is.na(Transmitted)]
d[, gene_set := as.character(gene_set)]

syn_rate <- d[as.character(class) == PARAMS$syn_class,
              .(p_syn = mean(as.logical(Transmitted)), n_syn = .N), by = gene_set]
pt <- d[as.character(class) == PARAMS$ptv_class]
alpha <- PARAMS$jackknife_recurrence_alpha
pt[, vkey := paste(CHROM, POS, REF, ALT, sep = ":")]
pt[, n_trios := .N, by = vkey]
## each CALL carries 1/k of the variant's weight, so a variant totals k^alpha
pt[, w := n_trios^alpha / n_trios]

scan <- pt[, .(n_ptv = sum(w), n_transmitted = sum(w * as.logical(Transmitted)),
               n_ptv_raw = .N, n_transmitted_raw = sum(as.logical(Transmitted))),
           by = .(gene_set, gene = vepSYMBOL)]
scan <- merge(scan, syn_rate, by = "gene_set")
scan <- scan[n_ptv_raw >= PARAMS$scan_min_ptv]
## binom.test needs integers, so the weighted counts are rounded for the test only
scan[, k_test := pmin(round(n_transmitted), round(n_ptv))]
scan[, n_test := pmax(round(n_ptv), 1)]
scan[, p_value := mapply(function(k, n, p)
  binom.test(k, n, p, alternative = "greater")$p.value, k_test, n_test, p_syn)]
## TWO CORRECTIONS, over the two families that mean something here.
##
##   q_tested   BH over the G2P genes tested -- dominant loss-of-function plus
##              dual-mode. This is the family the per-gene claims belong to.
##   q_exome    BH over every gene carrying a PTV call, the conservative bound.
##
## There is deliberately no correction within the constraint split: the genes
## tested are a curated set, not a constraint class.
tested_genes <- g2p_mechanism[mechanism %in% QQ_KEEP_MECHANISMS, gene]
scan[, tested := gene %in% tested_genes]
scan[, q_tested := NA_real_]
scan[tested == TRUE, q_tested := p.adjust(p_value, "BH")]
scan[, q_exome := p.adjust(p_value, "BH")]
scan[, transmission_rate := n_transmitted / n_ptv]
setorder(scan, q_exome, -n_ptv)

scan <- merge(scan, g2p_mechanism, by = "gene", all.x = TRUE)
scan[is.na(mechanism), mechanism := "Other or undetermined"]

# %% [markdown]
# ## Figure 4B: per-gene transmission excess against a simulated null
#
# Observed -log10(P) for the DDG2P genes tested, against the quantiles expected under each
# gene's synonymous transmission rate. Expected quantiles and the grey 95% envelope come
# from 2,000 simulations in which transmissions are redrawn as Binomial(*n* PTVs,
# synonymous rate) over exactly the genes drawn; the envelope is simultaneous, so it
# describes whole null curves rather than single ranks. Blue points are unconstrained genes
# (pLI < 0.9) with no DDG2P entry, matched one-to-one on number of PTV calls and shown as
# the median of 200 draws. The panel shows genes with at least three PTV calls; the seven
# genes with the smallest P are labelled, in bold if significant after BH correction over
# the genes tested.

# %%
## Expected quantiles come from SIMULATION rather than i/(m+1). These are
## one-sided BINOMIAL tests: a gene with 5 PTVs can produce only six possible
## p-values, so the null is discrete and conservative, and plotting against the
## uniform would show deflation that is an artefact of discreteness.
##
## The null is re-simulated over exactly the genes DRAWN. Order statistics
## spread with m, so a null simulated over every G2P gene has larger upper
## quantiles than one over the few hundred plotted.
qq_resimulate <- function(dd, n_sim = PARAMS$qq_n_sim, seed = PARAMS$seed) {
  m <- nrow(dd); if (m < 3L) return(NULL)
  p <- dd$p_syn; n <- pmax(as.integer(round(dd$n_ptv)), 1L)
  set.seed(seed)
  null_mat <- vapply(seq_len(n_sim), function(b) {
    k <- rbinom(m, n, p)
    sort(-log10(pmax(pbinom(k - 1L, n, p, lower.tail = FALSE), .Machine$double.xmin)))
  }, numeric(m))
  expected <- rowMeans(null_mat)

  ## A SIMULTANEOUS ENVELOPE, NOT A POINTWISE BAND. A pointwise 95% band on
  ## these ranks is a hairline, because order statistics of a fixed-size sample
  ## are pinned very tightly. The pointwise level is widened until 95% of
  ## complete null curves lie entirely inside, so "the observed curve leaves the
  ## envelope" is a 5%-level statement about the CURVE.
  simultaneous_env <- function(M, target = 0.95, alpha0 = 0.05) {
    a <- alpha0
    lo <- apply(M, 1, quantile, a / 2); hi <- apply(M, 1, quantile, 1 - a / 2)
    for (i in seq_len(80)) {
      if (mean(colSums(M < lo | M > hi) == 0) >= target) break
      a <- a / 1.15
      lo <- apply(M, 1, quantile, a / 2); hi <- apply(M, 1, quantile, 1 - a / 2)
    }
    list(lo = lo, hi = hi, pointwise_alpha = a,
         coverage = mean(colSums(M < lo | M > hi) == 0))
  }
  env <- simultaneous_env(null_mat)
  list(expected = expected, lo = env$lo, hi = env$hi,
       pointwise_alpha = env$pointwise_alpha, coverage = env$coverage,
       n_sim = n_sim, m = m)
}

#' Median observed curve for unconstrained non-G2P genes matched to `dd` on PTV
#' count. Matching matters: unconstrained genes carry more PTVs per gene, so an
#' unmatched draw would have more power, reach smaller p-values and extend
#' further along the expected axis -- the contrast would be about power. The
#' plotted curve is the MEDIAN over many matched draws, not one lucky draw.
qq_control_series <- function(dd, pool_scan, n_draws = PARAMS$qq_control_draws,
                              seed = PARAMS$seed) {
  pool <- pool_scan[gene_set == GS_LOW & !gene %in% dd$gene &
                      mechanism == "Other or undetermined"]
  if (nrow(pool) < 100) { msg("!! control pool is only %d genes; skipped", nrow(pool))
    return(NULL) }
  target <- dd$n_ptv_raw
  by_n <- split(seq_len(nrow(pool)), pool$n_ptv_raw)
  avail <- as.numeric(names(by_n))
  cand <- lapply(target, function(v) {
    g <- by_n[[as.character(v)]]
    if (is.null(g)) g <- by_n[[as.character(avail[which.min(abs(avail - v))])]]
    g
  })
  m <- length(target); set.seed(seed)
  curves <- matrix(NA_real_, m, n_draws)
  for (r in seq_len(n_draws)) {
    i <- vapply(cand, function(g) g[sample.int(length(g), 1L)], integer(1))
    x <- pool[i]
    n <- pmax(as.integer(round(x$n_ptv)), 1L)
    k <- pmin(as.integer(round(x$n_transmitted)), n)
    curves[, r] <- sort(-log10(pmax(pbinom(k - 1L, n, x$p_syn, lower.tail = FALSE),
                                    .Machine$double.xmin)))
  }
  list(observed = apply(curves, 1, median), n = m, n_draws = n_draws,
       n_pool = nrow(pool))
}

# %%
## The mechanism filter defines the genes TESTED; the PTV threshold defines the
## genes DRAWN. Multiple testing is corrected over the former: which points end
## up on the panel is a display decision and cannot move a q-value.
d_tested <- scan[mechanism %in% QQ_KEEP_MECHANISMS]
d_tested[, q_tested := p.adjust(p_value, "BH")]

d_drawn <- copy(d_tested[n_ptv_raw >= PARAMS$qq_min_ptv])
setorder(d_drawn, p_value)
d_drawn <- d_drawn[order(-p_value)]                    # least significant first
d_drawn[, observed := -log10(pmax(p_value, .Machine$double.xmin))]
sim <- qq_resimulate(d_drawn)
if (!is.null(sim))
  d_drawn[, expected := sim$expected][, lo := sim$lo][, hi := sim$hi]
ctrl <- qq_control_series(d_drawn, scan)

lab <- d_drawn[order(-observed)][seq_len(min(7L, .N))]
lab[, lab_expr := fifelse(q_tested < 0.1, sprintf('bolditalic("%s")', gene),
                          sprintf('italic("%s")', gene))]

pts <- rbind(
  if (!is.null(ctrl)) data.table(expected = d_drawn$expected,
                                 observed = ctrl$observed, series = QQ_CONTROL_LABEL),
  data.table(expected = d_drawn$expected, observed = d_drawn$observed,
             series = as.character(d_drawn$mechanism)), use.names = TRUE)
pts[, series := factor(series, levels = c(names(QQ_MECH_COLS), QQ_CONTROL_LABEL))]
pts <- pts[!is.na(series)]

fig4b <- ggplot(d_drawn, aes(expected, observed)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey78", alpha = 0.55) +
  geom_abline(slope = 1, intercept = 0, colour = "grey35", linewidth = 0.5)
if (!is.null(ctrl))
  ## the control is drawn first so the tested genes sit on top where they
  ## overlap, and with a border in its own colour: with colour = NA the points
  ## render in the legend key but not in the panel
  fig4b <- fig4b + geom_point(data = pts[series == QQ_CONTROL_LABEL],
                              aes(expected, observed, fill = series),
                              inherit.aes = FALSE, shape = 21,
                              colour = QQ_CONTROL_COL, stroke = 0.3, size = 1.5,
                              alpha = 0.8)
fig4b <- fig4b +
  geom_point(data = pts[series %in% names(QQ_MECH_COLS)],
             aes(expected, observed, colour = series), inherit.aes = FALSE,
             size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = QQ_SERIES_COLS, name = NULL, drop = TRUE) +
  scale_fill_manual(values = QQ_SERIES_COLS, name = NULL, drop = TRUE) +
  ## two guides rather than one two-row guide: guide_legend lays keys out on a
  ## grid, so a lone third key sits in the first column and reads as hanging off
  guides(colour = guide_legend(order = 1, nrow = 1,
                               override.aes = list(size = 2.6, alpha = 1)),
         fill = guide_legend(order = 2, nrow = 1,
                             override.aes = list(size = 2.6, alpha = 1, shape = 21,
                                                 colour = QQ_CONTROL_COL))) +
  labs(x = "Expected -log10(P) under the synonymous transmission rate",
       y = "Observed -log10(P)") +
  theme_pub(grid = "both") +
  theme(legend.position = "top", legend.box = "vertical", legend.box.just = "center",
        legend.justification = "center", legend.spacing.y = unit(2, "pt"))
fig4b <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  fig4b + ggrepel::geom_text_repel(data = lab, aes(label = lab_expr), parse = TRUE,
                                   size = FONT_ANNOT_SMALL, colour = "grey20",
                                   max.overlaps = 20, min.segment.length = 0.2,
                                   seed = 1)
} else {
  fig4b + geom_text(data = lab, aes(label = lab_expr), parse = TRUE,
                    size = FONT_ANNOT_SMALL, colour = "grey20", hjust = -0.15)
}

# %% [markdown]
# ## Figure 4C: cumulative share of the transmission excess
#
# Each gene's contribution is its leave-one-out influence on its gene set's transmission
# excess: the excess (transmitted unrescued PTVs above the number expected at the
# synonymous rate) with and without the gene, removing both its PTVs and its synonymous
# variants. Genes are ranked by contribution, and the cumulative share of the total is
# drawn for DD-associated and pLI ≥ 0.9 genes; dashed lines show the share expected if
# every PTV contributed equally (genes ranked from most to fewest PTVs).

# %%
## Leave-one-out is computed analytically from per-gene counts rather than by
## re-filtering once per gene: the estimator is a difference of two proportions,
## so dropping a gene is subtracting its four counts.
rate_of <- function(tp, np, ts, ns)
  fifelse(np <= 0 | ns <= 0, NA_real_, 100 * (tp / np - ts / ns))

jack_one <- function(gs) {
  x <- d[gene_set == gs]
  ptv <- x[as.character(class) == PARAMS$ptv_class]
  syn <- x[as.character(class) == PARAMS$syn_class]
  if (!nrow(ptv) || !nrow(syn)) return(NULL)
  Tp <- sum(as.logical(ptv$Transmitted)); Np <- nrow(ptv)
  Ts <- sum(as.logical(syn$Transmitted)); Ns <- nrow(syn)
  full_rate <- rate_of(Tp, Np, Ts, Ns); full_excess <- Np * full_rate / 100
  gp  <- ptv[, .(t_ptv = sum(as.logical(Transmitted)), n_ptv = .N),
             by = .(gene = vepSYMBOL)]
  gsn <- syn[, .(t_syn = sum(as.logical(Transmitted)), n_syn = .N),
             by = .(gene = vepSYMBOL)]
  out <- merge(gp, gsn, by = "gene", all.x = TRUE)
  out[is.na(t_syn), t_syn := 0L][is.na(n_syn), n_syn := 0L]
  out[, rate_loo_pp := rate_of(Tp - t_ptv, Np - n_ptv, Ts - t_syn, Ns - n_syn)]
  out[, excess_loo := (Np - n_ptv) * rate_loo_pp / 100]
  out[, `:=`(gene_set = gs, n_transmitted = t_ptv, rate_full_pp = full_rate,
             excess_full = full_excess, delta_pp = full_rate - rate_loo_pp,
             delta_excess = full_excess - excess_loo)]
  out[, pct_of_excess := 100 * delta_excess / excess_full]
  out[, .(gene_set, gene, n_ptv, n_transmitted, n_syn, rate_full_pp, rate_loo_pp,
          delta_pp, excess_full, excess_loo, delta_excess, pct_of_excess)]
}
FIG4C_SETS <- c(GS_DD, GS_HIGH)
jack <- rbindlist(lapply(FIG4C_SETS, jack_one))

per_set <- Filter(Negate(is.null), lapply(FIG4C_SETS, function(gs) {
  x <- jack[gene_set == gs]
  if (!nrow(x)) return(NULL)
  v <- sort(pmax(x$delta_excess, 0), decreasing = TRUE)
  if (!length(v) || sum(v) <= 0) return(NULL)
  cum <- data.table(gene_set = gs, rank = seq_along(v), cum = 100 * cumsum(v) / sum(v))
  ## expected under equal contribution per PTV: genes ranked from most to fewest
  ## PTVs, cumulative share of the PTV total
  nv <- sort(x$n_ptv, decreasing = TRUE)
  cum[, expected := 100 * cumsum(nv) / sum(nv)][, n_genes := length(v)][]
}))
cum <- rbindlist(per_set)
cum[, gene_set := factor(gene_set, levels = FIG4C_SETS)]
n50 <- cum[cum >= 50, .(n50 = min(rank), n_genes = n_genes[1]), by = gene_set]
marks <- rbindlist(lapply(FIG4C_SETS, function(gs) {
  x <- cum[gene_set == gs]
  if (!nrow(x)) return(NULL)
  x[rank %in% unique(c(1, 5, 20, n50[gene_set == gs, n50]))]
}))
marks[, above := gene_set == FIG4C_SETS[1]]
marks[, lab := sprintf("%d gene%s\n%.0f%%", rank, fifelse(rank == 1, "", "s"), cum)]
brks <- sort(unique(c(1, 5, 20, 100, n50$n_genes)))

fig4c <- ggplot(cum, aes(rank, cum, colour = gene_set)) +
  geom_line(aes(y = expected, group = gene_set), linetype = "dashed", linewidth = 0.5,
            alpha = 0.55, show.legend = FALSE) +
  geom_line(linewidth = 1.05, show.legend = FALSE) +
  geom_point(data = marks, size = 1.9, show.legend = FALSE) +
  geom_pointrange(data = data.frame(.x = NA_real_, .y = NA_real_,
                                    gene_set = factor(FIG4C_SETS,
                                                      levels = FIG4C_SETS)),
                  aes(x = .x, y = .y, ymin = .y, ymax = .y, colour = gene_set),
                  inherit.aes = FALSE, na.rm = TRUE, show.legend = TRUE) +
  geom_text(data = marks[above == TRUE], aes(label = lab), vjust = -0.55, hjust = 0.5,
            size = FONT_ANNOT_SMALL, lineheight = 0.9, show.legend = FALSE) +
  geom_text(data = marks[above == FALSE], aes(label = lab), vjust = 1.45, hjust = 0.5,
            size = FONT_ANNOT_SMALL, lineheight = 0.9, show.legend = FALSE) +
  scale_colour_manual(values = GENE_SET_COLOURS, name = NULL, drop = TRUE,
                      limits = FIG4C_SETS) +
  scale_x_continuous(trans = "log10", breaks = brks, labels = comma_format(accuracy = 1),
                     expand = expansion(mult = c(0.13, 0.13))) +
  scale_y_continuous(limits = c(-16, 116), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Genes, ranked by leave-one-out influence\non transmission excess (log scale)",
       y = "Cumulative share of\ntransmission excess") +
  theme_pub(grid = "both") +
  theme(legend.position = "top", legend.justification = "center",
        legend.background = element_blank(), legend.margin = margin(2, 4, 2, 4)) +
  guides(colour = guide_legend(nrow = 1))

# %% [markdown]
# ## Figure 4A: the DDG2P genes tested
#
# The two DDG2P sets tested, as one stacked bar: dominant loss-of-function genes (the
# DD-associated set used throughout the paper) and the genes that also carry a recessive
# entry. Within each set, genes carrying at least one unrescued parental PTV in the cohort
# (the genes drawn in Figure 4B) are striped, the rest solid.

# %%
GENE_BAR_SETS <- c(`Dominant loss of function` = "Dominant DDG2P genes (DD-associated)",
                   `Dominant LoF and recessive` = "Dominant and recessive DDG2P genes")
gene_bar_counts <- function() {
  rbindlist(lapply(names(GENE_BAR_SETS), function(m) {
    in_set <- g2p_mechanism[mechanism == m, unique(gene)]
    data.table(mechanism = m, label = GENE_BAR_SETS[[m]], n_genes = length(in_set),
               n_carrying = scan[mechanism == m & n_ptv_raw >= 1, uniqueN(gene)])
  }))
}
## stripes are drawn with clipped segments, so no pattern package is needed
hatch_in_rects <- function(r, slope, spacing) {
  rbindlist(lapply(seq_len(nrow(r)), function(i) {
    x <- r[i]
    a <- seq(x$x0 - x$ytop / slope, x$x1 - x$ybot / slope, by = spacing)
    xl <- pmax(x$x0, a + x$ybot / slope); xr <- pmin(x$x1, a + x$ytop / slope)
    ok <- xr > xl
    if (!any(ok)) return(NULL)
    data.table(x = xl[ok], xend = xr[ok], y = slope * (xl[ok] - a[ok]),
               yend = slope * (xr[ok] - a[ok]), colour = x$colour)
  }))
}
#' `width_in`: the panel's width in the assembled figure; `y_unit_in`: roughly
#' how many inches one y unit spans. Both only set the stripes' angle (~45 degrees).
panel_gene_sets <- function(cnt, width_in = 11, y_unit_in = 0.8) {
  cnt <- copy(cnt)[n_genes > 0][, colour := unname(QQ_MECH_COLS[mechanism])]
  if (!nrow(cnt)) return(NULL)
  bar_top <- 0.5; line_y <- 0.8           # the bar, and the line its total sits on
  seg <- rbindlist(lapply(seq_len(nrow(cnt)), function(i) {
    x <- cnt[i]
    data.table(mechanism = x$mechanism, colour = x$colour, carrying = c(TRUE, FALSE),
               n = c(x$n_carrying, x$n_genes - x$n_carrying))
  }))
  seg[, x1 := cumsum(n)][, x0 := x1 - n]
  seg[, `:=`(ybot = 0, ytop = bar_top)]
  total <- sum(cnt$n_genes)
  ## the stripes are lines of the set's colour over a pale tint of it
  tint <- function(col) colorRampPalette(c("white", col))(10)[3]
  seg[, fill := fifelse(carrying, vapply(colour, tint, ""), colour)]
  slope <- (width_in / total) / y_unit_in
  hatch <- hatch_in_rects(seg[carrying == TRUE & n > 0], slope = slope,
                          spacing = total / 120)
  ## counts inside each piece, bold at the axis-title size
  seg[, lab := fifelse(n > 0.04 * total, prettyNum(n, big.mark = ","), "")]
  ## each set's total centred over it, on a line in the set's colour that runs
  ## out to the set's edges; the ends are pulled in slightly so neighbouring
  ## sets' lines do not join
  sets <- seg[, .(x0 = min(x0), x1 = max(x1), colour = colour[1]), by = mechanism]
  sets <- merge(sets, cnt[, .(mechanism, n_genes)], by = "mechanism")
  sets[, title := sprintf("%s genes", prettyNum(n_genes, big.mark = ","))]
  inset <- 0.005 * total
  sets[, `:=`(xa = x0 + inset, xb = x1 - inset)]
  ggplot(seg) +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = ybot, ymax = ytop, fill = I(fill))) +
    geom_segment(data = hatch, aes(x = x, xend = xend, y = y, yend = yend, colour = I(colour)),
                 linewidth = 0.55) +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = ybot, ymax = ytop), fill = NA,
              colour = "black", linewidth = 0.35) +
    geom_label(aes(x = (x0 + x1) / 2, y = bar_top / 2, label = lab),
               data = seg[carrying == TRUE & lab != ""], size = FONT_MS / .pt,
               fontface = "bold", fill = "white", label.size = 0,
               label.padding = unit(1.5, "pt")) +
    geom_text(aes(x = (x0 + x1) / 2, y = bar_top / 2, label = lab),
              data = seg[carrying == FALSE & lab != ""], size = FONT_MS / .pt,
              colour = "white", fontface = "bold") +
    geom_segment(data = sets, aes(x = xa, xend = xb, y = line_y, yend = line_y,
                                  colour = I(colour)), linewidth = 0.7) +
    geom_label(data = sets, aes(x = (x0 + x1) / 2, y = line_y, label = title,
                                colour = I(colour)),
               size = FONT_MS / .pt, fill = "white", label.size = 0,
               label.padding = unit(4, "pt")) +
    scale_x_continuous(expand = expansion(mult = c(0.005, 0.005)), labels = comma) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
    labs(x = "DDG2P genes (striped: gene carries at least one unrescued parental PTV in our cohort)",
         y = NULL) +
    theme_pub(grid = "none") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
}
gb <- gene_bar_counts()
## Figure 4C's DD-associated genes are exactly Figure 4A's dominant genes with a PTV
n_dd_jack <- jack[gene_set == GS_DD, uniqueN(gene)]
if (gb[mechanism == "Dominant loss of function", n_carrying] != n_dd_jack)
  msg("!! Figure 4A counts %d dominant genes carrying a PTV, Figure 4C %d DD-associated genes",
      gb[mechanism == "Dominant loss of function", n_carrying], n_dd_jack)
fig4a <- panel_gene_sets(gb)

# %% [markdown]
# ## Figure 4: assembly

# %%
save_figure("fig4", (if (is.null(fig4a)) plot_spacer() else fig4a) / (fig4b | fig4c) +
              plot_layout(heights = c(0.24, 1)) +
              plot_annotation(tag_levels = "A") &
              theme(plot.tag = element_text(face = "bold", size = FONT_TAG)), 12.5, 7.3)

# %% [markdown]
# ## Fig. S17: DD-associated genes contributing most to transmission excess
#
# Every DD-associated gene whose own excess (transmitted unrescued PTVs minus the number
# expected at the gene set's synonymous transmission rate) is at least two PTVs, labelled
# with its transmitted and total PTVs. Colour gives the published evidence of PTVs
# transmitted by unaffected or subclinically affected parents, curated from the literature
# (the input file `gene_penetrance_consolidated.tsv`); a gene absent from that file is drawn
# as "No verified evidence identified".

# %%
## RANKED BY THE GENE'S OWN EXCESS AGAINST A FIXED BASELINE. Both leave-one-out
## quantities are contaminated for this purpose: dropping a gene removes its
## synonymous variants too, and that baseline is a noisy estimate around 0.5,
## which is enough to reorder genes. Holding it at the full-set synonymous rate
## gives excess = n_transmitted - n_ptv * p_syn, deterministic in the gene's own
## counts.
load_penetrance <- function() {
  f <- PATHS$penetrance
  if (!file.exists(f)) {
    msg("!! %s not found; %s is drawn in one colour", basename(f), fig_label("figS17"))
    return(NULL)
  }
  p <- fread(f)
  if (!all(c("gene", "classification") %in% names(p))) {
    msg("!! the penetrance table lacks gene / classification")
    return(NULL)
  }
  p[, tier := fcase(grepl("documented", tolower(classification)), "documented",
                    grepl("uncertain", tolower(classification)), "uncertain",
                    grepl("no verified|no published|absence", tolower(classification)),
                    "none", default = NA_character_)]
  p[]
}
pen <- load_penetrance()

infl <- jack[gene_set == GS_DD]
if (nrow(infl)) {
  Np <- sum(infl$n_ptv); Tp <- sum(infl$n_transmitted)
  p_syn_dd <- (Tp - infl$excess_full[1]) / Np
  infl <- copy(infl)[, excess_fixed := n_transmitted - n_ptv * p_syn_dd]
  ## the secondary keys only order the bars; the threshold decides membership
  setorder(infl, -excess_fixed, -n_ptv, gene)
  infl <- infl[excess_fixed >= PARAMS$s17_min_excess - 1e-9]
  msg("Fig. S17: %d DD-associated genes with an excess of at least %g transmitted PTVs: %s",
      nrow(infl), PARAMS$s17_min_excess, paste(infl$gene, collapse = ", "))
  if (!nrow(infl)) msg("!! %s is not drawn", fig_label("figS17"))
  if (!is.null(pen) && length(setdiff(infl$gene, pen$gene)))
    msg("!! not in the curated penetrance table, drawn as 'no verified evidence': %s",
        paste(setdiff(infl$gene, pen$gene), collapse = ", "))
}

# %%
if (nrow(infl)) {
  infl[, tier := if (is.null(pen)) "none"
                 else { t <- pen$tier[match(as.character(gene), pen$gene)]
                        ## a gene absent from the curated table has not been
                        ## searched, which is not the same as searched and found
                        ## wanting; it is drawn in the "no evidence" colour
                        fifelse(is.na(t), "none", t) }]
  infl[, tier_lab := factor(TIER_LABELS[tier], levels = TIER_LABELS)]
  infl[, gene_lab := if (HAVE_GGTEXT) paste0("*", gene, "*") else as.character(gene)]
  infl[, gene_lab := factor(gene_lab, levels = rev(gene_lab))]
  save_figure("figS17", ggplot(infl, aes(excess_fixed, gene_lab, fill = tier_lab)) +
                geom_col(width = 0.72) +
                geom_text(aes(label = sprintf("%d/%d transmitted", n_transmitted, n_ptv)),
                          hjust = -0.08, size = FONT_ANNOT_SMALL, colour = "grey30") +
                ## only the categories drawn: a key with no bar reads as a gene missing
                scale_fill_manual(values = setNames(TIER_COLS[names(TIER_LABELS)],
                                                    TIER_LABELS), name = NULL,
                                  drop = TRUE) +
                scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
                labs(x = "Excess transmitted variants (observed minus the synonymous expectation)",
                     y = NULL) +
                theme_pub(grid = "x") +
                theme(axis.text.y = element_md(face = if (HAVE_GGTEXT) "plain"
                                                      else "italic")),
              ## the same bar thickness whatever the number of genes
              8.5, max(3, 1.2 + 0.28 * nrow(infl)))
}

# %%
sessionInfo()

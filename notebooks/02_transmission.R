# %% [markdown]
# # 02 · Inherited protein-truncating variants contribute to DD risk
#
# Tests whether inherited PTVs are transmitted from parents to affected children more often than rare synonymous variants, and what that transmission excess implies: how many PTVs are rescued in population controls, the fraction of carrier probands whose DD is attributable to an inherited PTV, and how inherited PTV burden relates to diagnosis. DDD and RUMC families are analysed together unless stated.
#
# Figures drawn:
#
# - **Figure 2**: (A) transmission by LOEUF decile; (B, C) rescued PTVs and corrected relative risk in population controls (Samocha *et al.*); (D) attributable fraction among carrier probands; (E) sHet burden of inherited PTVs and diagnosis
# - **Fig. S2**: transmission model by study
# - **Fig. S3**: rescued PTVs among untransmitted parental alleles, by gnomAD allele frequency cap
# - **Fig. S4**: attributable fraction among carrier probands, with unrescued PTVs in NMD target regions
# - **Fig. S5**: clinical interpretation of inherited PTVs in DD-associated genes
# - **Fig. S6**: parent of origin of clinically interpreted PTVs
# - **Fig. S7**: parent of origin of PTVs that were not clinically reported
# - **Fig. S8**: transmission into undiagnosed versus *de novo* diagnosed probands
# - **Fig. S9**: number of inherited PTVs in constrained genes and diagnosis
# - **Fig. S10**: families segregating two PTVs in DD-associated genes

# %%
REPO_DIR <- Sys.getenv("RR_REPO", unset = if (file.exists("config.R")) "." else "..")
source(file.path(REPO_DIR, "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
  library(grid); library(readxl); library(svglite)
})

# %% [markdown]
# ## Inputs
#
# The released parental and proband variant tables, and the per-person files for the analyses that pool variants within a proband or family (Figures 2D and 2E, Figs. S4, S5, S9 and S10). Call QC, the exclusion of hypermutated trios, the heterozygous-in-proband check on transmitted variants, gene sets, gene constraint, and the child's clinical status and genetic ancestry were applied or attached before release (`data/README.md`). Gene-level annotation and Samocha *et al.*'s supplementary tables are public resources (`resources/README.md`).

# %%
PATHS <- list(
  parental       = file.path(DATA_DIR, "parental_variants.tsv.gz"),
  probands       = file.path(DATA_DIR, "proband_variants.tsv.gz"),
  burden         = file.path(DATA_DIR, "burden_per_proband.tsv"),
  af_carriers    = file.path(DATA_DIR, "attributable_fraction_carriers.tsv"),
  cohort_trios   = file.path(DATA_DIR, "cohort_trios.tsv"),
  interpretation = file.path(DATA_DIR, "interpretation_counts.tsv"),
  two_ptv        = file.path(DATA_DIR, "two_ptv_families.tsv"),
  genes         = file.path(RESOURCE_DIR, "gene_annotation.tsv.gz"),
  ## Samocha et al. 2024, medRxiv supplementary tables (media-2.xlsx, media-3.txt)
  samocha_xlsx  = file.path(RESOURCE_DIR, "samocha_2024_supp_tables.xlsx"),
  samocha_st10  = file.path(RESOURCE_DIR, "samocha_2024_supp_table10.txt")
)

# %% [markdown]
# ## Parameters and plotting helpers
#
# Shared constants, the plotting theme, and small helpers: `tidy_glm()` returns odds ratios with Wald 95% confidence intervals, and `add_sig()` adds Bonferroni-adjusted significance stars. In every transmission model the intercept is the reference stratum (rare synonymous variants in LOEUF deciles 3-9), so it carries no star.

# %%
STUDIES      <- c("DDD", "RUMC")
REF_LABEL    <- "3 - 9\nSynonymous"
BIN_LEVELS   <- c(GS_DD, "0", "1", "2", "3 - 9")
CLASS_LEVELS <- c("Unrescued PTVs", "Rescued PTVs", "Synonymous")
CONSTRAINED  <- c(GS_DD, "0", "1", "2")     # LOEUF < 0.58
N_MC   <- 2e5    # Monte Carlo draws, attributable fraction and Samocha correction
N_BOOT <- 1000   # bootstrap replicates, Figure 2E and Fig. S9
set.seed(1)      # the analyses that draw random numbers set their own seeds as well

## Samocha et al. case sets (Figures 2B and 2C). PRIMARY: cases with unaffected
## parents and no de novo diagnosis, known de novo variants removed. The de
## novo-inclusive set is carried for comparison only: its OR mixes de novo and
## inherited PTVs, while the rescued fractions describe inherited alleles.
SAMOCHA_SETS <- data.table(
  sheet    = c("SuppTable4", "SuppTable4", "SuppTable2"),
  case_set = c("unaffpar_noDNdiag_noDNM", "affpar_noDNM", "unaffpar_noDNdiag"),
  label    = c("Unaffected parents, de novos removed (primary)",
               "Affected parents, de novos removed",
               "Unaffected parents, de novos included (not valid)"),
  valid    = c(TRUE, TRUE, FALSE))
SAMOCHA_CONSEQUENCE <- "lof_HC"
## Samocha et al. define rare variants as gnomAD AF < 1e-5 (plus a within-
## dataset AF < 1e-4 in her joint DDD + INTERVAL call set, which we cannot
## reproduce and do not approximate). For THIS ANALYSIS ONLY, the rescued
## fraction in our probands uses her gnomAD threshold, so it is computed on the
## variants her odds ratios describe; every other analysis keeps MAF < 0.001.
## Variants absent from gnomAD (missing AF) are rare and kept.
SAMOCHA_MAX_GNOMAD_AF <- 1e-5
## gnomAD AF caps for the rescued share among untransmitted parental alleles
## (Fig. S3), from Samocha et al.'s threshold up to 0.1%
AF_CAPS <- c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3)
## Samocha et al.'s constrained gene count (Supplemental Note), for the log
SAMOCHA_N_GENES <- c(cons = 2699L, cons_dd = 422L, cons_notdd = 2277L)
## Gene sets: ours <-> Samocha categories.
GS_ALL  <- "All constrained genes"
SAM_MAP <- data.table(stratum  = c(GS_DD, GS_HIGH, GS_ALL),
                      category = c("cons_dd", "cons_notdd", "cons"),
                      label    = c("DD-associated", "Other\nconstrained",
                                   "All\nconstrained"))

## Significance annotation on the transmission forests: stars come from
## Bonferroni-adjusted p-values (add_sig).
SIG_SHOW_NS <- FALSE
SIG_OFFSET  <- 0.07
SIG_CEILING <- 1.97
BAR_WIDTH_E <- 0.90   # clinical interpretation stack (Fig. S5)
BAR_WIDTH_H <- 0.75
HATCH_SLOPE_FACTOR <- 0.45
HATCH_SPACING_X    <- 0.05
HATCH_COL <- "grey30"; HATCH_LWD <- 0.35; HATCH_FILL <- "grey92"

# %%
FONT_BASE <- 15; FONT_LEGEND <- 14; FONT_TITLE <- 15; FONT_SUBTITLE <- 12
FONT_TAG <- 20; ANNOT <- 4.6; ANNOT_SMALL <- 3.8

theme_ms <- function(base_size = FONT_BASE, legend_size = FONT_LEGEND) {
  theme(text = element_text(family = "sans", size = base_size, color = "black"),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1.2),
        axis.line.x.top = element_blank(), axis.line.y.right = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        axis.ticks.length = unit(3, "pt"),
        axis.text = element_text(size = base_size),
        legend.position = "bottom", legend.box = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size = legend_size), legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = FONT_TITLE, face = "bold"))
}

CLASS_COLS <- c("Unrescued PTVs" = "tomato4", "Rescued PTVs" = "darkgoldenrod",
                "Synonymous" = "grey50")
SYM_SHAPES <- c("All" = 19, "Clinically interpreted variants removed" = 17)
## whole numbers with thousands separators, for the log
n_fmt <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")

# %%
## `is_reference` marks the INTERCEPT. Every other row is a contrast against it,
## so the intercept's p-value tests a different thing and must not carry a star
## alongside the contrasts.
tidy_glm <- function(fit, intercept_name = NULL) {
  s  <- summary(fit)$coefficients
  ci <- suppressMessages(confint.default(fit))
  out <- data.table(term = rownames(s), estimate = exp(s[, "Estimate"]),
                    std.error = s[, "Std. Error"], statistic = s[, "z value"],
                    p.value = s[, "Pr(>|z|)"], conf.low = exp(ci[, 1]),
                    conf.high = exp(ci[, 2]))
  out[, is_reference := term == "(Intercept)"]
  if (!is.null(intercept_name)) out[term == "(Intercept)", term := intercept_name]
  out[, term := sub("^label", "", term)]
  out[]
}

split_bin_type <- function(dt) {
  dt[, c("bin", "type") := tstrsplit(term, "\n", fixed = TRUE)]
  dt[, bin := factor(bin, levels = BIN_LEVELS)]
  dt[]
}

star <- function(p) fcase(p <= 0.001, "***", p <= 0.01, "**", p <= 0.05, "*",
                          p <= 0.10, ".", default = "ns")

CI_SEP <- " - "
fmt_or <- function(est, lo, hi, d = 3)
  sprintf(paste0("%.", d, "f (%.", d, "f", CI_SEP, "%.", d, "f)"), est, lo, hi)
odds_to_rate <- function(o) o / (1 + o)

## `by` sets the correction family. NULL adjusts across the whole table, right
## for a single fit; by = "study" keeps independent cohort fits apart.
add_sig <- function(dt, by = NULL) {
  if (!"is_reference" %in% names(dt)) dt[, is_reference := FALSE]
  dt[is.na(is_reference), is_reference := FALSE]
  if (is.null(by)) dt[, p.adjusted := p.adjust(p.value, "bonferroni")]
  else             dt[, p.adjusted := p.adjust(p.value, "bonferroni"), by = by]
  dt[, sig.level := star(p.adjusted)]
  dt[, sig.level := fifelse(is_reference, "reference", sig.level)]
  lab <- as.character(dt$sig.level)
  lab[dt$is_reference] <- ""
  if (!SIG_SHOW_NS) lab[lab %in% c("ns", ".")] <- ""
  dt[, sig.label := lab]
  dt[, ref.label := fifelse(is_reference, "ref", "")]
  dt[, sig.y := pmin(conf.high + SIG_OFFSET, SIG_CEILING)]
  dt[]
}

sig_layer <- function(pos, size = 5.5, ref_size = 4.0) list(
  geom_text(aes(y = sig.y, label = sig.label), position = pos, size = size,
            show.legend = FALSE, vjust = 0.35),
  geom_text(aes(y = sig.y, label = ref.label), position = pos, size = ref_size,
            fontface = "italic", show.legend = FALSE, vjust = 0.25))

# %% [markdown]
# ## Gene annotation and LOEUF deciles
#
# Every variant record carries its gene's set (DD-associated, pLI ≥ 0.9 or pLI < 0.9) and gnomAD v2.1.1 LOEUF decile, matched by Ensembl gene ID; `gene_annotation.tsv.gz` gives the same for every gene (one row per row of the gnomAD constraint table), for the gene lists of Figures 2B-D and Fig. S3 and the decile ranges on the axes. For the transmission models, genes outside the DD-associated set are binned by LOEUF decile, with deciles 3-9 collapsed into one bin.

# %%
GENES <- fread(PATHS$genes, encoding = "UTF-8")
## a gene absent from the annotation is in the pLI < 0.9 set
gene_set_of <- function(symbol) {
  gs <- GENES$gene_set[match(symbol, GENES$gene)]
  fifelse(is.na(gs), GS_LOW, gs)
}

## axis labels: the LOEUF range of each decile over the gnomAD constraint table
loeuf <- GENES[!is.na(LOEUF_decile) & !is.na(LOEUF)]
decile_labels <- loeuf[, .(min_loeuf = min(LOEUF), max_loeuf = max(LOEUF)), keyby = LOEUF_decile]
decile_labels_vec <- setNames(sprintf("%s (%s - %s)", decile_labels$LOEUF_decile,
                                      round(decile_labels$min_loeuf, 2),
                                      round(decile_labels$max_loeuf, 2)),
                              decile_labels$LOEUF_decile)
decile_labels_vec["3 - 9"] <- loeuf[LOEUF_decile %in% 3:9,
                                    sprintf("3 - 9 (%s - %s)", round(min(LOEUF), 2),
                                            round(max(LOEUF), 2))]
decile_labels_vec[GS_DD] <- GS_DD

# %% [markdown]
# ## Parental variant table
#
# One row per parental variant: unrescued PTVs (LOFTEE HC, not predicted to be rescued), rescued PTVs (predicted to be rescued in haplotype context) and rare synonymous variants. The transmission analyses use the records with `strict_qc` (outside ENCODE blacklist regions, and for DDD depth > 10 in both parents). Transcript rows are collapsed to one record per variant, parent, gene and class, and variants carried by both parents (one `obs_id`, two parents) are removed because their parent of origin is undefined. Two further record sets are built from the same tables: parental LOFTEE HC PTVs with a rescue prediction (the Figure 1B definition), for Figures 2B-D and Fig. S3, and the probands' own inherited PTVs at Samocha *et al.*'s gnomAD threshold, for Figure 2B.

# %%
par <- fread(PATHS$parental, colClasses = list(character = "CHROM"), encoding = "UTF-8")
## RescueRanger's call on each transcript row (1 rescued, 0 not rescued, NA no
## prediction); a variant's records are distinct by it
par[, rescue_call := fcase(rescue == "not_rescued", 0L,
                           rescue %in% c("rescued", "not_rescued_new_ptv"), 1L,
                           default = NA_integer_)]

classify_parents <- function(dt) rbind(
  dt[is_ptv == TRUE & rescue %in% "rescued"][, class := "Rescued PTVs"],
  dt[is_ptv == TRUE & rescue %in% "not_rescued" & vepLoF %in% "HC"][, class := "Unrescued PTVs"],
  dt[is_ptv == FALSE][, class := "Synonymous"], fill = TRUE)

## The gnomAD gene of each row: its row of the gnomAD constraint table, matched
## by Ensembl gene ID (a few symbols are on more than one row), or by symbol
## where the ID is not in that table and the symbol is on a single row of it. A
## row whose gene matches no gnomAD row joins the matched rows of the same
## variant and symbol, so a variant is one record per gnomAD gene.
CONSTRAINT <- GENES[!is.na(gene_id)]
strip_gene_version <- function(x) sub("\\..*$", "", as.character(x))
constraint_index <- function(symbol, gene_id) {
  i <- match(strip_gene_version(gene_id), strip_gene_version(CONSTRAINT$gene_id))
  j <- match(symbol, CONSTRAINT$gene)
  j[!symbol %in% CONSTRAINT[, .N, by = gene][N == 1L, gene]] <- NA_integer_
  fifelse(is.na(i), j, i)
}
assign_gnomad_gene <- function(d, by) {
  d[, gnomad_row := constraint_index(vepSYMBOL, vepGene)]
  d[, gnomad_row := if (anyNA(gnomad_row) && !all(is.na(gnomad_row)))
    fifelse(is.na(gnomad_row), gnomad_row[!is.na(gnomad_row)][1L], gnomad_row) else gnomad_row,
    by = by]
  d[]
}

keep_cols <- c("CHROM", "POS", "REF", "ALT", "obs_id", "vepSYMBOL", "gnomad_row", "Relationship",
               "Transmitted", "class", "study", "rescue_call", "gene_set",
               "clinical_interpretation", "proband_diagnosis_status",
               "inherited_ptv_diagnosis")
otl <- assign_gnomad_gene(classify_parents(par[strict_qc == TRUE]),
                          c("CHROM", "POS", "REF", "ALT", "obs_id", "vepSYMBOL"))
otl <- unique(otl[, ..keep_cols])

## drop variants carried by both parents (parent of origin undefined); a
## variant in two gnomAD genes that share a symbol is a record in each
otl[, parent_count := .N, by = .(CHROM, POS, REF, ALT, obs_id, vepSYMBOL, gnomad_row)]
otl <- otl[parent_count == 1][, parent_count := NULL]
otl[, LOEUF_decile := CONSTRAINT$LOEUF_decile[gnomad_row]]

otl[, leouf_decile := fifelse(gene_set == GS_DD, GS_DD,
                              fifelse(LOEUF_decile %in% 3:9, "3 - 9",
                                      as.character(LOEUF_decile)))]
otl <- otl[!is.na(leouf_decile)]
otl[, label := factor(paste0(leouf_decile, "\n", class))]
otl[, label := relevel(label, ref = REF_LABEL)]
otl[, class := factor(class, levels = CLASS_LEVELS)]
overtransmission_table_leouf <- otl
rm(otl)

cat(sprintf(paste0("Transmission table: %s records (%s transmitted unrescued PTVs in ",
                   "DD-associated genes).\n"),
            n_fmt(nrow(overtransmission_table_leouf)),
            n_fmt(overtransmission_table_leouf[class == "Unrescued PTVs" & leouf_decile == GS_DD &
                                                 Transmitted == TRUE, .N])))

# %%
## ---- Samocha comparison: shared definitions ---------------------------------
## PTV records for Figures 2B and 2C use the Figure 1B definition (LOFTEE HC,
## four PTV consequences, rescue prediction present) plus Samocha's gnomAD
## threshold, in probands AND in parents. The same parental records at every
## frequency (sam_par_all) give the All LOFTEE HC PTVs estimates of Figure 2D
## and the gnomAD AF caps of Fig. S3.
SAM_PTV_CONSEQUENCES <- c("frameshift_variant", "splice_acceptor_variant",
                          "splice_donor_variant", "stop_gained")
## gnomAD AF as a number: missing, "" or "." means absent from gnomAD (AF 0);
## several values joined by "&" or "," (one per overlapping transcript or
## allele) take the maximum, the conservative choice for a "rare" filter.
parse_gnomad_af <- function(x) {
  if (is.numeric(x)) return(fifelse(is.na(x), 0, x))
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  multi <- is.na(out) & !is.na(x) & grepl("[&,]", x)
  if (any(multi)) out[multi] <- vapply(strsplit(x[multi], "[&,]"), function(v) {
    v <- suppressWarnings(as.numeric(v)); if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  }, numeric(1))
  fifelse(is.na(out), 0, out)
}

## PARENTAL PTV records for Figure 2B: transmitted and untransmitted alleles
## under the Figure 1B definition, with the transmission table's QC, each row
## with its gnomAD AF. Variants carried by both parents are dropped, as for the
## transmission table:
##   sam_par      Samocha's gnomAD threshold (Figures 2B and 2C): the cap is
##                applied first
##   sam_par_all  every frequency, with gnomad_af (Figure 2D, Fig. S3)
one_parent <- function(d) {
  d[, n_par := uniqueN(Relationship), by = .(CHROM, POS, REF, ALT, obs_id)]
  d[n_par == 1][, n_par := NULL][]
}
sam_par_all <- unique(par[strict_qc == TRUE & vepLoF %in% "HC" &
                            vepConsequence %in% SAM_PTV_CONSEQUENCES & !is.na(rescue),
                          .(CHROM, POS, REF, ALT, obs_id, vepSYMBOL, gene_set, Relationship,
                            Transmitted, study, proband_superpopulation_group,
                            rescued = rescue == "rescued",
                            gnomad_af = parse_gnomad_af(GNOMAD_AF))])
for (st in STUDIES)
  message(sprintf("%s: %d of %d parental PTV records at gnomAD AF < %g (Samocha comparison).",
                  st, sam_par_all[study == st & gnomad_af < SAMOCHA_MAX_GNOMAD_AF, .N],
                  sam_par_all[study == st, .N], SAMOCHA_MAX_GNOMAD_AF))
sam_par <- one_parent(unique(sam_par_all[gnomad_af < SAMOCHA_MAX_GNOMAD_AF][
  , gnomad_af := NULL]))
sam_par_all <- one_parent(sam_par_all)

# %%
## Inherited PTVs carried by PROBANDS, for the Samocha comparison (Figure 2B),
## under the Figure 1B definition: the rescued fraction must match the one
## reported for Figure 1B.
prob <- fread(PATHS$probands, colClasses = list(character = "CHROM"), encoding = "UTF-8")
samocha_inh <- prob[is_denovo == FALSE & vepLoF %in% "HC" &
                      vepConsequence %in% SAM_PTV_CONSEQUENCES & !is.na(rescue)]
rm(prob)
for (st in STUDIES)
  message(sprintf("%s: %d of %d inherited PTVs kept at gnomAD AF < %g (Samocha comparison only).",
                  st, samocha_inh[study == st & parse_gnomad_af(GNOMAD_AF) < SAMOCHA_MAX_GNOMAD_AF, .N],
                  samocha_inh[study == st, .N], SAMOCHA_MAX_GNOMAD_AF))
samocha_inh <- samocha_inh[parse_gnomad_af(GNOMAD_AF) < SAMOCHA_MAX_GNOMAD_AF,
                           .(vepSYMBOL, gene_set, study, proband_superpopulation_group,
                             rescued = rescue == "rescued")]

# %% [markdown]
# ## The child's diagnostic status
#
# Every parental record carries two attributes of the child: `proband_diagnosis_status` (diagnosed: DDD diagnostic tiers 4-5 or a RUMC report with an outcome other than a VUS; VUS: a variant of uncertain significance and no confident diagnosis; undiagnosed; or not_available, without a clinical record) and `inherited_ptv_diagnosis` (the child carries a transmitted unrescued PTV in a DD-associated gene clinically interpreted pathogenic, likely pathogenic or uncertain). Together they give Fig. S8's three groups: inherited PTV diagnostic, *de novo* diagnosis (a confident diagnosis or a VUS elsewhere in the genome) and undiagnosed; probands without a clinical record are in none of them. Transmitted unrescued PTVs in DD-associated genes also carry the clinical interpretation reported at their position (`clinical_interpretation`; Fig. S6).

# %%
## Fig. S8's groups. A clinically interpreted inherited PTV is SUFFICIENT for
## 'inherited' and does not additionally require a diagnosis: the
## interpretation and the diagnosis record come from different sources in both
## cohorts.
overtransmission_table_leouf[, dx_status := fcase(
  proband_diagnosis_status == "not_available", NA_character_,
  inherited_ptv_diagnosis, "inherited",
  proband_diagnosis_status %in% c("diagnosed", "VUS"), "other",
  default = "none")]

# %% [markdown]
# ## Figure 2A and Fig. S2: transmission by LOEUF decile
#
# A single logistic regression of transmission status on the interaction of variant class and LOEUF bin, with DD-associated genes as their own bin and rare synonymous variants in deciles 3-9 as the reference level. Points are odds ratios with Wald 95% confidence intervals; stars are Bonferroni-adjusted. Triangles refit the model after removing every proband carrying a clinically interpreted inherited PTV in a DD-associated gene (`inherited_ptv_diagnosis`: an unrescued PTV interpreted pathogenic, likely pathogenic or uncertain). Fig. S2 fits the same model separately in each cohort, with Bonferroni correction within each.

# %%
fit_panelA <- function(dat) {
  m_all  <- glm(Transmitted ~ label, binomial("logit"), dat)
  d_all  <- split_bin_type(tidy_glm(m_all, REF_LABEL))[, sym := "All"]
  m_clin <- glm(Transmitted ~ label, binomial("logit"),
                dat[inherited_ptv_diagnosis == FALSE])
  d_clin <- split_bin_type(tidy_glm(m_clin, REF_LABEL))
  d_clin <- d_clin[bin == GS_DD & type == "Unrescued PTVs"]
  d_clin[, sym := "Clinically interpreted variants removed"]
  list(all = d_all, clin = d_clin)
}

fa <- fit_panelA(overtransmission_table_leouf)
model_df <- rbind(fa$all, fa$clin, fill = TRUE)
model_df[, type := factor(type, levels = CLASS_LEVELS)]
model_df <- model_df[!is.na(bin) & !is.na(type)]
add_sig(model_df)

Y_LAB_A <- "Estimated transmission odds"
panelA_df <- model_df[sym == "All" |
                        (bin == GS_DD & sym == "Clinically interpreted variants removed")]

## Dodge order within each bin, derived from the data rather than hard-coded.
sym_order <- panelA_df[bin == GS_DD & type == "Unrescued PTVs"][
  order(-estimate), as.character(sym)]
sym_order <- c(sym_order, setdiff(unique(as.character(panelA_df$sym)), sym_order))
series_levels <- unlist(lapply(CLASS_LEVELS, function(tp)
  paste(tp, sym_order, sep = " | ")))
panelA_df[, series := factor(paste(as.character(type), sym, sep = " | "),
                             levels = series_levels)]
panelA_df[, series := droplevels(series)]
stopifnot(!anyNA(panelA_df$series))

## `preserve = "single"` keeps point spacing equal across bins; without it the
## DD-associated bin, which carries a fourth series, packs more tightly.
pdA <- position_dodge(width = 0.8, preserve = "single")
y_lim_A <- c(0.5, max(2.0, max(panelA_df$conf.high, na.rm = TRUE)))

panel_A <- ggplot(panelA_df, aes(x = bin, y = estimate, ymin = conf.low,
                                 ymax = conf.high, colour = type,
                                 shape = factor(sym, levels = names(SYM_SHAPES)),
                                 group = series)) +
  coord_flip() +
  geom_hline(yintercept = 1, colour = "black", linewidth = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.8, position = pdA, show.legend = FALSE) +
  geom_point(size = 5, position = pdA) +
  sig_layer(pdA) +
  scale_colour_manual(values = CLASS_COLS, drop = FALSE) +
  scale_shape_manual(values = SYM_SHAPES,
                     breaks = c("All", "Clinically interpreted variants removed")) +
  scale_y_continuous(limits = y_lim_A, breaks = c(0.5, 1.0, 1.5, 2.0),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = decile_labels_vec) +
  xlab("LOEUF decile") + ylab(Y_LAB_A) +
  ## A's legend is its own (two rows: PTV class, then subset), placed under the
  ## plot by the figure grid (Figure 2 assembly, below).
  guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 19)),
         shape = guide_legend(order = 2, nrow = 1)) +
  theme_ms() + theme(legend.spacing.y = unit(0.15, "cm"))

# %%
base_row <- model_df[is_reference == TRUE & type == "Synonymous"]
cat(sprintf(paste0("\nEmpirical transmission baseline (rare synonymous, LOEUF 3-9): ",
                   "rate = %.4f [%.4f - %.4f] (odds %.4f)\n"),
            odds_to_rate(base_row$estimate), odds_to_rate(base_row$conf.low),
            odds_to_rate(base_row$conf.high), base_row$estimate))

# %%
by_study <- rbindlist(lapply(STUDIES, function(st) {
  f <- fit_panelA(overtransmission_table_leouf[study == st])
  out <- rbind(f$all, f$clin, fill = TRUE)
  out[, study := st][]
}), fill = TRUE)
by_study[, type := factor(type, levels = CLASS_LEVELS)]
by_study <- by_study[!is.na(bin) & !is.na(type)]
add_sig(by_study, by = "study")

pd2 <- position_dodge(width = 0.8)
supp2_df <- by_study[sym == "All" |
                       (bin == GS_DD & sym == "Clinically interpreted variants removed")]
supp2_df[, series := factor(paste(as.character(type), sym, sep = " | "))]

fig_s2 <- ggplot(supp2_df, aes(x = bin, y = estimate, ymin = conf.low,
                                 ymax = conf.high, colour = type,
                                 shape = factor(sym, levels = names(SYM_SHAPES)),
                                 group = series)) +
  geom_hline(yintercept = 1, colour = "black", linewidth = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.8, position = pd2) +
  geom_point(size = 5, position = pd2) + sig_layer(pd2) +
  coord_flip() + facet_wrap(. ~ study, nrow = 1) +
  scale_colour_manual(values = CLASS_COLS, drop = FALSE) +
  scale_shape_manual(values = SYM_SHAPES,
                     breaks = c("All", "Clinically interpreted variants removed")) +
  scale_y_continuous(limits = c(0.5, 2.0), breaks = c(0.5, 1.0, 1.5, 2.0),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = decile_labels_vec) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 1)) +
  xlab("LOEUF decile") + ylab("Estimated transmission odds (95% CI)") +
  theme_ms() + theme(legend.box = "vertical", legend.spacing.y = unit(0.15, "cm"))
save_figure("figS2", fig_s2, 14, 7.5)

# %% [markdown]
# ## Fig. S8: transmission into undiagnosed versus *de novo* diagnosed probands
#
# Transmission odds of unrescued PTVs in three gene bins (DD-associated; constrained, LOEUF deciles 0-2 without DD-associated genes; unconstrained, deciles 3-9), relative to rare synonymous variants in unconstrained genes of the same probands. Each proband group is its own fit: all probands (circles, as Figure 2A), probands with a *de novo* diagnosis (squares), probands without one (diamonds), and probands without one after removing those whose inherited PTV in a DD-associated gene was clinically interpreted (triangles). Grey circles are rare synonymous variants. Stars are Bonferroni-adjusted across the estimates drawn.

# %%
BIN3 <- c(GS_DD, "Constrained", "Unconstrained")
lo3 <- GENES[LOEUF_decile %in% 3:9 & !is.na(LOEUF), min(LOEUF)]
BIN3_LAB <- setNames(c(GS_DD, sprintf("Constrained\n(LOEUF < %.2f)", lo3),
                       sprintf("Unconstrained\n(LOEUF \u2265 %.2f)", lo3)), BIN3)
REF3 <- paste0(BIN3[3], "\nSynonymous")
## transmitted / untransmitted counts per proband group, study, bin and class
## (every count below is a sum of these, so the large table is read once)
cnt_grp <- overtransmission_table_leouf[, .(
  T = sum(Transmitted == TRUE), U = sum(Transmitted == FALSE)),
  by = .(dx_status, study = as.character(study),
         bin = fcase(leouf_decile == GS_DD, BIN3[1],
                     leouf_decile %in% c("0", "1", "2"), BIN3[2],
                     leouf_decile == "3 - 9", BIN3[3]),
         class = as.character(class))]
cnt_grp <- cnt_grp[!is.na(bin)]

DXG_RM_LAB <- "Clinically interpreted inherited PTVs removed"
## status NULL: all probands
DXG <- list(
  all     = list(status = NULL,                    removed = FALSE, dd_only = FALSE),
  dn      = list(status = "other",                 removed = FALSE, dd_only = FALSE),
  nodn    = list(status = c("none", "inherited"),  removed = FALSE, dd_only = FALSE),
  nodn_rm = list(status = "none",                  removed = TRUE,  dd_only = TRUE))
## the counts of one proband group, within one study or both
dxg_cnt <- function(status, st = "All") {
  d <- if (is.null(status)) cnt_grp else cnt_grp[dx_status %in% status]
  if (st != "All") d <- d[study == st]
  d
}
## one fit per proband group: the panel A model on per-label counts
dxg_fit <- function(key) {
  g <- DXG[[key]]
  d <- dxg_cnt(g$status)
  cnt <- d[, .(T = sum(T), U = sum(U)), by = .(bin, class)]
  cnt[, label := relevel(factor(paste0(bin, "\n", class)), ref = REF3)]
  out <- tidy_glm(glm(cbind(T, U) ~ label, binomial("logit"), cnt), REF3)
  out[, c("bin", "type") := tstrsplit(term, "\n", fixed = TRUE)]
  out <- merge(out, cnt[, .(bin, type = class, T, U)], by = c("bin", "type"))
  out[, `:=`(group = key, removed = g$removed)]
  keep <- if (key == "all") c("Unrescued PTVs", "Synonymous") else "Unrescued PTVs"
  out <- out[type %in% keep]
  if (g$dd_only) out <- out[bin == GS_DD]
  out[]
}
dxg <- rbindlist(lapply(names(DXG), dxg_fit), use.names = TRUE)
dxg[type == "Synonymous", group := "syn"]

## the all-proband estimates are Figure 2A's (its per-variant decile model)
dxg_2a <- model_df[sym == "All" & bin == GS_DD & type == "Unrescued PTVs", estimate]
dxg_here <- dxg[group == "all" & bin == GS_DD & type == "Unrescued PTVs", estimate]
if (length(dxg_2a) != 1L || length(dxg_here) != 1L || abs(dxg_here / dxg_2a - 1) > 1e-6)
  stop("Transmission by diagnosis: the all-proband DD-associated estimate (",
       signif(dxg_here, 6), ") is not panel A's (", signif(dxg_2a, 6), ").")

# %%
## ---- figure: Figure 2A's layers, theme and stars -----------------------------
DXG_SERIES <- c("all", "nodn", "dn", "nodn_rm", "syn")   # bottom to top in a bin
DXG_TYPES  <- c("Unrescued PTVs", "Synonymous")
DXG_SET_OF <- c(all = "all", dn = "dn", nodn = "nodn", nodn_rm = "nodn_rm", syn = "all")
## all probands: panel A's circle; interpreted PTVs removed: panel A's
## triangle; the two diagnosis groups: a diamond (23, filled in the series
## colour) and a square
DXG_SHAPES <- c(all = SYM_SHAPES[["All"]], nodn = 23, dn = 15,
                nodn_rm = SYM_SHAPES[["Clinically interpreted variants removed"]])
DXG_SET_LABS <- c(all = expression("All probands"),
                  nodn = expression("No " * italic("de novo") * " diagnosis"),
                  dn = expression(italic("De novo") * " diagnosis"),
                  nodn_rm = as.expression(bquote("No " * italic("de novo") * " diagnosis, " *
                                                   .(sub("^C", "c", DXG_RM_LAB)))))
pdx <- copy(dxg)
pdx[, `:=`(bin = factor(bin, levels = BIN3), series = factor(group, levels = DXG_SERIES),
           colour_key = factor(type, levels = DXG_TYPES),
           set_key = factor(DXG_SET_OF[group], levels = names(DXG_SHAPES)))]
dxg_y_hi <- max(2.0, max(pdx$conf.high, na.rm = TRUE) + SIG_OFFSET + 0.08)
dxg_y_lim <- c(min(0.5, min(pdx$conf.low, na.rm = TRUE)), dxg_y_hi)
pdx <- add_sig(pdx)
## stars beyond the widest interval: this figure's own ceiling (the global
## SIG_CEILING stays as the other forests use it)
pdx[, sig.y := pmin(conf.high + SIG_OFFSET, dxg_y_hi - 0.03)]
## Series positions within a bin, set explicitly rather than by position_dodge:
## the bins carry different numbers of series (five in DD-associated, three
## elsewhere), and position_dodge(preserve = "single") does not dodge points in
## every ggplot2 version. Same spacing in every bin, centred on the bin.
DXG_SLOT <- 0.15
pdx[, slot := frank(match(as.character(series), DXG_SERIES), ties.method = "dense"), by = bin]
pdx[, xpos := as.integer(bin) + (slot - (max(slot) + 1) / 2) * DXG_SLOT, by = bin]
p_dx <- ggplot(pdx, aes(x = xpos, y = estimate, ymin = conf.low, ymax = conf.high,
                        colour = colour_key, fill = colour_key, shape = set_key,
                        group = series)) +
  coord_flip() +
  geom_hline(yintercept = 1, colour = "black", linewidth = 0.6) +
  geom_errorbar(width = 0.05, linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 5) +
  sig_layer(position_identity()) +
  scale_colour_manual(values = CLASS_COLS[DXG_TYPES], breaks = DXG_TYPES, drop = FALSE) +
  scale_fill_manual(values = CLASS_COLS[DXG_TYPES], guide = "none") +
  scale_shape_manual(values = DXG_SHAPES, breaks = names(DXG_SHAPES), labels = DXG_SET_LABS,
                     drop = FALSE) +
  scale_y_continuous(limits = dxg_y_lim,
                     breaks = seq(0.5, ceiling(dxg_y_hi * 2) / 2, by = 0.5),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = seq_along(BIN3), labels = unname(BIN3_LAB[BIN3]),
                     limits = c(0.5, length(BIN3) + 0.5), expand = c(0, 0)) +
  xlab(NULL) + ylab("Estimated transmission odds (95% CI)") +
  guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 19)),
         shape = guide_legend(order = 2, nrow = 2,
                              override.aes = list(colour = "black", fill = "black"))) +
  theme_ms() +
  theme(legend.spacing.y = unit(0.15, "cm"), panel.grid.minor = element_blank(),
        legend.text = element_text(size = FONT_LEGEND, hjust = 0))
save_figure("figS8", p_dx, 11, 8)

# %%
grp_records <- function(status) n_fmt(dxg_cnt(status)[, sum(T + U)])
cat(sprintf(paste0("\nTransmission by de novo diagnosis, records (transmitted and untransmitted): ",
                   "all probands %s; de novo diagnosis %s; no de novo diagnosis %s (%s of probands ",
                   "with a clinically interpreted inherited PTV, removed for the triangle); ",
                   "no clinical record %s.\n"),
            grp_records(NULL), grp_records("other"), grp_records(c("none", "inherited")),
            grp_records("inherited"), n_fmt(cnt_grp[is.na(dx_status), sum(T + U)])))

# %% [markdown]
# ## Fig. S5: clinical interpretation of inherited PTVs
#
# Probands carrying an inherited PTV in a DD-associated gene, classified by the most severe clinical interpretation among their PTVs of that class (pathogenic or likely pathogenic, uncertain, benign or likely benign). Where no PTV was reported, probands are split by whether any diagnostic finding was reported for them; this inclusive definition counts a VUS as a finding. Pooled and by study; probands without a clinical record are left out.
#
# > **Per-individual data:** each proband is counted once, by the most severe interpretation among their PTVs: read from `interpretation_counts.tsv` (probands by study, PTV class and category).

# %%
E_LEVELS <- c("Pathogenic/likely pathogenic", "Uncertain", "Benign/likely benign",
              "Not reported: diagnosed", "Not reported: undiagnosed")
E_COLS <- c("Pathogenic/likely pathogenic" = "tomato4", "Uncertain" = "darkorange3",
            "Benign/likely benign" = "skyblue3",
            "Not reported: diagnosed" = HATCH_FILL,
            "Not reported: undiagnosed" = "grey72")

interp_counts <- fread(PATHS$interpretation)
E_dat <- rbind(interp_counts[, .(n = sum(n_probands)), by = .(class, cat = category)][
                 , study := "All"],
               interp_counts[, .(study, class, cat = category, n = n_probands)])
E_dat[, cat := factor(cat, levels = E_LEVELS)]
E_dat[, class := factor(class, levels = c("Unrescued PTVs", "Rescued PTVs"))]
E_dat[, study := factor(study, levels = c("All", STUDIES))]
setorder(E_dat, study, class, cat)
cat(sprintf(paste0("\nProbands with a diagnosis record carrying a transmitted unrescued PTV in a ",
                   "DD-associated gene: %s.\n"),
            n_fmt(interp_counts[class == "Unrescued PTVs", sum(n_probands)])))

# %%
## Diagonal hatching for "Not reported: diagnosed", drawn with geom_segment
## clipped to the stacked rectangle: no ggpattern, no device pattern support.
stack_dt <- copy(E_dat)
setorder(stack_dt, study, class, -cat)
stack_dt[, ytop := cumsum(n), by = .(study, class)]
stack_dt[, ybot := ytop - n]
stack_dt[, xc := as.integer(class)]
stack_dt[, `:=`(x0 = xc - BAR_WIDTH_E / 2, x1 = xc + BAR_WIDTH_E / 2)]

hatch_segments <- function(rects, slope, spacing) {
  rbindlist(lapply(seq_len(nrow(rects)), function(i) {
    r <- rects[i]
    if (r$ytop <= r$ybot) return(NULL)
    a <- seq(r$x0 - r$ytop / slope, r$x1 - r$ybot / slope, by = spacing)
    xl <- pmax(r$x0, a + r$ybot / slope); xr <- pmin(r$x1, a + r$ytop / slope)
    keep <- xr > xl
    if (!any(keep)) return(NULL)
    data.table(study = r$study, x = xl[keep], xend = xr[keep],
               y = slope * (xl[keep] - a[keep]), yend = slope * (xr[keep] - a[keep]))
  }))
}
hatch_rects <- stack_dt[cat == "Not reported: diagnosed"]
ymax_by_study <- stack_dt[, .(ymax = max(ytop)), by = study]
hatch_dt <- rbindlist(lapply(levels(E_dat$study), function(s) {
  rr <- hatch_rects[study == s]
  if (nrow(rr) == 0) return(NULL)
  hatch_segments(rr, slope = ymax_by_study[study == s, ymax] * HATCH_SLOPE_FACTOR,
                 spacing = HATCH_SPACING_X)
}))
hatch_dt[, study := factor(study, levels = c("All", STUDIES))]

`%||%` <- function(a, b) if (is.null(a)) b else a
draw_key_maybe_hatch <- function(data, params, size) {
  base <- rectGrob(gp = gpar(col = data$colour %||% "black", fill = data$fill, lwd = 0.8))
  if (!identical(as.character(data$fill), HATCH_FILL)) return(base)
  cvals <- seq(-1, 1, by = 0.22); xs <- ys <- numeric(0); ids <- integer(0)
  for (k in seq_along(cvals)) {
    cc <- cvals[k]; x_lo <- max(0, cc); x_hi <- min(1, 1 + cc)
    if (x_hi <= x_lo) next
    xs <- c(xs, x_lo, x_hi); ys <- c(ys, x_lo - cc, x_hi - cc); ids <- c(ids, k, k)
  }
  grobTree(base, polylineGrob(xs, ys, id = ids,
                              gp = gpar(col = HATCH_COL, lwd = HATCH_LWD * 2)),
           rectGrob(gp = gpar(col = "black", fill = NA, lwd = 0.8)))
}

fig_s5 <- ggplot(E_dat, aes(class, n, fill = cat)) +
  geom_col(colour = "black", linewidth = 0.4, width = BAR_WIDTH_E,
           key_glyph = draw_key_maybe_hatch) +
  geom_segment(data = hatch_dt, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, colour = HATCH_COL, linewidth = HATCH_LWD) +
  facet_wrap(. ~ study, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = E_COLS, breaks = E_LEVELS, drop = FALSE) +
  scale_x_discrete(labels = c("Unrescued PTVs" = "Unrescued\nPTVs",
                              "Rescued PTVs" = "Rescued\nPTVs")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  guides(fill = guide_legend(nrow = 2)) +
  xlab("") + ylab("Number of probands") +
  theme_ms() + theme(panel.grid.major.x = element_blank(),
                     legend.position = "bottom",
                     legend.spacing.x = unit(0.4, "cm"),
                     legend.text = element_text(size = FONT_LEGEND))
save_figure("figS5", fig_s5, 10, 6)

# %% [markdown]
# ## Fig. S9: number of inherited PTVs and diagnosis
#
# Probands with a diagnostic label are analysed, leaving out probands whose most severe reported finding is a VUS and probands carrying a clinically reported inherited PTV in a DD-associated gene. Genetic ancestry is collapsed to European versus non-European. For each proband, distinct transmitted unrescued PTVs in constrained genes (LOEUF < 0.58, deciles 0-2 including DD-associated genes) are counted; the model is fitted among probands carrying at least one. Diagnostic status is modelled by logistic regression on that count, adjusted for study and ancestry; points are observed proportions with 95% Wilson intervals, and the line is the prediction standardised over the cohort's study and ancestry mix. The average marginal effect per additional PTV, in percentage points with a 95% interval from 1,000 bootstrap replicates (Methods), is printed after the model.
#
# > **Per-individual data:** the PTV count is a sum over each proband's variants, and the model needs one row per proband: read from `burden_per_proband.tsv` (study, ancestry group, diagnosis, sHet burden and number of PTVs in constrained genes).

# %%
## Each analysis that draws random numbers sets its own seed, so its intervals
## do not depend on what ran before it.
set.seed(12)
burden <- fread(PATHS$burden)
burden[, diagnosed := factor(diagnosed, levels = c("Diagnosed", "Undiagnosed"))]
cat(sprintf("\nFigure 2E cohort: %s probands, %s of them with an ancestry call (the models).\n",
            n_fmt(nrow(burden)), n_fmt(burden[!is.na(ancestry_group), .N])))

# %%
to_bin <- function(x) as.integer(as.character(x) == "Diagnosed")
STUDY_TERM <- if (length(STUDIES) > 1) "study" else character(0)
f_adj  <- function(x) reformulate(c(x, STUDY_TERM, "ancestry_group"), "y")
ame <- function(fit, term) {
  p <- predict(fit, type = "response")
  mean(p * (1 - p)) * coef(fit)[[term]] * 100
}
ame_boot <- function(formula, term, data, B = N_BOOT) {
  v <- replicate(B, {
    d <- data[sample(.N, .N, replace = TRUE)]
    f <- tryCatch(glm(formula, binomial, d), error = function(e) NULL)
    if (is.null(f)) NA_real_ else ame(f, term)
  })
  quantile(v, c(.025, .975), na.rm = TRUE)
}
std_curve <- function(fit, xname, xgrid, ref) {
  g <- ref[rep(seq_len(.N), length(xgrid))]
  set(g, j = xname, value = rep(xgrid, each = nrow(ref)))
  set(g, j = "p", value = predict(fit, newdata = g, type = "response"))
  g[, .(p = sum(p * w)), by = xname][order(get(xname))]
}
wilson <- function(k, n) {
  if (n == 0) return(c(NA, NA, NA))
  p <- k / n; z <- 1.96; d <- 1 + z^2 / n
  c(p, (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}

# %%
cnt <- burden[n_ptvs_constrained > 0, .(n_ptvs = n_ptvs_constrained, diagnosed, study,
                                        ancestry_group)]
cnt[, y := to_bin(diagnosed)]
cnt_a <- cnt[!is.na(ancestry_group)]
A_adj     <- glm(f_adj("n_ptvs"),  binomial, cnt_a)
A_ame     <- ame(A_adj, "n_ptvs")
A_ame_ci  <- ame_boot(f_adj("n_ptvs"), "n_ptvs", cnt_a)
cat(sprintf(paste0("\nAverage marginal effect per additional inherited PTV in a ",
                   "constrained gene: %.2f percentage points (95%% CI %.2f - %.2f).\n"),
            A_ame, A_ame_ci[1], A_ame_ci[2]))

# %%
ref_w_A <- cnt_a[, .(w = .N / nrow(cnt_a)), by = c("ancestry_group", STUDY_TERM)]
cnt_a[, kcap := pmin(n_ptvs, 6)]
obs_cnt <- cnt_a[, { w <- wilson(sum(y), .N); .(n = .N, prop = w[1], lo = w[2], hi = w[3]) },
                 by = kcap][order(kcap)]
obs_cnt[, lab := fifelse(kcap == 6, "6+", as.character(kcap))]
baseline_cnt <- mean(cnt_a$y)   # mean among EXPOSED probands, not the cohort mean
curve_cnt <- std_curve(A_adj, "n_ptvs", seq(1, 6, 0.05), ref_w_A)

p_burden_count <- ggplot(obs_cnt, aes(kcap, prop)) +
  geom_hline(yintercept = baseline_cnt, linetype = "dotted", colour = "grey40",
             linewidth = .7) +
  geom_line(data = curve_cnt, aes(n_ptvs, p), colour = "tomato4", linewidth = 1,
            alpha = .85) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "black", linewidth = .7) +
  geom_point(aes(size = n), colour = "tomato4") +
  geom_text(aes(y = hi + .035, label = paste0("n=", n)), size = ANNOT_SMALL,
            colour = "grey30") +
  scale_size_area(max_size = 8, guide = "none") +
  scale_x_continuous(breaks = obs_cnt$kcap, labels = obs_cnt$lab) +
  coord_cartesian(ylim = c(0, max(obs_cnt$hi) + .09)) +
  ylab("Proportion diagnosed") +
  xlab("Number of inherited PTVs in\nconstrained genes (LOEUF < 0.58)") +
  theme_ms() + theme(legend.position = "none", panel.grid.major.x = element_blank())
save_figure("figS9", p_burden_count, 8, 6.5)

# %% [markdown]
# ## Figure 2E: sHet burden of inherited PTVs and diagnosis
#
# The selection coefficient sHet is summed across the genes carrying a transmitted unrescued PTV in each proband (zero for probands with none), in Fig. S9's cohort. Diagnostic status is modelled by logistic regression on this burden, adjusted for study and ancestry. The line is the ancestry- and study-standardised prediction with a 1,000-replicate bootstrap band; points are observed proportions with 95% Wilson intervals in the zero-burden bin and in quintiles of positive burden.
#
# > **Per-individual data:** the sHet burden is a sum over each proband's variants, and the model needs one row per proband: read from `burden_per_proband.tsv`.

# %%
burden[, y := to_bin(diagnosed)]
ba <- burden[!is.na(ancestry_group)]
B_adj   <- glm(f_adj("shet_burden"),  binomial, ba)
ref_w  <- ba[, .(w = .N / nrow(ba)), by = c("ancestry_group", STUDY_TERM)]
grid_b <- seq(0, quantile(ba$shet_burden, 0.99), length.out = 120)
adj_curve <- std_curve(B_adj, "shet_burden", grid_b, ref_w)
bm <- matrix(NA_real_, N_BOOT, length(grid_b))
for (b in seq_len(N_BOOT)) {
  bs <- ba[sample(.N, .N, replace = TRUE)]
  f  <- tryCatch(glm(f_adj("shet_burden"), binomial, bs), error = function(e) NULL)
  if (!is.null(f))
    bm[b, ] <- tryCatch(std_curve(f, "shet_burden", grid_b, ref_w)$p,
                        error = function(e) rep(NA_real_, length(grid_b)))
}
adj_curve[, `:=`(lo = apply(bm, 2, quantile, .025, na.rm = TRUE),
                 hi = apply(bm, 2, quantile, .975, na.rm = TRUE))]
qs <- unique(quantile(ba[shet_burden > 0, shet_burden], seq(0, 1, .2)))
ba[, bin := fifelse(shet_burden == 0, "0",
                    as.character(cut(shet_burden, qs, include.lowest = TRUE)))]
obs_bin <- ba[, { w <- wilson(sum(y), .N)
                  .(x = mean(shet_burden), n = .N, prop = w[1], lo = w[2], hi = w[3]) },
              by = bin]

panel_E <- ggplot() +
  geom_ribbon(data = adj_curve, aes(shet_burden, ymin = lo, ymax = hi),
              fill = "tomato4", alpha = .20) +
  geom_line(data = adj_curve, aes(shet_burden, p), colour = "tomato4", linewidth = 1.2) +
  geom_pointrange(data = obs_bin, aes(x, prop, ymin = lo, ymax = hi),
                  colour = "black", linewidth = .6, fatten = 2.5) +
  ylab("Proportion diagnosed") + xlab("Total sHet burden of inherited PTVs") +
  theme_ms()

# %% [markdown]
# ## Figure 2B and 2C: rescue in population controls
#
# Samocha *et al.* compared rare LOFTEE HC PTVs between DD cases and INTERVAL controls, counting rescued PTVs as real. Rescued PTVs are transmitted at the neutral rate (Figure 2A), so they occur at the same per-person rate in cases and controls and only the unrescued count differs. The rescued fraction in controls is estimated from untransmitted parental alleles, scaled by their PTV burden relative to INTERVAL (E, from Samocha *et al.*'s per-gene counts, calibrated against synonymous variants and against unconstrained genes): r_ctrl = r_unt × E. The published odds ratio is then corrected to unrescued PTVs as OR_u = (OR − r_ctrl) / (1 − r_ctrl). Rescued fractions use LOFTEE HC PTVs with gnomAD AF < 1e-05, consistent with Samocha *et al.*, in their constrained genes, in all probands (DDD and RUMC). Intervals are Jeffreys intervals for observed fractions and 2 × 10⁵ Monte Carlo draws for estimates.

# %%
set.seed(13)   # Monte Carlo intervals of Figures 2B and 2C
stopifnot(file.exists(PATHS$samocha_xlsx))
read_samocha_sheet <- function(sheet) {
  d <- as.data.table(read_excel(PATHS$samocha_xlsx, sheet = sheet))
  setnames(d, tolower(names(d)))
  d[, sheet := sheet]
}
sam_pub <- rbindlist(lapply(unique(SAMOCHA_SETS$sheet), read_samocha_sheet), fill = TRUE)
sam_pub <- sam_pub[consequence == SAMOCHA_CONSEQUENCE & category %in% SAM_MAP$category]
sam_pub <- merge(sam_pub, SAMOCHA_SETS, by = c("sheet", "case_set"))
sam_pub <- sam_pub[, .(sheet, case_set, label, valid, category,
                       OR = as.numeric(odds_ratio), OR_lo = as.numeric(lower_bound),
                       OR_hi = as.numeric(upper_bound))]
stopifnot(nrow(sam_pub) == nrow(SAMOCHA_SETS) * nrow(SAM_MAP))

# %%
## Our rescued fractions: the Figure 1B definition plus Samocha's gnomAD
## threshold (parental variant table, above), restricted to pLI >= 0.9 genes
## (and to the genes of Samocha et al.'s Supplementary Table 10 when
## available). Gene sets are OURS, as in every other analysis, matched to her
## categories by meaning (SAM_MAP).
sam_universe <- unique(GENES[!is.na(pLI) & pLI >= 0.9, gene])
if (file.exists(PATHS$samocha_st10)) {
  sam_universe <- intersect(sam_universe, fread(PATHS$samocha_st10)$symbol)
  message("Rescue rates restricted to genes in Samocha et al.'s Supplementary Table 10.")
}
cat(sprintf("Samocha comparison, constrained genes used: %d.\n", length(sam_universe)))

sam_inh <- samocha_inh[vepSYMBOL %in% sam_universe]
sam_inh[, stratum := fifelse(gene_set == GS_DD, GS_DD, GS_HIGH)]

## PRIMARY cohort: all probands (DDD and RUMC). Rescue rates do not differ by
## ancestry, so European-ancestry DDD probands (the population Samocha's cases
## were drawn from) are kept as a sensitivity analysis only.
COHORT_PRIMARY <- "All probands (DDD and RUMC)"
COHORT_EUR     <- "DDD European probands"
is_ddd_eur <- function(d) d$study == "DDD" & d$proband_superpopulation_group == "European"
sam_rate_tab <- function(d, cohort) {
  d2 <- rbind(d, copy(d)[, stratum := GS_ALL])
  d2[, .(n_ptv = .N, n_rescued = sum(rescued)), by = stratum][, cohort := cohort]
}
sam_rates <- rbind(sam_rate_tab(sam_inh, COHORT_PRIMARY),
                   sam_rate_tab(sam_inh[is_ddd_eur(sam_inh)], COHORT_EUR))
sam_rates[, `:=`(rate = n_rescued / n_ptv,
                 lo = qbeta(.025, n_rescued + .5, n_ptv - n_rescued + .5),
                 hi = qbeta(.975, n_rescued + .5, n_ptv - n_rescued + .5))]
sam_rates <- merge(sam_rates, SAM_MAP, by = "stratum")
print(sam_rates[order(cohort, stratum)])

# %%
## Parental alleles (parental variant table), same genes, strata and cohorts as
## the probands. One record per variant per family.
sam_par13 <- unique(sam_par[vepSYMBOL %in% sam_universe],
                    by = c("CHROM", "POS", "REF", "ALT", "obs_id"))
sam_par13[, stratum := fifelse(gene_set == GS_DD, GS_DD, GS_HIGH)]
sam_jeff <- function(k, n, p) qbeta(p, k + .5, n - k + .5)
sam_par_tab <- function(d, cohort) {
  d2 <- rbind(d, copy(d)[, stratum := GS_ALL])
  d2[, .(n_tx = sum(Transmitted), k_tx = sum(rescued & Transmitted),
         n_unt = sum(!Transmitted), k_unt = sum(rescued & !Transmitted)),
     by = stratum][, cohort := cohort]
}
sam_par_rates <- rbind(sam_par_tab(sam_par13, COHORT_PRIMARY),
                       sam_par_tab(sam_par13[is_ddd_eur(sam_par13)], COHORT_EUR))
sam_par_rates[, `:=`(r_tx  = k_tx / n_tx,   r_tx_lo  = sam_jeff(k_tx, n_tx, .025),
                     r_tx_hi  = sam_jeff(k_tx, n_tx, .975),
                     r_unt = k_unt / n_unt, r_unt_lo = sam_jeff(k_unt, n_unt, .025),
                     r_unt_hi = sam_jeff(k_unt, n_unt, .975))]
print(sam_par_rates[order(cohort, stratum)])

# %%
## E: untransmitted PTVs in Samocha's trios relative to INTERVAL controls,
## synonymous-calibrated, per gene set (her Supplementary Table 10, same genes
## as above), then relative to unconstrained genes (main) or not (synonymous
## calibration only).
SAM_ENR_COLS <- c("inh.u.ptv", "inh.u.syn", "control.ptv", "control.syn")
sam_st10 <- if (file.exists(PATHS$samocha_st10)) fread(PATHS$samocha_st10) else NULL
if (!is.null(sam_st10) && all(c("symbol", SAM_ENR_COLS) %in% names(sam_st10))) {
  st <- sam_st10[symbol %in% sam_universe, c("symbol", SAM_ENR_COLS), with = FALSE]
  st[, stratum := fifelse(gene_set_of(symbol) == GS_DD, GS_DD, GS_HIGH)]
  st <- rbind(st, copy(st)[, stratum := GS_ALL])
  sam_enr <- st[, .(n_genes = .N, u_ptv = sum(inh.u.ptv), u_syn = sum(inh.u.syn),
                    c_ptv = sum(control.ptv), c_syn = sum(control.syn)), by = stratum]
  ## unconstrained baseline: her Supplementary Table 10 genes with pLI < 0.9
  base_genes <- intersect(unique(GENES[!is.na(pLI) & pLI < 0.9, gene]), sam_st10$symbol)
  b <- sam_st10[symbol %in% base_genes,
                .(u_ptv = sum(inh.u.ptv), u_syn = sum(inh.u.syn),
                  c_ptv = sum(control.ptv), c_syn = sum(control.syn))]
  E_base    <- (b$u_ptv / b$u_syn) / (b$c_ptv / b$c_syn)
  E_base_se <- sqrt(1 / b$u_ptv + 1 / b$u_syn + 1 / b$c_ptv + 1 / b$c_syn)
  cat(sprintf("Untransmitted vs INTERVAL in unconstrained genes (%d): %.3f (technical offset).\n",
              length(base_genes), E_base))
  sam_enr[, `:=`(E_unadj = (u_ptv / u_syn) / (c_ptv / c_syn),
                 E_unadj_se = sqrt(1 / u_ptv + 1 / u_syn + 1 / c_ptv + 1 / c_syn))]
  ## main: synonymous-calibrated and relative to unconstrained genes
  sam_enr[, `:=`(E = E_unadj / E_base, E_se = sqrt(E_unadj_se^2 + E_base_se^2),
                 E_base = E_base)]
  ## sensitivity: synonymous calibration only
  sam_enr[, `:=`(E_raw = E_unadj, E_raw_se = E_unadj_se)]
  sam_enr[, `:=`(E_lo = E * exp(-1.96 * E_se), E_hi = E * exp(1.96 * E_se))]
} else {
  warning("Samocha et al.'s Supplementary Table 10 (with inh.u.* and control.* ",
          "columns) not found: ",
          "E set to 1, i.e. r_ctrl = untransmitted rescued fraction.")
  sam_enr <- data.table(stratum = c(GS_DD, GS_HIGH, GS_ALL), n_genes = NA_integer_,
                        u_ptv = NA_real_, u_syn = NA_real_, c_ptv = NA_real_, c_syn = NA_real_,
                        E_raw = 1, E_raw_se = 0, E = 1, E_se = 0, E_base = 1, E_lo = 1, E_hi = 1)
}
print(sam_enr)

# %%
## Correction, with a Monte Carlo interval over the published OR (log-normal
## from its CI) and the parental rescued fractions (Jeffreys beta). Every case
## set, cohort and sensitivity value is drawn, in this order, so that the random
## draws behind Figures 2B and 2C are reproducible; the figure uses the primary
## case set and cohort.
sam_or_u <- function(OR, r_ctrl) ifelse(r_ctrl < 1, (OR - r_ctrl) / (1 - r_ctrl), NA_real_)
sam_q <- function(v, p) unname(quantile(v, p, na.rm = TRUE))
sam_grid <- merge(sam_pub, merge(merge(sam_rates, sam_par_rates, by = c("cohort", "stratum")),
                                  sam_enr[, .(stratum, E, E_se, E_lo, E_hi, E_raw, E_raw_se)],
                                  by = "stratum"),
                  by = "category", allow.cartesian = TRUE)
sam_corr <- sam_grid[, {
  or_d <- exp(rnorm(N_MC, log(OR), (log(OR_hi) - log(OR_lo)) / (2 * qnorm(.975))))
  ru_d <- rbeta(N_MC, k_unt + .5, n_unt - k_unt + .5)
  rt_d <- rbeta(N_MC, k_tx + .5, n_tx - k_tx + .5)
  e_d  <- exp(rnorm(N_MC, log(E), E_se))
  er_d <- exp(rnorm(N_MC, log(E_raw), E_raw_se))
  ## main: r_ctrl = r_unt * E
  r_ctrl <- r_unt * E;      rc_d <- ru_d * e_d
  OR_u   <- sam_or_u(OR, r_ctrl);   ou_d  <- sam_or_u(or_d, rc_d)
  ## sensitivity values
  OR_u0  <- sam_or_u(OR, r_unt);    ou0_d <- sam_or_u(or_d, ru_d)            # E = 1
  r_raw  <- r_unt * E_raw;  rr_d <- ru_d * er_d
  OR_ur  <- sam_or_u(OR, r_raw);    our_d <- sam_or_u(or_d, rr_d)            # synonymous calibration only
  r_pred <- r_tx * OR;      rp_d <- rt_d * or_d
  OR_up  <- sam_or_u(OR, r_pred);   oup_d <- sam_or_u(or_d, rp_d)            # r_tx x OR
  .(r_case = rate, r_case_lo = lo, r_case_hi = hi,
    n_tx, k_tx, r_tx, r_tx_lo, r_tx_hi, n_unt, k_unt, r_unt, r_unt_lo, r_unt_hi,
    E, E_lo, E_hi, E_se, E_raw,
    r_ctrl, r_ctrl_lo = sam_q(rc_d, .025), r_ctrl_hi = sam_q(rc_d, .975),
    p_ctrl_ge_1 = mean(rc_d >= 1),
    OR_u, OR_u_lo = sam_q(ou_d, .025), OR_u_hi = sam_q(ou_d, .975),
    OR_u_E1 = OR_u0, OR_u_E1_lo = sam_q(ou0_d, .025), OR_u_E1_hi = sam_q(ou0_d, .975),
    r_raw, OR_u_raw = OR_ur, OR_u_raw_lo = sam_q(our_d, .025), OR_u_raw_hi = sam_q(our_d, .975),
    r_pred, r_pred_lo = sam_q(rp_d, .025), r_pred_hi = sam_q(rp_d, .975),
    OR_u_pred = OR_up, OR_u_pred_lo = sam_q(oup_d, .025), OR_u_pred_hi = sam_q(oup_d, .975),
    or_ratio = OR_u / OR, excess_ratio = (OR_u - 1) / (OR - 1))
}, by = .(sheet, case_set, label = label.x, valid, category, stratum,
         strat_label = label.y, cohort, OR, OR_lo, OR_hi)]
setorder(sam_corr, cohort, sheet, case_set, category)
print(sam_corr[, .(cohort, case_set, category, r_tx, r_unt, E, r_ctrl, OR,
                   OR_u, OR_u_lo, OR_u_hi, OR_u_E1, OR_u_raw, OR_u_pred, valid)])
sam_main <- sam_corr[cohort == COHORT_PRIMARY & case_set == SAMOCHA_SETS$case_set[1]]

# %%
## --- Panels B and C -----------------------------------------------------------
## Panel B: rescued fraction among transmitted parental alleles (the probands'
## inherited PTVs; dark) and the estimate for INTERVAL controls (r_unt x E,
## i.e. r_ctrl; light), filled by gene set with the Figure 1B colours. Panel C:
## the published OR (grey) and the corrected OR for unrescued PTVs. Both panels
## are rotated (coord_flip) and stacked B above C in the right column of
## Figure 2, so the gene sets run down the side as in A and D, with
## DD-associated at the bottom. Intervals: Jeffreys for the observed fraction,
## Monte Carlo for the estimates.
SAM_A_LABELS <- c("DD cases (observed)", "INTERVAL controls (estimated)")
SAM_B_LABELS <- list(expression(paste("Samocha ", italic("et al."), " (all PTVs)")),
                     "Unrescued PTVs (corrected)")
SAM_OR_COLS  <- c(published = "grey50", corrected = CLASS_COLS[["Unrescued PTVs"]])
SAM_GS_FILL  <- GENE_SET_COLOURS[c(GS_DD, GS_HIGH)]
## Percentage labels sit centred INSIDE each bar. A bar shorter than this
## fraction cannot hold its label, which then goes just past the bar's end in
## dark text rather than straddling the edge (4% = 8% of the 0-50% axis).
## Every current bar (>= 14%) clears it; the rule is a guard for smaller values.
SAM_BAR_LABEL_MIN <- 0.04
SAM_B_MAX <- 0.5   # panel B x-axis upper limit

## all constrained genes combined are not drawn
sam_fig <- sam_main[category != "cons"]
sam_fig[, x := factor(strat_label, levels = SAM_MAP$label[1:2])]

## Two bars per gene set: transmitted parental alleles (= the probands'
## inherited PTVs) and INTERVAL controls (untransmitted x E). The untransmitted
## fraction behind the estimate is in the tables and notes; the transmission
## excess that separates the two is shown in panel A.
dB <- rbind(sam_fig[, .(x, stratum, grp = SAM_A_LABELS[1], y = r_tx,
                        lo = r_tx_lo, hi = r_tx_hi)],
            sam_fig[, .(x, stratum, grp = SAM_A_LABELS[2], y = r_ctrl,
                        lo = r_ctrl_lo, hi = pmin(r_ctrl_hi, 1))])
dB[, grp := factor(grp, levels = SAM_A_LABELS)]
dB[, inside := y >= SAM_BAR_LABEL_MIN]
dB[, `:=`(lab_y  = fifelse(inside, y / 2, y + 0.01),
          lab_hj = fifelse(inside, 0.5, 0),
          txt_col = fifelse(inside & grp == SAM_A_LABELS[1], "white", "grey15"))]
pd_B <- position_dodge(width = 0.75)
panel_B <- ggplot(dB, aes(x, y, group = grp)) +
  geom_col(aes(fill = stratum, alpha = grp), width = 0.7, position = pd_B,
           colour = "grey20", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.8,
                position = pd_B, na.rm = TRUE) +
  geom_text(aes(y = lab_y, label = percent(y, accuracy = 1), colour = txt_col,
                hjust = lab_hj),
            position = pd_B, size = ANNOT, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = SAM_GS_FILL, guide = "none") +
  scale_alpha_manual(values = c(1, 0.4)) +
  scale_colour_identity() +
  ## 0-50%: every bar and interval fits (largest upper limit ~42%); squish
  ## guards against an interval running past the edge.
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, SAM_B_MAX),
                     breaks = seq(0, SAM_B_MAX, 0.1), oob = scales::squish,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "PTVs rescued") +
  theme_ms() + theme(panel.grid.major.y = element_blank(),
                     plot.margin = margin(6, 18, 6, 6)) +
  guides(alpha = guide_legend(nrow = 1, override.aes = list(fill = "grey45")))

## Pairs dodged 0.5 of a category apart; the grey (published) label goes
## BELOW its point and the corrected label ABOVE, so each pair's numbers point
## away from each other and clear both intervals. The corrected OR uses
## r_ctrl = r_unt * E; its interval carries the uncertainty in her OR, r_unt
## and E.
dC <- rbind(sam_fig[, .(x, grp = "published", y = OR, lo = OR_lo, hi = OR_hi)],
            sam_fig[, .(x, grp = "corrected", y = OR_u, lo = OR_u_lo, hi = OR_u_hi)])
dC[, grp := factor(grp, levels = c("published", "corrected"))]
dC[, vj := fifelse(grp == "published", 1.9, -0.9)]
pd_C <- position_dodge(width = 0.5)
panel_C <- ggplot(dC, aes(x, y, colour = grp)) +
  geom_hline(yintercept = 1, linewidth = 0.6) +
  ## bar-ends as panel A's, ~0.06 in long at the Figure 2 layout:
  ## 0.107 / 2 dodged series = 0.053 of the 2.2-unit axis over a 2.5 in panel
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.107, linewidth = 0.8,
                position = pd_C, na.rm = TRUE) +
  geom_point(size = 5, position = pd_C, na.rm = TRUE) +
  geom_text(aes(label = formatC(y, format = "f", digits = 1), vjust = vj),
            position = pd_C, size = ANNOT, show.legend = FALSE, na.rm = TRUE) +
  coord_flip() +
  scale_colour_manual(values = unname(SAM_OR_COLS), labels = SAM_B_LABELS) +
  scale_y_log10(breaks = c(1, 2, 3, 4, 6, 8)) +
  labs(x = NULL,
       y = "Relative risk per inherited PTV\n(DD cases vs INTERVAL controls, 95% CI)") +
  theme_ms() +
  theme(panel.grid.major.y = element_blank(),
        legend.text = element_text(size = FONT_LEGEND, hjust = 0)) +
  guides(colour = guide_legend(nrow = 1))

# %% [markdown]
# ## Figure 2D and Fig. S4: attributable fraction among carrier probands
#
# From the transmission excess relative to rare synonymous variants in LOEUF deciles 3-9 (the reference of Figure 2A), the fraction of transmitted PTVs attributable to DD is q = (OR − 1) / OR. For carrier probands, each counted once however many PTVs they inherited, the attributable percentage is the sum over carriers of 1 − (1 − q)^k, with k the number of PTVs the proband inherited, divided by the number of carriers. Gene sets are Samocha *et al.*'s constrained genes (pLI ≥ 0.9, matched by Ensembl gene ID): DD-associated, other constrained, and all constrained pooled. Black points are all LOFTEE HC PTVs, red points unrescued PTVs; beside each row, the number of carrier probands the unrescued estimate represents. Intervals come from 2 × 10⁵ Monte Carlo draws on the PTV and synonymous transmission rates. Fig. S4 repeats the calculation for unrescued PTVs in NMD target regions.
#
# > **Per-individual data:** the number of carrier probands, and for each the number of inherited PTVs k, are counted within probands: read from `attributable_fraction_carriers.tsv` (probands by study, PTV definition, gene set and k), with the number of trios per study from `cohort_trios.tsv`. Transmitted and untransmitted counts come from the variant table.

# %%
## ---- gene sets: Samocha et al.'s constrained genes ---------------------------
if (!file.exists(PATHS$samocha_st10))
  stop("Samocha et al.'s Supplementary Table 10 (", PATHS$samocha_st10, ") is needed ",
       "for Figure 2D's gene sets.")
af_st10 <- fread(PATHS$samocha_st10)
if (!"symbol" %in% names(af_st10)) stop(PATHS$samocha_st10, " has no symbol column.")
af_con_hi <- GENES[!is.na(pLI) & pLI >= 0.9]
af_id_col <- names(af_st10)[vapply(af_st10, function(v)
  is.character(v) && mean(grepl("^ENSG[0-9]+", v)) > 0.9, logical(1))][1]
if (!is.na(af_id_col)) {
  af_matched_on <- "Ensembl gene ID"
  af_st10[, gid := sub("\\..*$", "", get(af_id_col))]
  af_ids <- intersect(unique(af_st10$gid), sub("\\..*$", "", af_con_hi$gene_id))
  ## each gene matched to the table under both her symbol and gnomAD's
  af_names <- unique(rbind(af_st10[gid %in% af_ids, .(gid, name = symbol)],
                           af_con_hi[sub("\\..*$", "", gene_id) %in% af_ids,
                                     .(gid = sub("\\..*$", "", gene_id), name = gene)]))
} else {
  af_matched_on <- "symbol"
  warning("No Ensembl gene IDs in her Supplementary Table 10: matching ",
          "Samocha et al.'s genes to gnomAD v2.1.1 pLI by symbol, which can add a few genes.")
  af_ids <- intersect(unique(af_st10$symbol), af_con_hi$gene)
  af_names <- data.table(gid = af_ids, name = af_ids)
}
DD_GENES <- GENES[gene_set == GS_DD, unique(gene)]
af_dd_ids <- af_names[name %in% DD_GENES, unique(gid)]
AF_SETS <- list(`All constrained`   = unique(af_names$name),
                `DD-associated`     = unique(af_names[gid %in% af_dd_ids, name]),
                `Other constrained` = unique(af_names[!gid %in% af_dd_ids, name]))
AF_N_GENES <- c(`All constrained` = length(af_ids), `DD-associated` = length(af_dd_ids),
                `Other constrained` = length(af_ids) - length(af_dd_ids))
cat(sprintf(paste0("\nFigure 2D gene sets (Samocha et al.'s constrained genes, matched by %s): ",
                   "%s constrained (her %s), %s DD-associated, %s other constrained; ",
                   "%d of our DD-associated genes are outside them.\n"),
            af_matched_on, n_fmt(AF_N_GENES[["All constrained"]]),
            n_fmt(SAMOCHA_N_GENES[["cons"]]), n_fmt(AF_N_GENES[["DD-associated"]]),
            n_fmt(AF_N_GENES[["Other constrained"]]),
            length(setdiff(DD_GENES, af_names$name))))

# %%
## ---- records: Figure 2A's unrescued PTVs and reference, Figure 1B's HC PTVs ---
af_fc <- unique(sam_par_all, by = c("CHROM", "POS", "REF", "ALT", "obs_id"))
af_tab <- rbind(
  overtransmission_table_leouf[class %in% c("Unrescued PTVs", "Synonymous"),
                               .(vepSYMBOL, Transmitted = as.logical(Transmitted),
                                 study = as.character(study), class = as.character(class),
                                 leouf_decile, rescued = FALSE)],
  af_fc[, .(vepSYMBOL, Transmitted = as.logical(Transmitted),
            study = as.character(study), class = "HC PTVs", leouf_decile = NA_character_,
            rescued = as.logical(rescued))])
AF_REF  <- sub("\n.*$", "", REF_LABEL)          # panel A's reference bin
AF_PTVS <- c(`Unrescued PTVs` = "Unrescued PTVs", `All PTVs` = "HC PTVs")
stopifnot(AF_REF %in% af_tab$leouf_decile)

## carrier probands: for each study, PTV definition and gene set, how many
## probands inherited k of the PTVs (the pooled cohort sums the studies)
AF_CARRIERS  <- fread(PATHS$af_carriers)
COHORT_TRIOS <- fread(PATHS$cohort_trios)
af_k_pro <- function(co, pv, gn) {
  d <- AF_CARRIERS[ptvs == pv & gene_set == gn]
  if (co != "DDD + RUMC") d <- d[study == co]
  d[, .(n = sum(n_probands)), keyby = k]
}
af_probands <- function(co)
  if (co == "DDD + RUMC") sum(COHORT_TRIOS$n_trios) else COHORT_TRIOS[study == co, n_trios]

# %%
## Counts as panel A counts them (rows of the table); k_pro: how many of these
## PTVs each carrier proband inherited
af_counts <- function(d, gn, cls, k_pro) {
  p <- d[class == cls & vepSYMBOL %in% AF_SETS[[gn]]]
  s <- d[class == "Synonymous" & leouf_decile == AF_REF]
  list(T = p[Transmitted == TRUE, .N], U = p[Transmitted == FALSE, .N],
       sT = s[Transmitted == TRUE, .N], sU = s[Transmitted == FALSE, .N],
       carriers = sum(k_pro$n), k_pro = k_pro)
}
## expected probands with at least one attributable inherited PTV: q a value or
## a vector of draws; ktab: k and the number of probands n
af_units <- function(ktab, q) {
  q <- pmin(pmax(q, -0.99), 0.99)
  acc <- numeric(length(q))
  for (i in seq_len(nrow(ktab))) acc <- acc + ktab$n[i] * (1 - (1 - q)^ktab$k[i])
  acc
}
af_q <- function(v, p) quantile(v[is.finite(v)], p, names = FALSE, na.rm = TRUE)
## counts as doubles throughout: U x S_T overflows R's integers
af_direct <- function(T_, U_, sT, sU, N, k_pro = NULL, B = N_MC) {
  T_ <- as.numeric(T_); U_ <- as.numeric(U_); sT <- as.numeric(sT); sU <- as.numeric(sU)
  or <- (T_ / U_) / (sT / sU); ex <- T_ * (1 - 1 / or)
  or_ci <- or * exp(c(-1, 1) * qnorm(.975) * sqrt(1 / T_ + 1 / U_ + 1 / sT + 1 / sU))
  t  <- rbeta(B, T_ + .5, U_ + .5); ts <- rbeta(B, sT + .5, sU + .5)
  or_b <- (t / (1 - t)) / (ts / (1 - ts))
  ex_b <- (T_ + U_) * t * (1 - 1 / or_b)
  out <- data.table(transmitted = T_, untransmitted = U_, syn_transmitted = sT,
                    syn_untransmitted = sU, or = or, or_lo = or_ci[1], or_hi = or_ci[2],
                    excess = ex, excess_lo = af_q(ex_b, .025), excess_hi = af_q(ex_b, .975),
                    pct = 100 * ex / N, pct_lo = 100 * af_q(ex_b, .025) / N,
                    pct_hi = 100 * af_q(ex_b, .975) / N,
                    pct_carrier_proband = NA_real_, pct_carrier_proband_lo = NA_real_,
                    pct_carrier_proband_hi = NA_real_)
  if (!is.null(k_pro)) {
    v <- 100 * af_units(k_pro, c(ex, ex_b) / T_) / sum(k_pro$n)
    out[, `:=`(pct_carrier_proband = v[1], pct_carrier_proband_lo = af_q(v[-1], .025),
               pct_carrier_proband_hi = af_q(v[-1], .975))]
  }
  out
}

# %%
## every gene set, PTV definition and cohort (the panel uses DDD and RUMC
## together; the cohorts are logged, and all are drawn so the random draws
## behind the panel are reproducible)
AF_COHORTS <- c("DDD + RUMC", STUDIES)
set.seed(20241)
af_res <- rbindlist(lapply(names(AF_PTVS), function(pv) rbindlist(lapply(AF_COHORTS, function(co) {
  in_co <- function(d) if (co == "DDD + RUMC") d else d[as.character(study) == co]
  n_pro <- af_probands(co)   # carriers or not
  d <- in_co(af_tab)
  rbindlist(lapply(names(AF_SETS), function(gn) {
    k <- af_counts(d, gn, AF_PTVS[[pv]], af_k_pro(co, pv, gn))
    cbind(data.table(ptvs = pv, cohort = co, gene_set = gn, genes = AF_N_GENES[[gn]],
                     probands = n_pro, carrier_probands = k$carriers),
          af_direct(k$T, k$U, k$sT, k$sU, n_pro, k_pro = k$k_pro))
  }))
}))))
## excess = T - U S_T / S_U exactly, and the two gene sets add up to the pooled
## one (they share the reference)
af_chk <- af_res[, transmitted - untransmitted * syn_transmitted / syn_untransmitted]
stopifnot(max(abs(af_res$excess - af_chk)) < 1e-6 * max(1, abs(af_chk)),
          af_res[, abs(sum(excess[gene_set != "All constrained"]) -
                         excess[gene_set == "All constrained"]) <
                   1e-6 * max(1, abs(excess)), by = .(ptvs, cohort)]$V1)

# %%
## ---- panel D: carrier probands, by gene set, both PTV definitions ----------------
AFD_SETS <- c(`DD-associated` = "DD-associated", `Other constrained` = "Other constrained",
              `All constrained` = "Constrained (pooled)")
AFD_LABS <- c(all = "All LOFTEE HC PTVs", unres = "Unrescued PTVs")
AFD_COLS <- c(all = "black", unres = CLASS_COLS[["Unrescued PTVs"]])
AFD_X    <- "Percentage of carrier probands whose DD is\nattributable to an inherited PTV"
AFD_OFF  <- 0.17   # each gene set's pair sits this far either side of its row
AFD_PAD  <- 0.6    # space beyond the outer rows (row units)
pd_af <- af_res[cohort == "DDD + RUMC" & gene_set %in% names(AFD_SETS)]
pd_af[, `:=`(who = fifelse(ptvs == "All PTVs", "all", "unres"),
             row = factor(AFD_SETS[gene_set], levels = AFD_SETS),
             x = pct_carrier_proband, x_lo = pct_carrier_proband_lo,
             x_hi = pct_carrier_proband_hi)]
pd_af[, y := match(as.character(row), rev(AFD_SETS)) + fifelse(who == "all", AFD_OFF, -AFD_OFF)]
afd_neg <- any(pd_af$x_lo < 0, na.rm = TRUE)
panel_D <- ggplot(pd_af, aes(x, y, colour = who)) +
  (if (afd_neg) geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5)) +
  ## room for the labels right of the intervals
  geom_blank(data = data.table(x = c(min(0, pd_af$x_lo), 1.2 * max(pd_af$x_hi))),
             aes(x = x), inherit.aes = FALSE) +
  ## bar-ends as panel A's, ~0.06 in long at the Figure 2 layout:
  ## 0.066 of the 3.2-unit axis over a 2.9 in panel
  geom_errorbar(aes(xmin = x_lo, xmax = x_hi), orientation = "y", width = 0.066,
                linewidth = 0.8) +
  geom_point(size = 5) +
  geom_text(aes(x = x_hi, label = sprintf("  %.1f%%", x)), hjust = 0, size = ANNOT_SMALL,
            colour = "grey20", fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = AFD_COLS, breaks = names(AFD_LABS), labels = unname(AFD_LABS),
                      name = NULL) +
  scale_y_continuous(breaks = seq_along(AFD_SETS), labels = rev(unname(AFD_SETS)),
                     limits = c(1 - AFD_PAD, length(AFD_SETS) + AFD_PAD),
                     expand = expansion(0)) +
  scale_x_continuous(breaks = breaks_pretty(5),
                     expand = expansion(mult = c(if (afd_neg) 0.02 else 0, 0))) +
  labs(x = AFD_X, y = NULL) +
  theme_ms() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        plot.margin = margin(6, 14, 6, 6))

## Beside each gene set's row (centred on its axis label): a small crowd of
## proband figures and the number of probands in our cohorts the red
## (unrescued PTVs) estimate represents: the estimate times the carrier
## probands, each counted once. The figure is the proband silhouette of the
## pedigree icons, traced (height 1, centred), so it stays vector in the SVG.
PROBAND_ICON <- data.table(
  x = c(0.063, 0.042, 0.032, 0.019, 0.017, 0.013, 0.011, 0.007, -0.000, -0.004, -0.007, -0.009, -0.013, -0.021, -0.037, -0.053, -0.069, -0.084, -0.093, -0.102, -0.106, -0.110, -0.113, -0.115, -0.118, -0.121, -0.127, -0.138, -0.154, -0.169, -0.180, -0.188, -0.193, -0.192, -0.186, -0.166, -0.143, -0.127, -0.080, -0.021, 0.079, 0.132, 0.148, 0.164, 0.181, 0.187, 0.192, 0.193, 0.192, 0.182, 0.164, 0.143, 0.127, 0.120, 0.116, 0.113, 0.110, 0.106, 0.102, 0.097, 0.092, 0.079, 0.063, -0.070, -0.066, -0.060, -0.043, -0.027, 0.016, 0.032, 0.042, 0.060, 0.069, 0.071, 0.064, 0.054, 0.042, 0.021, -0.011, -0.037, -0.060, -0.070, -0.070),
  y = c(0.000, 0.004, 0.011, 0.028, 0.038, 0.086, 0.308, 0.366, 0.378, 0.371, 0.345, 0.080, 0.028, 0.016, 0.005, 0.002, 0.004, 0.012, 0.022, 0.049, 0.154, 0.694, 0.688, 0.657, 0.477, 0.435, 0.426, 0.418, 0.416, 0.418, 0.426, 0.440, 0.472, 0.757, 0.778, 0.810, 0.827, 0.833, 0.837, 0.838, 0.837, 0.833, 0.828, 0.817, 0.794, 0.778, 0.709, 0.487, 0.450, 0.429, 0.417, 0.417, 0.430, 0.445, 0.646, 0.688, 0.694, 0.191, 0.059, 0.028, 0.017, 0.004, 0.000, 0.942, 0.963, 0.974, 0.991, 0.997, 1.000, 0.997, 0.990, 0.969, 0.947, 0.921, 0.900, 0.884, 0.872, 0.861, 0.858, 0.866, 0.889, 0.921, 0.942),
  id = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
## n: the count, as text; drawn "n probands" with n in bold, in the crowd's colour
afd_crowd <- function(n, col, k = 5, step = 0.27, fs = FONT_BASE * 0.72) {
  x0 <- PROBAND_ICON$x - min(PROBAND_ICON$x); w <- max(x0)
  figs <- lapply(seq_len(k), function(i)
    polygonGrob(x = unit(0, "npc") + unit(x0 + (i - 1) * step, "snpc"),
                y = unit(0.5, "npc") + unit(PROBAND_ICON$y - 0.5, "snpc"),
                id = PROBAND_ICON$id, gp = gpar(fill = col, col = "white", lwd = 0.9)))
  txt <- textGrob(bquote(bold(.(n)) * " probands"),
                  x = unit(0, "npc") + unit(w + (k - 1) * step + 0.3, "snpc"),
                  y = unit(0.5, "npc"), hjust = 0,
                  gp = gpar(fontsize = fs, col = col, fontfamily = "sans"))
  do.call(grobTree, c(figs, list(txt)))
}
afd_icons <- pd_af[who == "unres", .(row, y = match(as.character(row), rev(AFD_SETS)),
                                     n = round(pct_carrier_proband / 100 * carrier_probands))]
## its own cell beside panel D (Figure 2 assembly), with panel D's y range;
## each crowd is 2 x AFD_ICON_H rows tall
AFD_ICON_H <- 0.24
panel_D_icons <- ggplot() +
  scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
  scale_y_continuous(limits = c(1 - AFD_PAD, length(AFD_SETS) + AFD_PAD), expand = expansion(0)) +
  coord_cartesian(clip = "off") +   # a long count may run into the column gap
  theme_void() + theme(plot.margin = margin(6, 6, 6, 0))
for (i in seq_len(nrow(afd_icons)))
  panel_D_icons <- panel_D_icons +
  annotation_custom(afd_crowd(n_fmt(afd_icons$n[i]), AFD_COLS[["unres"]]),
                    xmin = 0, xmax = 1, ymin = afd_icons$y[i] - AFD_ICON_H,
                    ymax = afd_icons$y[i] + AFD_ICON_H)

# %%
cat("\nAttributable fraction (Figure 2D):\n")
af_line <- function(r) sprintf(paste0(
  "  %-15s %-17s %s genes | OR %.2f [%.2f - %.2f] | excess %.0f [%.0f - %.0f] = PAR %.2f%% ",
  "[%.2f - %.2f] of %s probands | carrier probands %.1f%% [%.1f - %.1f] of %s (~%s probands)"),
  r$ptvs, r$gene_set, n_fmt(r$genes), r$or, r$or_lo, r$or_hi, r$excess, r$excess_lo,
  r$excess_hi, r$pct, r$pct_lo, r$pct_hi, n_fmt(r$probands), r$pct_carrier_proband,
  r$pct_carrier_proband_lo, r$pct_carrier_proband_hi, n_fmt(r$carrier_probands),
  n_fmt(round(r$pct_carrier_proband / 100 * r$carrier_probands)))
for (co in AF_COHORTS) {
  cat(" ", co, "\n")
  for (i in which(af_res$cohort == co)) cat(af_line(af_res[i]), "\n")
}

# %%
## ---- Fig. S4: panel D with unrescued PTVs in NMD target regions ---------------
## Panel D's two estimates (af_res, above) without the proband counts, and a
## third: the same calculation (af_counts(), af_direct(), Figure 2A's reference,
## the same genes) on the unrescued PTVs whose variant lies in an NMD target
## region. A variant takes the NMD region of one of the transcript rows that
## make it an unrescued PTV: an NMD target row if it has one, then long exon,
## distal and start-proximal, and unplaced rows last. It counts if that region
## is NMD target, i.e. not start-proximal, long exon or distal. Carrier
## probands: those who inherited at least one unrescued PTV in an NMD target
## region, a subset of panel D's.
AFN_KEEP  <- "NMD-target"
AFN_CLASS <- "Unrescued PTVs, NMD target"
REGION_RANK <- c(`NMD-target` = 4, `Long exon` = 3, Distal = 2, `Start-proximal` = 1)
afn_ids <- unique(overtransmission_table_leouf[class == "Unrescued PTVs", obs_id])
afn_region <- par[obs_id %in% afn_ids & is_ptv == TRUE & rescue %in% "not_rescued" &
                    vepLoF %in% "HC",
                  .(nmd_region = nmd_region[which.max(fcoalesce(REGION_RANK[nmd_region], 0))]),
                  by = obs_id]
message(sprintf("Fig. S4: unrescued PTVs (variant x family): %s; %s", n_fmt(nrow(afn_region)),
                paste(afn_region[, .N, by = nmd_region][
                  order(-N), sprintf("%s %s (%.1f%%)", nmd_region, n_fmt(N),
                                     100 * N / nrow(afn_region))], collapse = ", ")))
afn_tab <- rbind(af_tab, overtransmission_table_leouf[
  class == "Unrescued PTVs" & obs_id %in% afn_region[nmd_region %in% AFN_KEEP, obs_id],
  .(vepSYMBOL, Transmitted = as.logical(Transmitted), study = as.character(study),
    class = AFN_CLASS, leouf_decile, rescued = FALSE)])
set.seed(20246)
afn_res <- rbindlist(lapply(AF_COHORTS, function(co) {
  in_co <- function(d) if (co == "DDD + RUMC") d else d[as.character(study) == co]
  n_pro <- af_probands(co)
  d <- in_co(afn_tab)
  rbindlist(lapply(names(AF_SETS), function(gn) {
    k <- af_counts(d, gn, AFN_CLASS, af_k_pro(co, AFN_CLASS, gn))
    cbind(data.table(ptvs = AFN_CLASS, cohort = co, gene_set = gn, genes = AF_N_GENES[[gn]],
                     probands = n_pro, carrier_probands = k$carriers),
          af_direct(k$T, k$U, k$sT, k$sU, n_pro, k_pro = k$k_pro))
  }))
}))
afn_all <- rbind(af_res, afn_res)
## the NMD target PTVs are a subset of the unrescued ones
afn_chk <- merge(af_res[ptvs == "Unrescued PTVs", .(cohort, gene_set, T_u = transmitted,
                                                    U_u = untransmitted, c_u = carrier_probands)],
                 afn_res[, .(cohort, gene_set, T_n = transmitted, U_n = untransmitted,
                             c_n = carrier_probands)], by = c("cohort", "gene_set"))
stopifnot(nrow(afn_chk) == nrow(afn_res), afn_chk[, all(T_n <= T_u & U_n <= U_u & c_n <= c_u)])

# %%
AFN_SERIES <- c(all = "All PTVs", unres = "Unrescued PTVs", nmd = AFN_CLASS)
AFN_LABS   <- c(all = "All LOFTEE HC PTVs", unres = "Unrescued PTVs",
                nmd = "Unrescued PTVs, NMD target region")
AFN_COLS   <- c(all = "black", unres = CLASS_COLS[["Unrescued PTVs"]],
                nmd = CLASS_COLS[["Unrescued PTVs"]])
## diamond: the triangle already means "clinically interpreted variants
## removed" in panel A
AFN_SHAPES <- c(all = 19, unres = 19, nmd = 18)
AFN_SIZES  <- c(all = 5, unres = 5, nmd = 6.2)   # a diamond draws smaller than a circle
AFN_OFF    <- c(all = 0.25, unres = 0, nmd = -0.25)   # within a gene set, top to bottom
AFN_PAD    <- 0.62
pd_afn <- afn_all[cohort == "DDD + RUMC" & gene_set %in% names(AFD_SETS)]
pd_afn[, who := factor(names(AFN_SERIES)[match(ptvs, AFN_SERIES)], levels = names(AFN_SERIES))]
pd_afn[, `:=`(y = match(AFD_SETS[gene_set], rev(AFD_SETS)) + AFN_OFF[as.character(who)],
              x = pct_carrier_proband, x_lo = pct_carrier_proband_lo,
              x_hi = pct_carrier_proband_hi)]
afn_neg <- any(pd_afn$x_lo < 0, na.rm = TRUE)
p_af_nmd <- ggplot(pd_afn, aes(x, y, colour = who, shape = who)) +
  (if (afn_neg) geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5)) +
  geom_blank(data = data.table(x = c(min(0, pd_afn$x_lo), 1.2 * max(pd_afn$x_hi))),
             aes(x = x), inherit.aes = FALSE) +
  ## bar-ends as panel A's, ~0.06 in long: 0.062 of the 3.24-unit axis over
  ## this figure's 3.15 in panel
  geom_errorbar(aes(xmin = x_lo, xmax = x_hi), orientation = "y", width = 0.062,
                linewidth = 0.8) +
  geom_point(aes(size = who)) +
  geom_text(aes(x = x_hi, label = sprintf("  %.1f%%", x)), hjust = 0, size = ANNOT_SMALL,
            colour = "grey20", fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = AFN_COLS, breaks = names(AFN_LABS), labels = unname(AFN_LABS),
                      name = NULL) +
  scale_shape_manual(values = AFN_SHAPES, breaks = names(AFN_LABS), labels = unname(AFN_LABS),
                     name = NULL) +
  scale_size_manual(values = AFN_SIZES, guide = "none") +
  scale_y_continuous(breaks = seq_along(AFD_SETS), labels = rev(unname(AFD_SETS)),
                     limits = c(1 - AFN_PAD, length(AFD_SETS) + AFN_PAD), expand = expansion(0)) +
  scale_x_continuous(breaks = breaks_pretty(5),
                     expand = expansion(mult = c(if (afn_neg) 0.02 else 0, 0))) +
  labs(x = AFD_X, y = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = unname(AFN_SIZES))),
         shape = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_ms() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        plot.margin = margin(6, 14, 6, 6))
save_figure("figS4", p_af_nmd, 10, 4.9)
cat("\nFig. S4, carrier probands whose DD is attributable to an inherited PTV (DDD and RUMC):\n")
print(pd_afn[order(match(gene_set, names(AFD_SETS)), who),
             .(gene_set, ptvs, carrier_probands, transmitted, untransmitted,
               OR = sprintf("%.2f [%.2f - %.2f]", or, or_lo, or_hi),
               attributable = sprintf("%.1f%% [%.1f - %.1f]", x, x_lo, x_hi),
               probands = round(x / 100 * carrier_probands))], row.names = FALSE)

# %% [markdown]
# ## Fig. S3: rescued PTVs among untransmitted parental alleles, by gnomAD allele frequency cap
#
# The percentage of LOFTEE HC PTVs predicted to be rescued among untransmitted parental alleles, as the gnomAD allele frequency cap rises from Samocha *et al.*'s threshold (0.001%, as in Figure 2B) to 0.1%. Variants absent from gnomAD are counted at every cap. Genes are Samocha *et al.*'s analysed genes: DD-associated and other constrained (pLI ≥ 0.9), and unconstrained (pLI < 0.9, not DD-associated). Error bars are 95% Jeffreys intervals.

# %%
## ---- Fig. S3: rescued share among untransmitted alleles, by gnomAD AF cap --------
## Constrained genes: Figure 2B's records (sam_par13, above) at each cap: genes
## of Samocha et al.'s Supplementary Table 10 with pLI >= 0.9 (sam_universe),
## one record per variant per family. The cap is applied before the
## one-record step, as the records behind Figure 2B were built.
## Unconstrained genes: the same records in the genes Figure 2B calibrates
## against, her Supplementary Table 10 genes with pLI < 0.9 (not
## DD-associated), counted on their own so a variant annotated to genes in both
## sets cannot move between them.
AF_VK <- c("CHROM", "POS", "REF", "ALT", "obs_id")
af_uncons <- setdiff(intersect(unique(GENES[!is.na(pLI) & pLI < 0.9, gene]), af_st10$symbol),
                     DD_GENES)
af_cap_records <- function(gene_list, cap)
  unique(sam_par_all[vepSYMBOL %in% gene_list & gnomad_af < cap], by = AF_VK)
af_cap <- rbindlist(lapply(AF_CAPS, function(cap) {
  d <- af_cap_records(sam_universe, cap)
  d[, stratum := fifelse(gene_set == GS_DD, GS_DD, GS_HIGH)]
  d <- rbind(d, copy(d)[, stratum := GS_ALL],
             af_cap_records(af_uncons, cap)[, stratum := GS_LOW])
  d[Transmitted == FALSE, .(n_unt = .N, k_unt = sum(rescued)), by = stratum][, cap := cap]
}))
af_cap[, `:=`(r_unt = k_unt / n_unt, r_unt_lo = sam_jeff(k_unt, n_unt, .025),
              r_unt_hi = sam_jeff(k_unt, n_unt, .975))]
## at Samocha et al.'s threshold the constrained-gene fractions are Figure 2B's
af_cap_chk <- merge(af_cap[cap == SAMOCHA_MAX_GNOMAD_AF],
                    sam_par_rates[cohort == COHORT_PRIMARY, .(stratum, n2 = n_unt, k2 = k_unt)],
                    by = "stratum")
if (nrow(af_cap_chk) != 3L || af_cap_chk[, any(n_unt != n2 | k_unt != k2)])
  warning("Fig. S3 at gnomAD AF < ", SAMOCHA_MAX_GNOMAD_AF, " differs from Figure 2B's ",
          "untransmitted parental alleles.")
af_above <- unique(sam_par_all[vepSYMBOL %in% c(sam_universe, af_uncons) & Transmitted == FALSE &
                                 gnomad_af >= max(AF_CAPS)], by = AF_VK)[, .N]
if (af_above > 0)
  message(sprintf("Fig. S3: %d untransmitted PTV records at gnomAD AF >= %g are above the top cap.",
                  af_above, max(AF_CAPS)))
AF_CAP_LAB <- paste0(formatC(100 * AF_CAPS, format = "fg", digits = 1), "%")
## one panel, coloured by gene set (as Figures 1B and 2B); the three series
## dodged so their intervals do not overlap
AF_CAP_GS <- setNames(c("DD-associated", "Other constrained", "Unconstrained"),
                      c(GS_DD, GS_HIGH, GS_LOW))
sn1 <- af_cap[stratum %in% names(AF_CAP_GS)][
  , `:=`(xi = match(cap, AF_CAPS), gs = factor(stratum, levels = names(AF_CAP_GS)))]
pd_sn1 <- position_dodge(width = 0.45)
p_af_cap <- ggplot(sn1, aes(xi, 100 * r_unt, ymin = 100 * r_unt_lo, ymax = 100 * r_unt_hi,
                            colour = gs, group = gs)) +
  geom_line(linewidth = 0.8, position = pd_sn1) +
  geom_errorbar(width = 0.3, linewidth = 0.8, position = pd_sn1) +
  geom_point(size = 4, position = pd_sn1) +
  scale_colour_manual(values = GENE_SET_COLOURS[names(AF_CAP_GS)], breaks = names(AF_CAP_GS),
                      labels = unname(AF_CAP_GS), name = NULL) +
  scale_x_continuous(breaks = seq_along(AF_CAPS), labels = AF_CAP_LAB,
                     limits = c(0.6, length(AF_CAPS) + 0.4), expand = expansion(0)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "gnomAD allele frequency cap",
       y = "LOFTEE HC PTVs rescued\n(untransmitted alleles, 95% CI)") +
  theme_ms() +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
save_figure("figS3", p_af_cap, 9, 6)
cat("\nFig. S3, LOFTEE HC PTVs rescued among untransmitted parental alleles (%):\n")
print(dcast(af_cap[, .(stratum, cap = factor(AF_CAP_LAB[match(cap, AF_CAPS)], levels = AF_CAP_LAB),
                       pct = round(100 * r_unt, 1))], stratum ~ cap, value.var = "pct"))

# %% [markdown]
# ## Fig. S10: families segregating two PTVs in DD-associated genes
#
# Families in which the parents carry two unrescued PTVs in DD-associated genes, arranged by which parent carried them and how many were transmitted, and coloured by the child's diagnostic status (diagnosis attributed to the inherited PTV, *de novo* diagnosis, or undiagnosed). Here an inherited PTV counts as diagnostic only when it was clinically reported *and* the proband has a diagnosis record. Squares are lighter where both PTVs are in the same parent on the same chromosome, so their transmissions may not be independent.
#
# > **Per-individual data:** each family's PTVs are grouped (number carried, which parent, how many transmitted, same chromosome) and the child's status joined: read from `two_ptv_families.tsv` (families by these categories, as in the published figure).

# %%
lev_src <- c("Both carried by mother", "Both carried by father",
             "One carried by each parent")
lev_tx  <- c("Both transmitted", "One transmitted", "Neither transmitted")
lev_dx  <- c("Inherited PTV is diagnostic", "De novo diagnosis", "Undiagnosed")
## one row per family
two <- fread(PATHS$two_ptv)
two <- two[rep(seq_len(.N), n), .(source, transmission, independent, dx_cat = dx_category)]
two[, source := factor(source, levels = lev_src)]
two[, transmission := factor(transmission, levels = lev_tx)]
two[, dx_cat := factor(dx_cat, levels = lev_dx)]
grp_n <- two[, .(N = .N), by = .(transmission, source)]

# %%
dx_cols <- c("Inherited PTV is diagnostic" = "#1D3557",
             "De novo diagnosis" = "#E9A44C", "Undiagnosed" = "grey78")
ncol_w <- 5
waffle <- copy(two)
waffle[, indep_lab := factor(
  fifelse(independent, "Independent transmissions",
          "Same parent & chromosome (may not be independent)"),
  levels = c("Independent transmissions",
             "Same parent & chromosome (may not be independent)"))]
waffle[, ord := as.integer(dx_cat)]
setorder(waffle, transmission, source, -ord, independent)
waffle[, i := seq_len(.N) - 1L, by = .(transmission, source)]
waffle[, `:=`(wx = i %% ncol_w, wy = i %/% ncol_w)]
n_lab <- copy(grp_n)[, `:=`(wx = (ncol_w - 1) / 2, wy = max(waffle$wy) + 1.1)]

fig_s10 <- ggplot(waffle, aes(wx, wy, fill = dx_cat, alpha = indep_lab)) +
  geom_tile(colour = "white", linewidth = 1.1, width = .92, height = .92) +
  geom_text(data = n_lab, aes(wx, wy, label = paste0("n = ", N)),
            inherit.aes = FALSE, size = ANNOT_SMALL, fontface = "italic") +
  facet_grid(transmission ~ source,
             labeller = labeller(source = setNames(sub(" (mother|father|each parent)$", "\n\\1",
                                                       lev_src), lev_src))) +
  scale_fill_manual(values = dx_cols, breaks = lev_dx) +
  scale_alpha_manual(values = c("Independent transmissions" = 1,
                                "Same parent & chromosome (may not be independent)" = 0.42),
                     drop = FALSE) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 1)), alpha = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  coord_equal(clip = "off") +
  labs(x = NULL, y = NULL,
       title = "Diagnostic status of probands in families with n = 2 PTVs in a DD-associated gene") +
  theme_ms() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), panel.spacing = unit(12, "pt"),
        strip.clip = "off",
        strip.text = element_text(size = FONT_TITLE - 1, face = "bold"),
        strip.text.y = element_text(angle = -90),
        plot.title = element_text(size = FONT_TITLE - 1, face = "bold", hjust = 0),
        plot.title.position = "plot", legend.background = element_blank(),
        legend.text = element_text(size = FONT_LEGEND - 2, margin = margin(r = 10)),
        plot.margin = margin(10, 24, 6, 6))

## Square tiles at a fixed size (TILE_IN inches per tile), with the panels flush
## left and any spare width to the right of the row labels
TILE_IN <- 0.357
pp <- ggplot_build(fig_s10)$layout$panel_params[[1]]
fig_s10 <- ggplotGrob(fig_s10)
pan <- fig_s10$layout[grepl("^panel", fig_s10$layout$name), ]
fig_s10$widths[unique(pan$l)]  <- unit(diff(pp$x.range) * TILE_IN, "in")
fig_s10$heights[unique(pan$t)] <- unit(diff(pp$y.range) * TILE_IN, "in")
fig_s10$respect <- FALSE
fig_s10 <- gtable::gtable_add_cols(fig_s10, unit(1, "null"))
fig_s10 <- gtable::gtable_add_rows(fig_s10, unit(1, "null"))
save_figure("figS10", fig_s10, 8, 7.3)

# %% [markdown]
# ## Fig. S6 and Fig. S7: parent of origin
#
# Fig. S6: among clinically interpreted unrescued PTVs in DD-associated genes (pathogenic or likely pathogenic, and uncertain), the number transmitted by each parent, against the 50:50 expectation (dashed line). Fig. S7: for variants that were not clinically reported, a logistic regression of transmission on the interaction of variant class, LOEUF bin and transmitting parent, excluding probands with a clinically interpreted inherited PTV (`inherited_ptv_diagnosis`, as for Figure 2A's triangles). The mother-to-father odds ratio in the reference stratum is the baseline, and the remaining estimates are relative to it; stars are Bonferroni-adjusted.

# %%
H_LEVELS <- c("Pathogenic/likely pathogenic", "Uncertain")
## one record per variant and family
poh <- unique(overtransmission_table_leouf[
  class == "Unrescued PTVs" & leouf_decile == GS_DD & Transmitted == TRUE &
    clinical_interpretation %in% H_LEVELS,
  .(study, obs_id, interp = clinical_interpretation, parent = Relationship)], by = "obs_id")
if (poh[is.na(parent), .N] > 0)
  stop("Records without a transmitting parent in the parent-of-origin set.")
poh[, interp := factor(interp, levels = H_LEVELS)]

H_dat <- rbind(poh[, .(n = .N), by = .(interp, parent)][, study := "All"],
               poh[, .(n = .N), by = .(study, interp, parent)])
setorder(H_dat, study, interp, parent)
half_lines <- H_dat[, .(half = sum(n) / 2), by = .(study, interp)]

fig_s6 <- ggplot(H_dat[study == "All"], aes(interp, n, fill = parent)) +
  geom_col(width = BAR_WIDTH_H) +
  geom_errorbar(data = half_lines[study == "All"],
                aes(x = interp, ymin = half, ymax = half), inherit.aes = FALSE,
                linetype = "dashed", colour = "black", linewidth = 0.8,
                width = BAR_WIDTH_H) +
  scale_fill_manual(values = c("Father" = "steelblue4", "Mother" = "pink3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c("Pathogenic/likely pathogenic" =
                                "Pathogenic/\nlikely pathogenic",
                              "Uncertain" = "Uncertain")) +
  ylab("Number of clinically interpreted\ninherited PTVs") +
  xlab("Clinical variant interpretation") +
  theme_ms() + theme(panel.grid.major.x = element_blank())
save_figure("figS6", fig_s6, 7, 6)

# %%
poo_all <- overtransmission_table_leouf[class %in% c("Unrescued PTVs", "Synonymous") &
                                          inherited_ptv_diagnosis == FALSE]
poo_all[, parent := Relationship]
poo_all <- poo_all[!is.na(parent)]
poo_all[, label := relevel(droplevels(label), ref = REF_LABEL)]
d_poo <- tidy_glm(glm(Transmitted ~ label * parent, binomial("logit"), poo_all))
d_poo <- d_poo[grepl("parentMother", term)]
d_poo[term == "parentMother", term := REF_LABEL]
d_poo[, term := gsub(":parentMother", "", term, fixed = TRUE)]
d_poo[, term := sub("^label", "", term)]
d_poo <- split_bin_type(d_poo)
d_poo[, type := factor(type, levels = CLASS_LEVELS)]
d_poo <- d_poo[!is.na(bin) & !is.na(type)]
d_poo[, type := droplevels(type)]
add_sig(d_poo)

pd5 <- position_dodge(width = 0.8)
fig_s7 <- ggplot(d_poo, aes(bin, estimate, ymin = conf.low, ymax = conf.high,
                              colour = type, group = type)) +
  geom_hline(yintercept = 1, colour = "black", linewidth = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.8, position = pd5) +
  geom_point(size = 5, position = pd5) + sig_layer(pd5) + coord_flip() +
  scale_colour_manual(values = CLASS_COLS, drop = FALSE) +
  scale_y_continuous(limits = c(0.5, 2.0), breaks = c(0.5, 1.0, 1.5, 2.0),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = decile_labels_vec) +
  xlab("LOEUF decile") +
  ylab("Estimated transmission odds from mothers versus fathers (95% CI)") +
  theme_ms()
save_figure("figS7", fig_s7, 9, 7)

# %% [markdown]
# ## Figure 2: assembly
#
# Panels A-E in one grid: A spans the rows of B and C, D sits beside its proband counts, and E shares the right column with B and C.

# %%
## ONE grid for the whole figure. Patchwork aligns plot areas only within a
## grid, so in one grid:
##   - A spans the rows of B and C, so its x axis lines up with C's;
##   - A spans the columns of D and its proband counts, so A and D share one
##     left edge, and B, C and E share the right column's;
##   - D, its proband counts and E share one row, so the counts sit level with
##     D's gene sets.
## A and D carry their own legends.
## Patchwork orders the plots by their design letters, alphabetically, and
## tags them in that order: E is panel D's proband counts (untagged), F is
## panel E.
FIG_DESIGN <- "
AAB
AAC
DEF
"
COL_WIDTHS  <- c(1.15, 1)          # left column (A; D and its counts) : right column (B, C, E)
D_SPLIT     <- c(4, 1.7)           # panel D : its proband counts, within the left column
TOP_H       <- 1.5                 # rows 1-2 (A | B, C)
BOTTOM_H    <- 0.91                # row 3 (D | E)
BC_HEIGHTS  <- c(0.8, 0.85)        # B : C within the top rows
ROW_HEIGHTS <- c(TOP_H * BC_HEIGHTS / sum(BC_HEIGHTS), BOTTOM_H)
FIG_W <- 17; FIG_H <- 13.2         # inches: A ~6.3 in tall, B and C ~2.4, D and E ~2.9

## Sharing the right column's left edge reserves room for B and C's gene-set
## labels, which would push E's y-axis title away from its tick labels.
## free(type = "label") keeps E's panel aligned but lets its title sit next to
## its own tick labels (patchwork >= 1.3.0; older versions leave the gap).
panel_E_cell <- if (packageVersion("patchwork") >= "1.3.0")
  free(panel_E, type = "label", side = "l") else panel_E

## Tags in design-letter order; the proband counts are untagged.
tag_seq <- c("A", "B", "C", "D", "", "E")
main_figure <- wrap_plots(A = panel_A, B = panel_B, C = panel_C, D = panel_D,
                          E = panel_D_icons, F = panel_E_cell, design = FIG_DESIGN,
                          widths = c(COL_WIDTHS[1] * D_SPLIT / sum(D_SPLIT), COL_WIDTHS[2]),
                          heights = ROW_HEIGHTS) +
  plot_annotation(tag_levels = list(tag_seq)) &
  theme(plot.tag = element_text(size = FONT_TAG, face = "bold"))
save_figure("fig2", main_figure, FIG_W, FIG_H)

# %%
sessionInfo()

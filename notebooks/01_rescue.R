# %% [markdown]
# # 01 · Inherited PTVs are frequently rescued in constrained genes
#
# Rescue of inherited protein-truncating variants (PTVs) in DD probands, as
# predicted by RescueRanger from phased sequence data, by gene set and VEP
# consequence; rescue in NMD target and escape regions; and how often all
# possible canonical splice variants with compensatory cryptic splice sites are
# seen in gnomAD.
#
# Figures drawn:
#
# - **Figure 1B**: proportion of LOFTEE HC inherited PTVs rescued, by gene set and consequence
#   (Figure 1A is a schematic and is not drawn here)
# - **Fig. S1**: the same by study, by gene set (A) and constrained vs unconstrained (B)
# - **Fig. S19**: canonical splice variants observed in gnomAD, by rescue probability
# - **Fig. S20**: rescue in NMD target vs NMD escape regions of transcripts
# - **Figure N2**: RNA sequencing read support for the splice junctions behind
#   rescue predictions (A), and rescue predictions against cell depletion in
#   *DDX3X* saturation genome editing (B)

# %%
REPO_DIR <- Sys.getenv("RR_REPO", unset = if (file.exists("config.R")) "." else "..")
source(file.path(REPO_DIR, "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(svglite)
  library(pROC)
})

# %% [markdown]
# ## Inputs
#
# - `proband_variants.tsv.gz` (EGA release): the probands' PTV and synonymous
#   calls, one row per variant and transcript, with the RescueRanger prediction
#   (`rescue`), the NMD region of each row's transcript and the gene set.
#   Hypermutated trios, call quality filters and NMD placement were applied
#   before release (`data/README.md`).
# - `gene_annotation.tsv.gz` (included in `resources/`): the gene set of every
#   gene, for Fig. S19.
# - Rescue probabilities for all possible canonical splice site SNVs, and their
#   gnomAD v2.1.1 allele counts, for Fig. S19 (supplementary data; no cohort data).
# - For Figure N2 (supplementary data; no cohort data): splice junction read
#   counts in Geuvadis RNA sequencing for heterozygous splice PTVs in 1000 Genomes
#   individuals, and rescue predictions for *DDX3X* saturation genome editing
#   variants.

# %%
PATHS <- list(
  probands        = file.path(DATA_DIR, "proband_variants.tsv.gz"),
  gene_annotation = file.path(RESOURCE_DIR, "gene_annotation.tsv.gz"),
  splice_rescue   = file.path(DATA_DIR, "splice_rescue_probabilities.tsv"),
  gnomad_ac       = file.path(DATA_DIR, "gnomad_ACs.tsv"),
  junctions       = file.path(DATA_DIR, "junction_table.tsv"),
  ddx3x_sge       = file.path(DATA_DIR, "DDX3X_splice_consequences.tsv.gz")
)

# %% [markdown]
# ## Theme, helpers and vocabulary

# %%
theme_rr <- function(base_size = 16) {
  theme(
    text = element_text(family = "sans", size = base_size, colour = "black"),
    panel.background = element_blank(), plot.background = element_blank(),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    axis.line.x.top = element_blank(), axis.line.y.right = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(3, "pt"),
    axis.text = element_text(size = base_size),
    legend.position = "bottom", legend.box = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = base_size - 2), legend.key = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size - 3, face = "bold"))
}

## normal-approximation 95% binomial interval half-width
binom_ci  <- function(prop, total) 1.96 * sqrt((prop * (1 - prop)) / total)
fmt_count <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)

significance_label <- function(p)
  fifelse(p <= 0.001, "***", fifelse(p <= 0.01, "**",
          fifelse(p <= 0.05, "*", fifelse(p <= 0.1, ".", "ns"))))

PTV_CONSEQUENCES <- c("frameshift_variant", "splice_acceptor_variant",
                      "splice_donor_variant", "stop_gained")
CSQ_FIGURE_LABELS <- c(splice_acceptor_variant = "Canonical splice\nacceptor variants",
                       splice_donor_variant    = "Canonical splice\ndonor variants",
                       frameshift_variant      = "Frameshift\nvariants",
                       stop_gained             = "Stop gained\nvariants")

## gene sets as ordered factors, so bars and legends follow GENE_SETS
gene_set_factor  <- function(x) factor(x, levels = GENE_SETS)
constraint_factor <- function(gene_set)
  factor(fifelse(as.character(gene_set) == GS_LOW, "Unconstrained", "Constrained"),
         levels = CONSTRAINT_LEVELS)

# %% [markdown]
# ## The probands' inherited PTVs
#
# Inherited calls (`is_denovo` FALSE) from the proband table. PTVs without a
# rescue prediction (`rescue` NA) are left out. A PTV counts as rescued when
# RescueRanger predicts rescue (rescue probability ≥ 0.50) without creating
# another PTV (`rescue == "rescued"`); `not_rescued_new_ptv` counts as not
# rescued.

# %%
probands <- fread(PATHS$probands, colClasses = list(character = "CHROM"),
                  encoding = "UTF-8")
probands[, gene_set := gene_set_factor(gene_set)]
inherited <- probands[is_denovo == FALSE & is_ptv == TRUE & !is.na(rescue)]
inherited[, rescued := rescue == "rescued"]

# %% [markdown]
# ## Figure 1B and Fig. S1: rescue by gene set and consequence
#
# The proportion of LOFTEE HC inherited PTVs predicted to be rescued, for each
# canonical PTV consequence and for all four pooled. Denominators are
# variant–transcript records. Error bars are 95% binomial intervals; brackets
# compare each constrained gene set against pLI < 0.9 by two-sided Fisher's
# exact test with Bonferroni correction (within study in Fig. S1). Fig. S1B
# pools the two constrained sets against unconstrained genes.

# %%
inherited_ptvs <- inherited[vepLoF == "HC" & vepConsequence %in% PTV_CONSEQUENCES]
inherited_ptvs[, consequence := CSQ_FIGURE_LABELS[vepConsequence]]
inherited_ptvs[, constraint  := constraint_factor(gene_set)]

CSQ_FIGURE_LEVELS <- c("All\nPTVs",
                       CSQ_FIGURE_LABELS[c("splice_acceptor_variant",
                                           "splice_donor_variant",
                                           "frameshift_variant", "stop_gained")])

## Rescue proportions with binomial CIs, plus an "All PTVs" aggregate over the
## four consequence classes. Shared by Figure 1B and Fig. S1 so the two cannot
## disagree.
summarise_rescue <- function(dt, by_cols) {
  by_csq  <- dt[, .(total = .N, rescued = sum(rescued)), by = c(by_cols, "consequence")]
  all_ptv <- dt[, .(total = .N, rescued = sum(rescued)),
                by = by_cols][, consequence := "All\nPTVs"]
  out <- rbind(all_ptv, by_csq, use.names = TRUE)
  out[, prop := rescued / total]
  out[, ci := binom_ci(prop, total)]
  out[, `:=`(ci_low = prop - ci, ci_high = prop + ci)]
  out[, consequence := factor(consequence, levels = CSQ_FIGURE_LEVELS)]
  out[]
}

RESCUE_REFERENCE <- GS_LOW
RESCUE_COMPARISONS <- list(
  list(label = paste(GS_DD, "vs", GS_LOW),   a = GS_DD,             b = RESCUE_REFERENCE),
  list(label = paste(GS_HIGH, "vs", GS_LOW), a = GS_HIGH,           b = RESCUE_REFERENCE),
  list(label = "Constrained vs unconstrained", a = c(GS_DD, GS_HIGH), b = RESCUE_REFERENCE))
PANEL_BRACKET_COMPARISONS <- RESCUE_COMPARISONS[1:2]

## Fisher's exact test of each comparison, per consequence (and per study)
rescue_comparison_table <- function(dt, by_cols = character(0)) {
  counts <- summarise_rescue(dt, c("gene_set", by_cols))
  group_cols <- c("consequence", by_cols)
  res <- rbindlist(lapply(RESCUE_COMPARISONS, function(cmp) {
    a <- counts[gene_set %in% cmp$a, .(total_a = sum(total), rescued_a = sum(rescued)),
                by = group_cols]
    b <- counts[gene_set %in% cmp$b, .(total_b = sum(total), rescued_b = sum(rescued)),
                by = group_cols]
    d <- merge(a, b, by = group_cols)
    d[, comparison := cmp$label][]
  }))
  res[, fisher_p := {
    m <- matrix(c(rescued_a, total_a - rescued_a,
                  rescued_b, total_b - rescued_b), nrow = 2, byrow = TRUE)
    fisher.test(m)$p.value
  }, by = c(group_cols, "comparison")]
  ## Bonferroni within each study; the pooled rows are their own family
  if (length(by_cols) > 0)
    res[, fisher_p_bonferroni := p.adjust(fisher_p, "bonferroni"), by = by_cols]
  else
    res[, fisher_p_bonferroni := p.adjust(fisher_p, "bonferroni")]
  res[, comparison := factor(comparison, levels = vapply(RESCUE_COMPARISONS,
                                                         `[[`, character(1), "label"))]
  setorderv(res, c("consequence", by_cols, "comparison"))
  res[]
}

DODGE_WIDTH <- 0.7
bar_offsets <- function(levels, dodge_width = DODGE_WIDTH)
  setNames((seq_along(levels) - (length(levels) + 1) / 2) * dodge_width /
             length(levels), levels)
GENE_SET_BAR_OFFSET <- bar_offsets(GENE_SETS)

comparison_bracket <- function(rescue_data, stats,
                               comparisons = PANEL_BRACKET_COMPARISONS,
                               by_cols = character(0),
                               gap = 0.02, tick = 0.008, step = 0.05) {
  tops <- rescue_data[, .(y_top = max(ci_high)), by = c("consequence", by_cols)]
  br <- rbindlist(lapply(seq_along(comparisons), function(i) {
    cmp <- comparisons[[i]]
    d <- merge(tops, stats[comparison == cmp$label,
                           c("consequence", by_cols, "fisher_p_bonferroni"),
                           with = FALSE], by = c("consequence", by_cols))
    d[, x_start := as.integer(consequence) + mean(GENE_SET_BAR_OFFSET[cmp$a])]
    d[, x_end   := as.integer(consequence) + mean(GENE_SET_BAR_OFFSET[cmp$b])]
    d[, level   := i][]
  }))
  br[, y_bar   := y_top + gap + (level - 1) * step]
  br[, y_tick  := y_bar - tick]
  br[, y_label := y_bar + tick]
  br[, label   := significance_label(fisher_p_bonferroni)]
  br[]
}

add_significance_brackets <- function(p, brackets, expand_mult = 0.18) {
  p +
    geom_segment(data = brackets, aes(x = x_start, xend = x_end, y = y_bar,
                                      yend = y_bar), inherit.aes = FALSE,
                 linewidth = 0.5) +
    geom_segment(data = brackets, aes(x = x_start, xend = x_start, y = y_bar,
                                      yend = y_tick), inherit.aes = FALSE,
                 linewidth = 0.5) +
    geom_segment(data = brackets, aes(x = x_end, xend = x_end, y = y_bar,
                                      yend = y_tick), inherit.aes = FALSE,
                 linewidth = 0.5) +
    geom_text(data = brackets, aes(x = (x_start + x_end) / 2, y = y_label,
                                   label = label), inherit.aes = FALSE,
              vjust = 0, size = 5) +
    scale_y_continuous(expand = expansion(mult = c(0, expand_mult)))
}

# %%
## Figure 1B
fig1_data <- summarise_rescue(inherited_ptvs, "gene_set")
setorder(fig1_data, consequence, gene_set)

fig1_denominators <- fig1_data[, .(n = sum(total)), by = consequence]
setorder(fig1_denominators, consequence)
FIG1_AXIS_LABELS <- setNames(
  paste0(levels(fig1_data$consequence), "\n(n = ",
         fmt_count(fig1_denominators$n), ")"), levels(fig1_data$consequence))

constraint_stats_pooled <- rescue_comparison_table(inherited_ptvs)
fig1_brackets <- comparison_bracket(fig1_data, constraint_stats_pooled)

fig1b <- ggplot(fig1_data, aes(consequence, prop, fill = gene_set)) +
  geom_col(position = position_dodge(width = DODGE_WIDTH), width = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2,
                position = position_dodge(width = DODGE_WIDTH)) +
  scale_fill_manual(values = GENE_SET_COLOURS) +
  scale_x_discrete(labels = FIG1_AXIS_LABELS) +
  xlab("VEP variant effect prediction") +
  ylab("Proportion of LOFTEE HC PTVs\nrescued in haplotype context") +
  theme_rr() +
  theme(legend.position = "top", panel.grid.major.x = element_blank())
fig1b <- add_significance_brackets(fig1b, fig1_brackets)
fig1b <- fig1b + labs(tag = "B") +
  theme(plot.tag = element_text(size = 22, face = "bold"),
        plot.tag.position = "topleft")

save_figure("fig1", fig1b, width = 12, height = 6)

# %%
## Fig. S1: panel A by gene set, panel B constrained vs unconstrained, by study
fig1_counts <- rescue_comparison_table(inherited_ptvs, "study")

f1a_data <- summarise_rescue(inherited_ptvs, c("gene_set", "study"))
setorder(f1a_data, study, consequence, gene_set)
f1a_brackets <- comparison_bracket(f1a_data, fig1_counts, by_cols = "study")

big_text <- theme(legend.position = "top", panel.grid.major.x = element_blank(),
                  axis.title = element_text(size = 24),
                  axis.text.y = element_text(size = 20),
                  axis.text.x = element_text(size = 16),
                  strip.text = element_text(size = 22, face = "bold"),
                  legend.text = element_text(size = 20))

supp_f1a <- ggplot(f1a_data, aes(consequence, prop, fill = gene_set)) +
  geom_col(position = position_dodge(width = DODGE_WIDTH), width = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2,
                position = position_dodge(width = DODGE_WIDTH)) +
  facet_wrap(. ~ study) + scale_fill_manual(values = GENE_SET_COLOURS) +
  xlab("VEP variant effect prediction") +
  ylab("Proportion of LOFTEE HC PTVs\nrescued in haplotype context") +
  theme_rr() + big_text
supp_f1a <- add_significance_brackets(supp_f1a, f1a_brackets)

f1b_data <- summarise_rescue(inherited_ptvs, c("constraint", "study"))
setorder(f1b_data, study, consequence, constraint)
f1b_labels <- merge(
  f1b_data[, .(y = max(ci_high)), by = .(study, consequence)],
  fig1_counts[comparison == "Constrained vs unconstrained",
              .(study, consequence, label = significance_label(fisher_p_bonferroni))],
  by = c("study", "consequence"))

supp_f1b <- ggplot(f1b_data, aes(consequence, prop, fill = constraint)) +
  geom_col(position = position_dodge(width = DODGE_WIDTH), width = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2,
                position = position_dodge(width = DODGE_WIDTH)) +
  geom_text(data = f1b_labels, aes(consequence, y + 0.03, label = label),
            inherit.aes = FALSE, size = 5) +
  facet_wrap(. ~ study) + scale_fill_manual(values = CONSTRAINT_COLOURS) +
  xlab("VEP variant effect prediction") +
  ylab("Proportion of LOFTEE HC PTVs\nrescued in haplotype context") +
  theme_rr() + big_text

supp_f1 <- supp_f1a / supp_f1b + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

save_figure("figS1", supp_f1, width = 20, height = 13)

# %% [markdown]
# ## Fig. S20: rescue in NMD target vs NMD escape regions
#
# Inherited PTVs with a rescue prediction and no RescueRanger UTR flag
# (`utr_flag`), with no LOFTEE filter, by the NMD region of each transcript row
# (`nmd_region`, from the per-transcript annotation of Blakes et al.). NMD
# target regions are the `NMD-target` rows; escape regions are the
# start-proximal, long exon and distal rows (distal: the last exon and the 3′
# end of the penultimate exon). Rows that could not be placed (`non_coding`,
# `unplaced`) are left out. A variant is counted once in each region class its
# transcript rows fall in (once per rescue prediction, if its rows differ). The
# proportion rescued in NMD target regions is compared with distal regions and
# with all escape regions (Fisher's exact test, Bonferroni across the
# comparisons shown).

# %%
NMD_REGION_LEVELS  <- c("All NMD target regions", "All NMD escape regions",
                        "Last exon")
NMD_REGION_COLOURS <- c("All NMD target regions" = "darkolivegreen3",
                        "All NMD escape regions" = "plum4",
                        "Last exon"              = "slategray3")
NMD_DODGE_WIDTH <- 0.9
PLACED_REGIONS <- c("NMD-target", "Start-proximal", "Long exon", "Distal")

region_ptvs <- inherited[utr_flag == FALSE & nmd_region %in% PLACED_REGIONS]
region_ptvs[, NMD_all := fifelse(nmd_region == "NMD-target", "All NMD target regions",
                                 "All NMD escape regions")]

## one record per variant (obs_id) in each region class
count_nmd_rescue <- function(dt)
  unique(dt[, .(obs_id, rescue, rescued, NMD_all, gene_set)])[
    , .(total = .N, rescued = sum(rescued)), by = .(gene_set, NMD_all)]

nmd_counts <- rbind(count_nmd_rescue(region_ptvs),
                    count_nmd_rescue(region_ptvs[nmd_region == "Distal"][
                      , NMD_all := "Last exon"]), use.names = TRUE)
nmd_counts[, prop_rescued := rescued / total]
nmd_counts[, ci := binom_ci(prop_rescued, total)]
nmd_counts[, `:=`(ci_low = prop_rescued - ci, ci_high = prop_rescued + ci)]
nmd_counts[, gene_set := factor(gene_set, levels = GENE_SETS)]
nmd_counts[, NMD_all  := factor(NMD_all, levels = NMD_REGION_LEVELS)]
setorder(nmd_counts, gene_set, NMD_all)

# %%
NMD_COMPARISONS <- c("Last exon" = "NMD target versus distal NMD escape (last exon)",
                     "All NMD escape regions" = "NMD target versus all NMD escape")

## Fisher's exact test of NMD target against each escape class, per gene set
nmd_fisher_table <- function(counts) {
  d <- copy(counts)
  d[, grp := gene_set]
  d <- d[, .(total = sum(total), rescued = sum(rescued)), by = .(grp, NMD_all)]
  res <- rbindlist(lapply(levels(droplevels(d$grp)), function(g) {
    target <- d[grp == g & NMD_all == "All NMD target regions"]
    rbindlist(lapply(names(NMD_COMPARISONS), function(escape_class) {
      escape <- d[grp == g & NMD_all == escape_class]
      row <- data.table(grp = g, escape_class = escape_class,
                        comparison = NMD_COMPARISONS[[escape_class]])
      if (nrow(target) == 0 || nrow(escape) == 0) {
        warning("No records for ", g, "/", escape_class, "; reported as NA.")
        return(cbind(row, data.table(fisher_p = NA_real_)))
      }
      m <- matrix(c(target$rescued, target$total - target$rescued,
                    escape$rescued, escape$total - escape$rescued),
                  nrow = 2, byrow = TRUE)
      cbind(row, data.table(fisher_p = fisher.test(m)$p.value))
    }))
  }))
  res[, fisher_p_bonferroni := p.adjust(fisher_p, "bonferroni")]
  res[, label := significance_label(fisher_p_bonferroni)]
  res[]
}
nmd_stats_gene_set <- nmd_fisher_table(nmd_counts)

NMD_BAR_OFFSET <- bar_offsets(NMD_REGION_LEVELS, NMD_DODGE_WIDTH)
nmd_brackets <- merge(nmd_counts[, .(y_top = max(ci_high)), by = gene_set],
                      nmd_stats_gene_set[, .(gene_set = factor(grp, levels = GENE_SETS),
                                             escape_class, label)], by = "gene_set")
nmd_brackets[, level   := match(escape_class, names(NMD_COMPARISONS))]
nmd_brackets[, x_start := as.integer(gene_set) + NMD_BAR_OFFSET[["All NMD target regions"]]]
nmd_brackets[, x_end   := as.integer(gene_set) + NMD_BAR_OFFSET[escape_class]]
nmd_brackets[, y_bar   := y_top + 0.02 + (level - 1) * 0.035]
nmd_brackets[, y_tick  := y_bar - 0.008]
nmd_brackets[, y_label := y_bar + 0.004]

fig_s20 <- ggplot(nmd_counts, aes(gene_set, prop_rescued, fill = NMD_all)) +
  geom_col(position = position_dodge(width = NMD_DODGE_WIDTH), width = 0.85) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2,
                position = position_dodge(width = NMD_DODGE_WIDTH)) +
  scale_fill_manual(values = NMD_REGION_COLOURS) + xlab("") +
  ylab("Proportion of PTVs\nrescued in haplotype context") +
  theme_rr() + theme(panel.grid.major.x = element_blank())
save_figure("figS20", add_significance_brackets(fig_s20, nmd_brackets), 12, 6)

# %% [markdown]
# ## Fig. S19: canonical splice variants in gnomAD, by rescue probability
#
# Every possible SNV at a canonical splice donor or acceptor site, with its
# rescue probability averaged across transcripts of the gene (Methods), is
# matched to gnomAD v2.1.1 (PASS sites only). Variants are binned into ten
# equal rescue probability intervals, and the proportion seen at least once in
# gnomAD is shown per gene set with 95% normal-approximation intervals. This
# uses no cohort data.

# %%
genes <- fread(PATHS$gene_annotation, encoding = "UTF-8")
## gene set of a symbol; genes absent from the annotation are pLI < 0.9
gene_set_of <- function(symbol) {
  gs <- genes$gene_set[match(symbol, genes$gene)]
  gene_set_factor(fifelse(is.na(gs), GS_LOW, gs))
}

srp <- fread(PATHS$splice_rescue)
srp[, gene_set := gene_set_of(gene)]
srp <- srp[, .(avg_rescue = mean(rescue)),
           by = c("gene_set", "chr", "pos", "ref", "variant", "var_pos")]
gnomad <- merge(srp, fread(PATHS$gnomad_ac),
                by.x = c("chr", "var_pos", "ref", "variant"),
                by.y = c("CHROM", "POS", "REF", "ALT"),
                all.x = TRUE, all.y = FALSE)
gnomad[, present_in_gnomad := !(is.na(AC) | AC == 0)]
gnomad[, AC := fifelse(is.na(AC), 0, AC)]
gnomad <- gnomad[is.na(FILTER) | FILTER == "PASS"]
gnomad[, p_rescue_bin := cut(avg_rescue, breaks = seq(0, 1, length.out = 11),
                             include.lowest = TRUE, labels = FALSE)]
BIN_LABELS_10 <- c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4", "0.4 - 0.5",
                   "0.5 - 0.6", "0.6 - 0.7", "0.7 - 0.8", "0.8 - 0.9", "0.9 - 1")

observed <- unique(gnomad[, .(total = .N, n = sum(present_in_gnomad)),
                          by = .(p_rescue_bin, gene_set)])
observed[, prop := n / total]
observed[, ci := binom_ci(prop, total)]
observed[, `:=`(min = prop - ci, max = prop + ci)]
setorder(observed, gene_set, p_rescue_bin)

fig_s19 <- ggplot(observed, aes(factor(p_rescue_bin, levels = 1:10), prop,
                                colour = gene_set, group = gene_set)) +
  geom_point(position = position_dodge(width = 0.2), size = 5) +
  geom_line(position = position_dodge(width = 0.2), linetype = "dashed",
            linewidth = 0.6) +
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.2,
                position = position_dodge(width = 0.2)) +
  scale_x_discrete(labels = BIN_LABELS_10) +
  scale_colour_manual(values = GENE_SET_COLOURS) +
  xlab("Rescue probability bin") +
  ylab("Proportion of canonical splice site\nvariants observed in gnomAD") +
  theme_rr()
save_figure("figS19", fig_s19, 12, 6)

# %% [markdown]
# ## Figure N2: validation of rescue predictions
#
# A: heterozygous LOFTEE HC canonical splice PTVs in 1000 Genomes individuals with
# Geuvadis RNA sequencing. For each variant, transcript and affected exon, three
# junctions are counted: the canonical junction, the competitor junction (skipping
# the affected exon) and, where predicted, the cryptic junction. Sites where no
# junction has a read are left out (the transcript is not expressed), as are
# variants RescueRanger flags as UTR. Bars count junctions with at least one read,
# by read count (pseudo-log axis), for unrescued and rescued PTVs.
#
# B: *DDX3X* saturation genome editing SNVs in NMD target regions of
# ENST00000457138, with RescueRanger's predicted splice consequences. A predicted
# junction is in frame when its offset from the canonical site is a multiple of
# three. For each window (50, 100, 250 or 500 bp), a variant is predicted rescued
# when its highest in-frame alternative splice probability within that distance of
# the canonical site exceeds 0.2. The ROC curves compare these predictions with the
# day 15 combined log fold change.

# %%
## A: junction read support
junctions <- fread(PATHS$junctions)
junctions[, site_id := paste(chr, pos, ref, alt, transcript_id, exon_id, sep = "_")]
junctions <- junctions[site_id %in% junctions[n_reads >= 1, site_id] &
                         !grepl("3utr|5utr", flags)]
junctions[, rescue := factor(rescue, levels = c(0, 1),
                             labels = c("Unrescued PTVs", "Rescued PTVs"))]
junctions[, junction_type := factor(junction_type, levels = c("canonical", "competitor", "cryptic"),
                                    labels = c("Canonical", "Competitor", "Cryptic"))]

fig_n2a <- ggplot(junctions[n_reads > 0], aes(n_reads, fill = junction_type)) +
  geom_histogram(colour = "white") +
  facet_wrap(~rescue, ncol = 1) +
  scale_x_continuous(trans = "pseudo_log") +
  scale_fill_manual(values = c(Canonical = "maroon", Competitor = "olivedrab4",
                               Cryptic = "steelblue4")) +
  xlab("Number of reads supporting splice junction") +
  ylab("Number of splice junctions") +
  theme_rr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 13, face = "bold"))

# %%
## B: DDX3X saturation genome editing
ALT_THRESHOLD <- 0.2
WINDOWS <- c(50, 100, 250, 500)
sge <- fread(PATHS$ddx3x_sge)
## SNVs, in-frame predicted junctions
sge <- sge[ref_len - alt_len == 0 & splice_offset %% 3 == 0]
## per variant and splice site: rescue predicted within each window
sge <- rbindlist(lapply(WINDOWS, function(w)
  sge[abs(splice_offset) < w, .(prediction_label = max(alt > ALT_THRESHOLD)),
      by = .(variant, site, D15_combined_LFC)][, window_size := w]))
sge_roc <- sge[, {
  r <- roc(prediction_label, D15_combined_LFC, quiet = TRUE)
  .(fpr = 1 - r$specificities, tpr = r$sensitivities, auc = as.numeric(auc(r)))
}, by = window_size]
sge_roc[, window := sprintf("window=%d (AUC = %.2f)", window_size, auc)]
sge_roc[, window := factor(window, levels = unique(window[order(window_size)]))]
setorder(sge_roc, fpr, tpr)
print(unique(sge_roc[, .(window_size, auc)]))

fig_n2b <- ggplot(sge_roc, aes(fpr, tpr, colour = window, group = window)) +
  geom_point() +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  scale_colour_manual(values = c("#440154", "#31688E", "#35B779", "#FDE725")) +
  labs(x = "False positive rate (FPR)", y = "True positive rate (TPR)", colour = NULL) +
  theme_rr() +
  theme(legend.position = c(0.98, 0.02), legend.justification = c("right", "bottom"))

# %%
fig_n2 <- (fig_n2a | fig_n2b) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))
save_figure("figN2", fig_n2, 16, 7)

# %%
sessionInfo()

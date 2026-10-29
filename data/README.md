# Data

The notebooks read every file you download from `DATA_DIR` (default: this folder), and
the resources included with the code from `RESOURCE_DIR` (`config.R`;
`../resources/README.md`). Files with individual-level data are released under managed
access at the European Genome-phenome Archive (EGA). The aggregate counts, and the
non-cohort data generated in this work, are released openly as supplementary data with
the paper. Put both sets in this folder. No released file carries an individual
identifier.

## 1. Files, by notebook

| File | Contents | Released | 01 | 02 | 03 | 04 |
|---|---|---|:-:|:-:|:-:|:-:|
| `parental_variants.tsv.gz` | parental PTV and synonymous calls, one row per variant × transcript × parent | EGA | | ✓ | ✓ | ✓ |
| `proband_variants.tsv.gz` | proband PTV and synonymous calls, inherited or *de novo*, one row per variant × transcript | EGA | ✓ | ✓ | ✓ | |
| `burden_per_proband.tsv` | per proband: sHet burden and number of inherited PTVs (Figure 2E, Fig. S9) | EGA | | ✓ | | |
| `within_carrier_summary.tsv` | per carrier: inherited NMD-scaled PEXT and *de novo* outcomes (Figure 3D, Fig. S16) | EGA | | | ✓ | |
| `attributable_fraction_carriers.tsv`, `cohort_trios.tsv` | carrier probands by number of inherited PTVs, and trios per study (Figure 2D, Fig. S4) | supplementary data | | ✓ | | |
| `interpretation_counts.tsv` | carrier probands by clinical interpretation (Fig. S5) | supplementary data | | ✓ | | |
| `two_ptv_families.tsv` | families with two PTVs in DD-associated genes (Fig. S10) | supplementary data | | ✓ | | |
| `columns.tsv` | every column of the two variant tables | supplementary data | | | | |

Four supplementary data files, read by notebook 01, contain no cohort data:

| File | Contents | Figure |
|---|---|---|
| `splice_rescue_probabilities.tsv` | RescueRanger rescue probabilities for every possible canonical splice site SNV: `gene, chr, pos, ref, variant, var_pos, rescue` | Fig. S19 |
| `gnomad_ACs.tsv` | gnomAD v2.1.1 allele counts at those sites: `CHROM, POS, REF, ALT, AC, FILTER` | Fig. S19 |
| `junction_table.tsv` | Geuvadis RNA sequencing read counts for the canonical, competitor (exon skipping) and predicted cryptic junctions of heterozygous LOFTEE HC canonical splice PTVs in 1000 Genomes individuals, with the RescueRanger call: `chr, pos, ref, alt, transcript_id, exon_id, junction_type, n_reads, rescue, flags` | Figure N2A |
| `DDX3X_splice_consequences.tsv.gz` | RescueRanger's predicted splice consequences for *DDX3X* saturation genome editing SNVs in NMD target regions of ENST00000457138, with the day 15 combined log fold change: `variant, site, splice_offset, alt` (alternative splice probability), `ref_len, alt_len, D15_combined_LFC` | Figure N2B |

## 2. Processing applied before release

The cohort files were built from the cohort variant tables, the clinical records and genetic
ancestry while trio identifiers were still available (Methods).

- **Trios.** Trios with more than 20 filtered *de novo* calls are removed from every file.
- **Call quality.** Filters that every figure applies were applied, and their columns are
  not released.
  - Parents: VAF > 0.3, depth > 10 in parent and proband, QUAL > 300, heterozygous, and
    protein-coding transcripts. A transmitted record is kept only where the proband is
    heterozygous for the variant.
  - Probands: inherited or *de novo* (`is_denovo`), VAF > 0.3, heterozygous, outside
    ENCODE blacklist regions, QUAL ≥ 300 (SNVs) or ≥ 500 (indels), and depth > 10. DDD
    probands also need depth > 10 in both parents.
- **Flags.** Filters that only some figures apply are flags:
  - `strict_qc` (parents): outside ENCODE blacklist regions and, for DDD, depth > 10 in
    both parents. Used by the transmission analyses of notebook 02.
  - `parental_depth_ok` and `utr_flag` (probands): the *de novo* filters of notebook 03;
    `utr_flag` also applies to Fig. S20.
- ***De novo* calls.** `is_denovo` is TRUE for DeNovoGear posterior probability > 0.75
  (DDD), and for a DeNovoCNN *de novo* call with probability > 0.75 (RUMC).
- **Rescue.** `rescue` holds the RescueRanger prediction:
  - `rescued`
  - `not_rescued`
  - `not_rescued_new_ptv`: the original PTV is rescued, but the rescue creates another PTV
  - NA: no prediction
- **NMD.** `nmd_region` is the NMD region of the variant on that row's transcript: NMD
  target, start-proximal, long exon or distal (Blakes *et al.* 2026).
  - `nmd_scaled_pext` is PEXT multiplied by the region's NMD efficiency relative to NMD
    target regions (Methods).
  - Variants that cannot be placed are labelled `non_coding` (splice sites of untranslated
    exons) or `unplaced`, and have no NMD-scaled PEXT.
- **Genes.** Each row carries the gene symbol (`vepSYMBOL`) and Ensembl gene ID (`vepGene`),
  the gene's set, gnomAD v2.1.1 pLI, LOEUF and LOEUF decile, and sHet, as in
  `resources/gene_annotation.tsv.gz`. Constraint is matched by Ensembl gene ID, because a
  few symbols are on more than one row of the gnomAD table. Where the ID is not in that
  table, the symbol is used if it is on one row.
- **Child-level attributes.** Clinical data are released only for DD-associated genes, and
  only as far as the figures need them. The child's clinical status and genetic ancestry
  are attached to the records of every variant (table below).
- **No identifiers.** `obs_id` is a random label shared by the rows of one variant in one
  family: its transcript rows, and the mother's and father's records. Parental and proband
  records are labelled separately. Nothing links two variants of the same child.

| Column | Definition | Records | Used by |
|---|---|---|---|
| `clinical_interpretation` | most severe interpretation reported at the variant's position (P/LP, uncertain, B/LB) | DD-associated genes: transmitted unrescued HC PTVs (parental); *de novo* PTVs passing the *de novo* filters (proband) | Figs. S6, S18 |
| `proband_diagnosis_status` | diagnosed (DDD tier 4–5; RUMC report other than a VUS) / VUS / undiagnosed / not_available (no clinical record) | all parental | Fig. S8 groups; Fig. S13C outcome |
| `inherited_ptv_diagnosis` | the child carries a transmitted unrescued PTV in a DD-associated gene interpreted P/LP or uncertain | all parental | Fig. S8; Figure 2A and Fig. S2 triangles; Figs. S7, S13C |
| `proband_superpopulation_group` | European / Non-European (DDD: superpopulation; RUMC: largest ancestry component) / not_available (no ancestry call) | all parental and proband | Figures 2B and 2C (DDD-European sensitivity analysis, printed) |

## 3. Computed in the notebooks

- **One value per variant from its transcript rows.** Each notebook does this itself.
  - The analyses on NMD-scaled PEXT drop distal, `non_coding` and `unplaced` rows (rescued
    PTVs keep distal rows). Each variant then takes its highest NMD-scaled PEXT.
  - The analyses by NMD region place each variant in one region: NMD target first, then
    long exon, distal and start-proximal.
- **PEXT** without NMD scaling (Figs. S12, S13A) is NMD-scaled PEXT divided by the
  region's NMD efficiency.
- **Variants carried by both parents** are removed by `obs_id`, because their parent of
  origin is undefined.

## 4. Analyses that pool variants within an individual

These analyses need a quantity summed or counted over each child's (or family's) variants,
so they cannot run from variant records without identifiers. Each reads a file of
de-identified counts or summaries instead.

| Figure | What is pooled | File (columns) |
|---|---|---|
| **Figure 2E**, **Fig. S9** | sHet summed over each proband's transmitted unrescued PTVs (zero if none), and the number of those PTVs in constrained genes (LOEUF < 0.58); logistic models adjusted for study and ancestry, bootstrapped over probands | `burden_per_proband.tsv`: one row per proband, no identifiers: `study, ancestry_group, diagnosed, shet_burden, n_ptvs_constrained` |
| **Figure 2D**, **Fig. S4** | the number of PTVs each carrier proband inherited (a proband inheriting *k* PTVs counts as 1 − (1 − *q*)^*k*), and the number of probands per study | `attributable_fraction_carriers.tsv`: `study, ptvs, gene_set, k, n_probands`; `cohort_trios.tsv`: `study, n_trios` |
| **Fig. S5** | each carrier proband classified by the most severe interpretation among their inherited PTVs in DD-associated genes | `interpretation_counts.tsv`: `study, class, category, n_probands` |
| **Fig. S10** | families with exactly two unrescued parental PTVs in DD-associated genes: which parent carried them, how many were transmitted, whether they are on the same chromosome, and the child's status | `two_ptv_families.tsv`: `source, transmission, independent, dx_category, n`; the counts of the published figure |
| **Figure 3D**, **Fig. S16** | for each carrier, the sum and maximum NMD-scaled PEXT of their transmitted unrescued PTVs, with the same proband's *de novo* outcomes | `within_carrier_summary.tsv`: one row per carrier × gene set, no identifiers: `stratum, study, ancestry_group, index_pext_sum, index_pext_max, n_index, dnm_ptv_any, dnm_ptv_dd_any, dnm_synonymous_n, dnm_shet_burden` |

`burden_per_proband.tsv` and `within_carrier_summary.tsv` are still one row per person:
the regressions on continuous exposures cannot be rebuilt from aggregated counts. Their
rows are ordered by value rather than by proband.

## 5. Reproducing the figures

Every figure was checked on mock data, drawn once from the cohort tables with the code used
for the paper and once from release files built from the same tables. All are identical,
except:

- **Fig. S20** places variants with the released NMD regions.
- **Fig. S13C** leaves out the children flagged by `inherited_ptv_diagnosis`, the
  definition used for Figure 2A. Notebook 03's own list differed by two children in the
  cohort data.
- **Random draws.** These results depend on row order as well as the seed, and agree up to
  Monte Carlo error:
  - the bootstrap band of Figure 2E;
  - the dose-response bands of Figure 3C and Figure N3;
  - the permutation P values of Fig. S18.
- **Draw order.** Figs. S5, S6 and S14 draw the same bars in a different order, which
  changes a few pixels at bar edges.
- **Notebooks 03 and 04** keep about 1,800 parental variants (0.06%) that the code used
  for the paper left out. Their transcript rows differed only in PEXT values these analyses
  do not read (the original `PEXT` column, or PEXT on rows outside the NMD regions), which
  that code treated as two records.

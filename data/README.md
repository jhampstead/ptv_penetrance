# Data

The notebooks read every file you download from `DATA_DIR` (default: this folder), and
the resources included with the code from `RESOURCE_DIR` (`config.R`;
`../resources/README.md`). Files with individual-level data are released under managed
access at the European Genome-phenome Archive (EGA). The aggregate counts, and the
non-cohort data generated in this work, are released openly as supplementary data with
the paper. Put both sets in this folder.

## Data files

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

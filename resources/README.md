# Resources

Gene-level and other resources the notebooks read from `RESOURCE_DIR` (default: this
folder). All are included in the repository to reproduce figures and analyses from the manuscript. None contain individual-level data.
Data files associated with this manuscript, whether individual-level from EGA or aggregate from the paper's supplementary data, go in `data/`
(`../data/README.md`).

| Path | Contents | Used by |
|---|---|---|
| `gene_annotation.tsv.gz` | One row per row of the gnomAD v2.1.1 constraint table, plus the DDG2P and sHet genes it lacks: `gene, gene_id, gene_set, pLI, LOEUF, LOEUF_decile, sHet, ddg2p_mechanism`. | all |
| `roulette_ptv_grid.tsv.gz` | Expected *de novo* PTVs: Roulette mutation rates summed by gene and NMD-scaled PEXT (`GENE, PEXT_pct, MR, n_sites`). See below. | 03 |
| `samocha_2024_supp_tables.xlsx`, `samocha_2024_supp_table10.txt` | Samocha *et al.* (2024), *medRxiv*, [doi:10.1101/2024.08.28.24312746](https://doi.org/10.1101/2024.08.28.24312746): Supplementary Tables 2-8 and 12 (media-2.xlsx; sheets SuppTable2 and SuppTable4 are read) and Supplementary Table 10 (media-3.txt), unmodified apart from the file names, under the preprint's CC-BY 4.0 license | 02 |
| `gene_penetrance_consolidated.tsv` | Published evidence of reduced penetrance for the genes in Fig. S17, curated from the literature: `gene, classification, evidence, citation, reference_detail` | 04 |

## `gene_annotation.tsv.gz`

- **`gene_set`** (Methods). DD-associated genes are in the Developmental Disorder
  Genotype2Phenotype (DDG2P) panel, exported from
  [G2P](https://www.ebi.ac.uk/gene2phenotype) on 10 August 2026. A gene qualifies if it
  has all of the following:
  - a monoallelic autosomal entry with loss-of-function evidence (molecular mechanism "loss
    of function", or variant consequence "absent gene product" or "decreased gene product
    level");
  - no monoallelic entry with a gain-of-function, dominant-negative or undetermined
    non-loss-of-function mechanism, and no biallelic entry (disputed and refuted entries
    are excluded);
  - pLI ≥ 0.9.

  The pLI ≥ 0.9 set is the remaining genes with pLI ≥ 0.9, and pLI < 0.9 is all other
  genes.
- **`pLI`, `LOEUF` (`oe_lof_upper`) and `LOEUF_decile` (`oe_lof_upper_bin`)** are from
  gnomAD v2.1.1 (`gnomad.v2.1.1.lof_metrics.by_gene.txt`).
- **`sHet`** is the mean of the positive `s_het_det` values per gene in Weghorn *et al.*
  (2019), *Mol Biol Evol* 36:1701–1710.
- **`ddg2p_mechanism`** is the DDG2P mechanism class drawn in Figure 4A.

The variant tables carry the same gene-level values on every record. A few gene symbols are on
more than one row of the gnomAD table (different Ensembl genes), so constraint is matched by
Ensembl gene ID (`gene_id`, `vepGene` in the variant tables). Gene sets, sHet and the DDG2P
mechanism are matched by symbol.

## `roulette_ptv_grid.tsv.gz`

Every possible PTV SNV, from Roulette [(Seplyarskiy *et al.* 2023)](https://www.nature.com/articles/s41588-023-01562-0), is annotated with VEP and
PEXT and placed in the NMD regions of each of its transcripts [(Blakes *et al.* 2026)](https://www.nature.com/articles/s41467-026-69983-z), as
the cohort variants are. Rows in distal regions are removed. Each site's PEXT is then
multiplied by the NMD efficiency of its most efficient remaining region to generate
NMD-scaled PEXT.

Mutation rates (`MR`) and site counts (`n_sites`) are summed by gene and by NMD-scaled
PEXT and is rounded to a percentage (`PEXT_pct`).

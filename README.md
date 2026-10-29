# Incompletely penetrant protein-truncating variants confer developmental disorder risk through partial reductions in gene expression

Analysis code for the figures of Hampstead *et al.*, *medRxiv* (2026). The notebooks
regenerate every main-text figure and every supplementary figure from the variant-level
data. RescueRanger itself, which made the rescue predictions these analyses use, is
described in Note S1 and its Methods section.

## Layout

```
config.R            data paths, gene-set labels and colours, figure output
figures.tsv         every figure: key, label and output file name
notebooks/          one notebook per results section, as .ipynb (rendered on GitHub)
                    and as an identical plain R script (.R)
data/README.md      the released files the notebooks read, and how they were made
resources/          gene-level and other resources included with the code (see resources/README.md)
tools/              keeps each .ipynb in sync with its .R script
install.R           R packages
```

## Figures

| Notebook | Figures |
|---|---|
| `01_rescue` | Figure 1B; Figs. S1, S19, S20; Figure N2 |
| `02_transmission` | Figure 2 (A–E); Figs. S2–S10 |
| `03_pext` | Figure 3 (A–D); Figs. S11–S16, S18; Figure N3 |
| `04_burden_by_gene` | Figure 4 (A–C); Fig. S17 |

Figure 1A is a schematic. Supplementary tables and data files are not regenerated: the
notebooks contain only the code the figures need.

## Running

1. Install the packages: `Rscript install.R`.
2. Put every file you download in `data/` (`data/README.md`): the individual-level files
   from EGA (managed access) and the paper's supplementary data files. Every other input
   is included in `resources/` (`resources/README.md`). Set `RR_DATA_DIR`,
   `RR_RESOURCE_DIR` or `RR_FIG_DIR` to read files or write figures elsewhere.
3. Run a notebook, in any order:
   - **Jupyter:** open `notebooks/0X_*.ipynb` with the R kernel (`IRkernel::installspec()`)
     and run all cells.
   - **Command line:** `cd notebooks && Rscript 02_transmission.R`.

Figures are written to `figures/` as PNG (300 dpi) and SVG, named as in `figures.tsv`.

Random seeds are fixed in the code. Results drawn from random numbers (bootstrap and
Monte Carlo intervals, permutation P values) also depend on the order of the input rows,
so from the released files they agree with the paper up to Monte Carlo error
(`data/README.md`, section 5). `RR_BOOT` sets the number of bootstrap replicates for Fig. S12 (default 500). Figure 4B
sets gene names in italics with `ggtext` if it is installed, as published, and in plain
italics otherwise.

## Editing the notebooks

The `.R` files are the source; each `.ipynb` is generated from its `.R` file
([jupytext](https://jupytext.readthedocs.io) "percent" format, where `# %%` starts a code
cell and `# %% [markdown]` a text cell). After editing an `.R` file, regenerate the
notebooks:

```
pip install jupytext nbformat
python tools/sync_notebooks.py notebooks
```

## Data

Individual-level data are available through managed access at the European
Genome-phenome Archive (EGA): two variant tables without individual identifiers, and two
de-identified per-person summaries for the analyses that pool variants within a child.
Aggregate counts, the files describing the variant tables, and the non-cohort data
generated in this work (RescueRanger predictions, Geuvadis junction read counts, *DDX3X*
splice consequences) are released openly as supplementary data with the paper.
`data/README.md` describes each file, where it is released, which notebook reads it, and
the processing applied before release.

## Packages used by the notebooks (tested with R 4.3).
## Run once: Rscript install.R
## R.utils lets data.table read the .tsv.gz files directly
pkgs <- c("data.table", "R.utils", "ggplot2", "patchwork", "svglite", "scales", "readxl",
          "dplyr", "tidyr", "gtable", "ggrepel", "MASS", "pROC")
## optional: ggtext (italic gene labels in Figure 4B, as published),
## IRkernel + IRdisplay (to run the notebooks in Jupyter)
optional <- c("ggtext", "IRkernel", "IRdisplay")
install.packages(setdiff(c(pkgs, optional), rownames(installed.packages())))
## register the R kernel with Jupyter: IRkernel::installspec()

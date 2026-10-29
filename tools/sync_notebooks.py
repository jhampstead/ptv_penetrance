"""Regenerate notebooks/*.ipynb from the paired notebooks/*.R scripts.

The .R files (jupytext percent format) are the source. Each .ipynb gets the R
kernel and no outputs.

    python tools/sync_notebooks.py notebooks
"""
import pathlib
import sys

import jupytext
import nbformat

KERNEL = {"display_name": "R", "language": "R", "name": "ir"}
LANG = {"name": "R", "codemirror_mode": "r", "file_extension": ".r",
        "mimetype": "text/x-r-source", "pygments_lexer": "r"}

folder = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "notebooks")
for src in sorted(folder.glob("[0-9][0-9]_*.R")):
    nb = jupytext.read(src, fmt="R:percent")
    nb.metadata["kernelspec"] = KERNEL
    nb.metadata["language_info"] = LANG
    nb.metadata.pop("jupytext", None)
    for cell in nb.cells:
        cell.metadata = {}
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    nbformat.validate(nb)
    nbformat.write(nb, src.with_suffix(".ipynb"))
    print(f"{src.with_suffix('.ipynb').name}: {len(nb.cells)} cells")

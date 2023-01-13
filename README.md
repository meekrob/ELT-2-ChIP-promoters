# ELT-2-ChIP-revision

Analysis and data accompanying the following publication:
"Genome-wide characterization of the Caenorhabditis elegans intestine GATA transcription factor ELT-2" by Williams et al.

The analysis was performed using Rstudio and Rmarkdown documents. HTML renderings of these analyses are provided.

# Repository structure

- `\DATA`: contains input data for ELT-2 ChIP-seq aquired for the modERN database and FACS-isolated intestine RNA-seq data
- `\David`: contains analysis related to ELT-2 ChIP-seq, including peak-to-gene assignment
- `\Rob`: contains analysis related to FACS-isolated intestine RNA-seq, and ELT-2 target response to elt-2(-)

Each analysis directory contains three sub-directories: `01_input`, `02_scripts`, and `03_output`. In `01_input` additional data in stored necessary for the specific analysis. `02_scripts` contains the Rmarkdown files used to perform the analysis. `03_output` contains and files that were generated from the analysis including plots and text files.

# Usage

To perform the analysis in this repository, clone the repository. Ensure you have installed R (version 4.1.0) and Rstudio (Version 1.4.1106). Open the relevant `.Rproj` file. Package installation commands for the relevant Rmarkdown analysis file are included in the first code chunk.
# ibrutinib_swath.R
A custom R script for automatic data preprocessing, analysis and visualization of ibrutinib_swath dataset

## How to start
1. Download and install R (3.4.4 or later) and RStudio (1.1.453 or later)
2. Open RStudio
3. Install package dependencies 

```
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "ggrepel", "reshape2", "FactoMineR", "pheatmap"))

source("https://bioconductor.org/biocLite.R")

biocLite(c("biomaRt", "preprocessCore"))
```

4. Download the dataset (ibrutinib_swath.xlsx) onto the user's desktop. The dataset is available via ProteomeXchange repository (PXD013402) and also downloadable as the supplementary dataset 1 once the dataset published
5. Download ibrutinib_swath_final.R (or R_Notebook_ibrutinib_swath.Rmd) and open in RStudio
6. Run the code.

## Tutorial
A tutorial based on the R_Notebook_ibrutinib_swath.Rmd is available - see the Tutorial.pdf file for details

## Author
Somchai Chutipongtanate

## License
schuti/ibrutinib_swath.R is licensed under the GNU General Public License v3.0 - see the LICENSE file for details







# collagen_pca
Expression, biological and clinical relevance of the collagen pathway genes in prostate carcinoma

## Summary and results

To see the analysis report, please follow the links:

* _analysis report_: still pending

* PDF files with [analysis report figures](https://github.com/PiotrTymoszuk/collagen_pca/tree/main/report/figures)

* PDF files with [analysis report supplementary figures](https://github.com/PiotrTymoszuk/collagen_pca/tree/main/report/supplementary%20figures)

* an Excel file with [analysis report supplementary tables](https://github.com/PiotrTymoszuk/collagen_pca/blob/main/report/manuscript_supplementary_tables.xlsx)

To access manuscript files, please use the following links:

* a Word file with [main manuscript tables, figures and parts of methods](https://github.com/PiotrTymoszuk/collagen_pca/blob/main/report/report_supplementary_tables.xlsx)

* a Word file with _manuscript supplementary material_: still pending

* PDF files with [main manuscript figures](https://github.com/PiotrTymoszuk/collagen_pca/tree/main/report/manuscript%20figures)

* PDF files with _supplementary manuscript figures_: still pending

* an Excel file with [manuscript supplementary tables](https://github.com/PiotrTymoszuk/collagen_pca/blob/main/report/manuscript_supplementary_tables.xlsx)

## Terms of use
Publicly available data sets were analyzed. 
To reference and use analysis results, please cite our GitHub repository and our publication if available. In any questions, please contact [Isabel Heidegger](mailto:Isabel-Maria.Heidegger@i-med.ac.at) or [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Usage

To make sure to install required development packages prior to runing the pipeline:

```r

devtools::install_github('PiotrTymoszuk/ExDA')
devtools::install_github('PiotrTymoszuk/kmOptimizer')
devtools::install_github('PiotrTymoszuk/coxExtensions')
devtools::install_github('PiotrTymoszuk/trafo')
devtools::install_github('PiotrTymoszuk/clustTools')
devtools::install_github('PiotrTymoszuk/soucer')
devtools::install_github('PiotrTymoszuk/figur')
devtools::install_github('PiotrTymoszuk/caretExtra')
devtools::install_github('PiotrTymoszuk/xena')
devtools::install_github('PiotrTymoszuk/gseaTools')

```
To launch the entire pipeline, source the `exec.R` file:

```r

source('exec.R')

```

## Contact

The maintainer of the repository is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

# collagen_pca
Expression, biological and clinical relevance of the collagen pathway genes in prostate carcinoma

## Summary and results

To see the analysis report, please follow the links:

* _analysis report_: still pending

* [analysis report figures](https://github.com/PiotrTymoszuk/collagen_pca/tree/main/report/figures)

* [analysis report supplementary figures](https://github.com/PiotrTymoszuk/collagen_pca/tree/main/report/supplementary%20figures)

To access manuscript files, please use the following links:

* a Word file with _main manuscript tables, figures and parts of methods_

* a Word file with _manuscript supplementary material_

* PDF files with the _main manuscript figures_

* PDF files with the _supplementary manuscript figures_

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

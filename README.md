

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Quantifying Levallois: a 3D geometric morphometrics approach to Nubian technology

[![Last-changedate](https://img.shields.io/badge/last%20change-2024--09--03-brightgreen.svg)](https://github.com/jmcascalheira/LGMIberiaCluster/commits/master)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.2.4-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)

This repository contains the data and code for our paper:

> Hallinan, E. & Cascalheira, J. (2024). *Quantifying Levallois: a 3D
> geometric morphometrics approach to Nubian technology*.

### How to cite

Please cite this compendium as:

> Hallinan, E. & Cascalheira, J. (2024). *Compendium of R code and data
> for Quantifying Levallois: a 3D geometric morphometrics approach to
> Nubian technology*. Accessed 03 Sep 2024. Online at
> <https://doi.org/xxx/xxx>

## Contents

The **analysis** directory contains:

- [:file_folder: scripts](./analysis/scripts): Series of scripts with
  code to reproduce the analysis. It also has an annotated R Markdown
  version with a HTML rendered version suitable for reading the complete
  methods.
- [:file_folder: data](./analysis/data): Data used in the analysis.
- [:file_folder:
  supplementary-materials](./analysis/supplementary-materials):
  Supplementary materials including notes and other documents prepared
  and collected during the analysis.

## How to download and run locally

This research compendium has been developed using the statistical
programming language R. To work with the compendium, you will need
installed on your computer the [R
software](https://cloud.r-project.org/) itself and [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

You can download the compendium as a zip from from this URL:
[master.zip](./archive/master.zip). After unzipping: - open the `.Rproj`
file in RStudio - run `devtools::install()` to ensure you have the
packages this analysis depends on (also listed in the
[DESCRIPTION](./DESCRIPTION) file). - finally, go to [:file_folder:
scripts](./analysis/scripts), open each of the scripts in order and run
the code to produce the results. Alternatively, look for the
`/analysis/scripts/GEM_Nubian_methods.Rmd` and render it. ATTENTION:
both options will take several minutes.

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.

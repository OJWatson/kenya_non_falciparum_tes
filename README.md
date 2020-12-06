
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.0.3-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)

## Research compendium for non-falciparum in Kenya

This is a working R compendium (think R package but for reproducible
analysis). The analysis directory contains R scripts used to generate
the results.

### Installation

    git clone https://github.com/OJWatson/kenya_non_falciparum.git
    cd kenya_non_falciparum
    open kenya_non_falciparum.Rproj
    devtools::install_deps()

### Overview

The structure within analysis is as follows:

    |
    |──R/                      # R functions used throughout analysis
    |
    |──analysis/
       |
       ├── 00_<>.R             # analysis scripts
       |
       ├── data/
       │   ├── DO-NOT-EDIT-ANY-FILES-IN-HERE-BY-HAND
       │   ├── raw/            # data obtained from elsewhere
       │   └── derived/        # data generated during the analysis
       |
       ├── plots/              # plots generated by analysis
       |
       ├── tables/             # table generated by analysis

### Compendium DOI:

<http://doi.org/10.5281/zenodo.4308404>

The files at the URL above will generate the results as found in the
publication. The files hosted at
github.com/OJWatson/kenya\_non\_falciparum are the development versions
and may have changed since the report was published

### Overview of contents

This repository is our research compendium for our analysis of
<https://github.com/OJWatson/kenya_non_falciparum>. The compendium
contains all data, code, and text associated with the publication. The
`R` files in the `analysis` directory contain details of how all the
analyses reported in the paper were conducted, as well as instructions
on how to rerun the analysis to reproduce the results. The `data/`
directory in the `analysis/` directory contains all the raw and derived
data generated.

### The R package

This repository is organized as an R package. There are only a few R
functions exported in this package - the majority of the R code is in
the analysis directory. The R package structure is here to help manage
dependencies, to take advantage of continuous integration, and so we can
keep file and data management simple.

To download the package source as you see it on GitHub, for offline
browsing, use this line at the shell prompt (assuming you have Git
installed on your computer):

``` r
git clone https://github.com/OJWatson/kenya_non_falciparum.git
```

Once the download is complete, open the `kenya_non_falciparum.Rproj` in
RStudio to begin working with the package and compendium files. We will
endeavour to keep all package dependencies required listed in the
DESCRIPTION. This has the advantage of allowing
`devtools::install_dev_deps()` to install the required R packages needed
to run the code in this repository

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2020, copyright
holder: OJ Watson

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

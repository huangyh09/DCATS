[![Linux Build Status](https://travis-ci.org/huangyh09/DCATS.svg?branch=master)](https://travis-ci.org/huangyh09/DCATS)
[![codecov.io](https://codecov.io/github/huangyh09/DCATS/coverage.svg?branch=master)](https://codecov.io/github/huangyh09/DCATS/?branch=master)

# DCATS: Differential Composition Analysis Transformed by a Similarity matrix

This R package contains methods to detect the differential composition 
abundances between two conditions in singel-cell RNA-seq experiments, with or 
without replicates.

## Installation

### From R

The **latest** `DCATS` package can be conveniently installed using the
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) package thus:

```{R}
devtools::install_github("huangyh09/DCATS", build_vignettes = TRUE)
```

#### For development

Download this repository to your local machine and open it in Rstudio as a 
project, and build it by install and restart.


## Getting started

The best place to start are the vignettes. From inside an R session, load 
`DCATS` and then browse the vignettes:

```{r}
library(DCATS)
browseVignettes("DCATS")
```




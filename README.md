## What is `twc`?

An open-access tool/framework that constitutes the core functions to analyze terrestrial water cycle data across various spatio-temporal scales.

## Features

* Homogenized database currently encompassing 27 different precipitation data products.
* Exploratory data analysis.
* Spatiotemporal analysis.
* Time series analysis.
* Simple and aesthetic data visualization.

## Installation

```r
# Get the latest version of pRecipe

if (!require('devtools')) {install.packages('devtools'); library(devtools)}

install_github('imarkonis/twc')

# Alternatively, the latest CRAN release

install.packages("twc")

#

library(twc)
?`twc-package`
```

## References

Vargas Godoy, M.R. and Markonis, Y., 2023. pRecipe: A global precipitation climatology toolbox and database. Environmental Modelling & Software, 165, p.105711.

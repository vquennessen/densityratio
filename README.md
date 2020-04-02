
<!-- README.md is generated from README.Rmd. Please edit that file -->

# densityratio

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/vquennessen/densityratio.svg?branch=master)](https://travis-ci.org/vquennessen/densityratio)
<!-- badges: end -->

The goal of densityratio is to assess the performance of density ratio
control rules (DRCR) in terms of biomass, yield, spawning stock
abundance, age structure, or some other biological reference. Different
DRCRs may sample different individuals (all vs. only mature
individuals), different areas (all areas vs. only areas far from the
reserve), different amounts of time (1 year or multiple), may have a
different final target density ratio, and may use a static or a
transient final target density ratio. Alternatively, it can be used to
calculate simpler life history traits for a given species, such as
length at age, weight at age, proportion mature at age, selectivity at
age, the stable age distribution, and more.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vquennessen/densityratio")
```

## Example 1: Length at Age of Black Rockfish (2015)

This is a basic example which shows you how to solve a common problem:

``` r
library(densityratio)
## basic example code
```

## Example 2: Relative yield after reserve implementation given static and transient control rules (Black Rockfish 2015)

This is a basic example which shows you how to solve a common problem:

``` r
library(densityratio)
## basic example code
```

Don’t forget to commit and push the resulting figure files, so they
display on GitHub\!


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
remotes::install_github("vquennessen/densityratio")
#> Downloading GitHub repo vquennessen/densityratio@master
#>   
  
  
   checking for file 'C:\Users\Vic\AppData\Local\Temp\Rtmpc3iCB0\remotes76045baa5c53\vquennessen-densityratio-0838d53/DESCRIPTION' ...
  
   checking for file 'C:\Users\Vic\AppData\Local\Temp\Rtmpc3iCB0\remotes76045baa5c53\vquennessen-densityratio-0838d53/DESCRIPTION' ... 
  
v  checking for file 'C:\Users\Vic\AppData\Local\Temp\Rtmpc3iCB0\remotes76045baa5c53\vquennessen-densityratio-0838d53/DESCRIPTION' (416ms)
#> 
  
  
  
-  preparing 'densityratio':
#> 
  
   checking DESCRIPTION meta-information ...
  
   checking DESCRIPTION meta-information ... 
  
v  checking DESCRIPTION meta-information
#> 
  
  
  
-  checking for LF line-endings in source and make files and shell scripts
#> 
  
  
  
-  checking for empty or unneeded directories
#> 
  
  
  
-  building 'densityratio_2.2.0.tar.gz'
#> 
  
   
#> 
#> Installing package into 'C:/Users/Vic/Documents/R/win-library/3.5'
#> (as 'lib' is unspecified)
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

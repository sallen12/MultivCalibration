# MultivCalibration

[![DOI](https://zenodo.org/badge/596475593.svg)](https://zenodo.org/doi/10.5281/zenodo.10201288)

This repository contains R code to reproduce the results presented in the preprint  

> Allen, S., Ziegel, J. and Ginsbourger, D. (2023). 
> Assessing the calibration of multivariate probabilistic forecasts.
> ArXiv preprint.
> [arxiv.2307.05846](https://arxiv.org/abs/2307.05846)

Scripts to reproduce the results therein can be found in the _scripts_ repository. 

## Forecast calibration

Loosely speaking, probabilistic forecasts are _calibrated_ if they align statistically with the corresponding outcomes. Calibration is a necessary property for forecasts to be considered trustworthy.

Several notions of calibration exist for univariate forecasts, and univariate calibration is often assessed in practice using rank histograms (more generally probability integral transform (PIT) histograms). This paper discusses how multivariate forecasts and observations can be transformed to univariate objects, from which univariate rank histograms can be constructed. There is considerable flexibility in how to transform the multivariate forecasts and observations, and several transformations can be employed to gain a more complete understanding of how the multivariate forecasts behave.

This repository provides the functionality to implement these multivariate rank histograms in practice. Functionality is currently available for transformations that have previously been proposed in the literature, as well as custom user-specified transformations. The vignette thoroughly documents the usage of the repository, and reproduces the results in the above paper.

## Installation and development

This package has not been submitted to CRAN, and can therefore be installed in R using devtools
```r
# install.packages("devtools")
library(devtools)
install_github("sallen12/MultivCalibration")
```
The package is still in active development, and the vignette lists several possible extensions that could be implemented. Comments, suggestions, and input are more than welcome.

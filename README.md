[![R-CMD-check](https://github.com/shiyangm/LAVA-Knock/workflows/R-CMD-check/badge.svg)](https://github.com/shiyangm/LAVA-Knock/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# LAVA-Knock 
This is an R package of performing local genetic correlation via knockoffs. LAVA-Knock is a local genetic correlation method that builds off an existing genetic correlation method, LAVA, and augments it by leveraging synthetic/knockoff summary statistics.

## Description
The package contain functions for knockoff generation of LAVA-Knock univariate and bivariate testing, as well as the knockoff filter to detect significant local genetic correlation.

## Prerequisites
R (recommended version 4.3.1)

## Dependencies
LAVAKnock depends on R packages SKAT, Matrix, MASS, SPAtest, CompQuadForm, irlba, matrixsampling, corpcor and GhostKnockoff. Make sure to install those packages before installing BIGKnock.

## Installation
library(devtools) 

devtools::install_github("shiyangm/LAVA-Knock")

The current version is 0.1 (May 21, 2024).

## Usage
Please see the LAVA-Knock <a href="https://github.com/shiyangm/LAVA-Knock/blob/master/LAVAKnock_0.1.pdf"> **user manual** </a> for detailed usage of LAVA-Knock package. 

## Contact
If you have any questions about LAVA-Knock please contact

- <mashiyang1991@outlook.com>

If you want to submit a issue concerning the software please do so using the **LAVA-Knock** [Github repository](https://github.com/shiyangm/LAVA-Knock/issues).


## Reference
* Ma, S., Wang, F., Border, R., Buxbaum, J., Zaitlen, N. and Ionita-Laza, I. (2024+) Local genetic correlation via knockoffs reduces confounding due to cross-trait assortative mating.

## License
This software is licensed under GPLv3.

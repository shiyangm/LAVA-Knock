% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LAVAKnock.R
\docType{data}
\name{Example.LAVAKnock_bivariate}
\alias{Example.LAVAKnock_bivariate}
\alias{LAVAKnock_bivariate.example}
\title{Data example for LAVA-Knock bivariate local genetic correlation analysis on one window.}
\format{
An object of class \code{list} of length 4.
}
\usage{
data("LAVAKnock_bivariate.example")
}
\description{
This example dataset contains a 100-Kb window on chr1:64,344,227-64,444,227, its genotype matrix G_window with n=20000 samples and p=10 variants, with the original and M=5 knockoff Z-scores of 10 variants and 2 phenotypes
}
\examples{
data("LAVAKnock_bivariate.example")
window=LAVAKnock_bivariate.example$window
window
G_window=LAVAKnock_bivariate.example$G_window
dim(G_window) #20000    10
Zscore_pheno1_window=LAVAKnock_bivariate.example$Zscore_pheno1_window
dim(Zscore_pheno1_window) #10  6
Zscore_pheno2_window=LAVAKnock_bivariate.example$Zscore_pheno2_window
}
\keyword{data}

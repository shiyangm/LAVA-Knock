% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LAVAKnock.R
\name{LAVAKnock}
\alias{LAVAKnock}
\title{Conduct LAVA-Knock knockoff filter for multiple windows.}
\usage{
LAVAKnock(
  M = 5,
  p0 = Results_LAVAKnock.example$pval.orginal,
  p_ko = cbind(Results_LAVAKnock.example$pval.knockoff1,
    Results_LAVAKnock.example$pval.knockoff2, Results_LAVAKnock.example$pval.knockoff3,
    Results_LAVAKnock.example$pval.knockoff4, Results_LAVAKnock.example$pval.knockoff5),
  fdr = 0.1,
  window_id = Results_LAVAKnock.example$window,
  Rej.Bound = 20000
)
}
\arguments{
\item{M}{number of knockoffs. The recommended number is 5.}

\item{p0}{p-values of the original Z-scores of multiple windows, which is a w*1 vector for w windows.}

\item{p_ko}{matrix of the p-values for M knockoffs, which is a w*M vector for w windows.}

\item{fdr}{target FDR level. The recommended level is 0.1.}

\item{window_id}{id of the windows considered for multiple tesing.}

\item{Rej.Bound}{calculate ratios for top Rej.Bound tau values. The recommended level is 20000.}
}
\value{
\item{W}{knockoff statistics for each window.}

\item{Qvaluel}{Qvalue for each window.}

\item{W.threshold}{threshold of the W statistics with target fdr level.}

\item{window_sign}{Significant windows with q-values less then the fdr threshold.}
}
\description{
This function conducts the knockoff filter of LAVA-Knock on multiple windows and select significant windows.
}
\examples{
Results_LAVAKnock=LAVAKnock(M=5,p0=Results_LAVAKnock.example$pval.orginal,
p_ko=cbind(Results_LAVAKnock.example$pval.knockoff1,
          Results_LAVAKnock.example$pval.knockoff2,
          Results_LAVAKnock.example$pval.knockoff3,
          Results_LAVAKnock.example$pval.knockoff4,
          Results_LAVAKnock.example$pval.knockoff5),
fdr = 0.1,window_id=Results_LAVAKnock.example$window,Rej.Bound=20000)
Results_LAVAKnock$W.threshold#1.071494
sum(Results_LAVAKnock$W>=Results_LAVAKnock$W.threshold) #278
sum(Results_LAVAKnock$Qvalue<=0.1) #278
print(paste0(length(Results_LAVAKnock$window_sign), ' detected windows under M=5 knockoffs'))
#278 detected windows under M=5 knockoffs

}

utils::globalVariables(c('Zscore_0','pf','X_window','Zscore_pheno1_window','Zscore_pheno2_window','pnorm','Results_LAVAKnock.example','median','Rej.Bound',"stats", "as.dist", "cutree", "hclust"))

#' Data example for LAVA-Knock univariate heritability test on one locus.
#'
#' This example dataset contains a locus on chr1:64,344,227-65,894,184, its genotype matrix G_locus with n=20000 samples and p=197 variants, with the original Z-score of 197 variants and 2 phenotypes
#' 
#' @name Example.LAVAKnock_univariate
#' @docType data
#' @keywords data
#' @usage data("LAVAKnock_univariate.example")
#' @examples
#' data("LAVAKnock_univariate.example")
#'
#'Zscore_0=LAVAKnock_univariate.example$Zscore_0
#'G_locus=LAVAKnock_univariate.example$G_locus
#'locus=LAVAKnock_univariate.example$locus
#'print(locus)
#'dim(G_locus) #20000  197
#'dim(Zscore_0) #197   2
"LAVAKnock_univariate.example"

#' Conduct loci-level univariate heritability test.
#'
#' This function conducts the loci-level univariate heritability test on one locus.
#'
#' @param chr chromosome.
#' @param locus_start locus start.
#' @param locus_end locus end.
#' @param G_locus The genotype matrix or the reference genotype for a locus, which is a n*p matrix.
#' @param Zscore The original Z-score of a locus, which is a p*2 matrix for p variants and 2 phenotypes.
#' @param n  sample size.
#' @param prune.thresh Pruning threshold of sungular value decomposition. The recommended level is to explain 99 percent of total variance.
#' @return \item{h2_phe1.orginal}{heritability estimate of phenotype 1.}
#' @return \item{h2_phe2.orginal}{heritability estimate of phenotype 2.}
#' @return \item{p_uni_phe1.orginal}{p-value of univariate heritability test of phenotype 1.}
#' @return \item{p_uni_phe2.orginal}{p-value of univariate heritability test of phenotype 2.}
#' @examples
#'data("LAVAKnock_univariate.example")
#'Zscore_0=LAVAKnock_univariate.example$Zscore_0
#'G_locus=LAVAKnock_univariate.example$G_locus
#'locus=LAVAKnock_univariate.example$locus
#'chr=locus$chr
#'locus_start=locus$start
#'locus_end=locus$stop
#'Results_univariate_locus=LAVAKnock_univariate(Zscore=Zscore_0,
#'G_locus=G_locus,chr,locus_start,locus_end,n=20000,prune.thresh=99)
#'Results_univariate_locus
#'
#' @import Matrix
#' @import SKAT
#' @import MASS
#' @import SPAtest
#' @import CompQuadForm
#' @import irlba
#' @export
LAVAKnock_univariate=function(Zscore=Zscore_0,G_locus=G_locus,chr,locus_start,locus_end,n=20000,prune.thresh=99){
  #missing genotype imputation
  G_locus[G_locus<0 | G_locus>2]<-0 #NA
  
  #sparse matrix operation
  MAF<-colMeans(G_locus)/2;MAC<-colSums(G_locus)
  s<-colMeans(G_locus^2)-colMeans(G_locus)^2
  SNP.index<-which(MAF>0 & MAC>=25 & s!=0 & !is.na(MAF))# & MAC>10
  G_locus<-G_locus[,SNP.index,drop=F]
  
  #get positions and reorder G_locus
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_locus)))
  G_locus<-G_locus[,order(pos),drop=F]
  
  G_locus<-Matrix(G_locus,sparse=T)
  #get positions
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_locus)))
  
  AF<-colMeans(G_locus)/2
  X_locus=(G_locus-matrix(2*AF, nrow=n, ncol=length(AF), byrow=TRUE))/matrix(sqrt(2*AF*(1-AF)), nrow=n, ncol=length(AF), byrow=TRUE)
  
  #conduct LAVA analyses
  svd_locus = svd(X_locus/sqrt(n-1))
  #eigenvectors and eigenvalues of LD matrix
  Q_locus = svd_locus$v 
  lambda_locus = svd_locus$d * svd_locus$d 
  
  cum.perc_locus = cumsum(lambda_locus / sum(lambda_locus) * 100)
  keep_locus = 1:min(which(cum.perc_locus >= prune.thresh))
  K_locus=length(keep_locus)
  
  Q_pruned_locus=Q_locus[,keep_locus]
  lambda_pruned_locus=lambda_locus[keep_locus] 
  
  r_locus = Zscore/ sqrt(Zscore^2 + n - 2)	# for continuous phenos, convert Z to r
  alpha_locus = Q_pruned_locus%*% diag(1/lambda_pruned_locus) %*% t(Q_pruned_locus) %*% r_locus #inverse of LD * beta from GWAS
  
  delta_locus = diag(c(sqrt(lambda_pruned_locus))) %*% t(Q_pruned_locus) %*% alpha_locus	# using keep to filter out dropped PCs (since this one is updated every time a PC is dropped)
  R2_locus = diag(t(r_locus) %*% alpha_locus) 
  eta2_locus = (n-1) / (n-K_locus-1) * (1-R2_locus) #residual variance
  sigma2_locus=eta2_locus/(n-1) #sampling variance
  h2_locus=1-eta2_locus  # get h2
  
  T_univariate =diag(t(delta_locus)%*%delta_locus)/sigma2_locus/K_locus
  p_univariate=pf(T_univariate, K_locus, n-K_locus-1, lower.tail=F)
  
  Results_univariate_locus<-data.frame(chr,locus_start,locus_end,
                                       h2_locus[1],h2_locus[2],
                                       p_univariate[1],p_univariate[2])
  
  colnames(Results_univariate_locus)=c('chr','locus.start','locus.end',
                                       'h2_phe1.orginal','h2_phe2.orginal','p_uni_phe1.orginal','p_uni_phe2.orginal')
  
  return(Results_univariate_locus)
}

#' Data example for LAVA-Knock bivariate local genetic correlation analysis on one window.
#'
#' This example dataset contains a 100-Kb window on chr1:64,344,227-64,444,227, its genotype matrix G_window with n=20000 samples and p=10 variants, with the original and M=5 knockoff Z-scores of 10 variants and 2 phenotypes
#' 
#' @name Example.LAVAKnock_bivariate
#' @docType data
#' @keywords data
#' @usage data("LAVAKnock_bivariate.example")
#' @examples
#' data("LAVAKnock_bivariate.example")
#' window=LAVAKnock_bivariate.example$window
#' window
#' G_window=LAVAKnock_bivariate.example$G_window
#' dim(G_window) #20000    10
#' Zscore_pheno1_window=LAVAKnock_bivariate.example$Zscore_pheno1_window
#' dim(Zscore_pheno1_window) #10  6
#' Zscore_pheno2_window=LAVAKnock_bivariate.example$Zscore_pheno2_window
"LAVAKnock_bivariate.example"

#' Conduct window-level LAVA-Knock bivariate local genetic correlation analysis.
#'
#' This function conducts the window-level bivariate local genetic correlation test on one window.
#'
#' @param chr chromosome.
#' @param window_start window start.
#' @param window_end window end.
#' @param G_window The genotype matrix for a window, which is a n*p matrix.
#' @param Zscore_pheno1_window The original and M=5 knockoffs Z-scores of a window under phenotype 1, which is a p*6 matrix for p variants.
#' @param Zscore_pheno2_window The original and M=5 knockoffs Z-scores of a window under phenotype 2, which is a p*6 matrix for p variants.
#' @param n  sample size.
#' @param prune.thresh Pruning threshold of sungular value decomposition. The recommended level is to explain 99 percent of total variance.
#' @return \item{rg.orginal}{local genetic correlation estimate of original Z-score.}
#' @return \item{pval.orginal}{bivariate p-value of local genetic correlation test for original Z-score.}
#' @return \item{rg.knockoff}{local genetic correlation estimate of knockoff Z-score.}
#' @return \item{pval.knockoff}{bivariate p-value of local genetic correlation test for knockoff Z-score.}
#' 
#' @examples
#'data("LAVAKnock_bivariate.example")
#'window=LAVAKnock_bivariate.example$window
#'G_window=LAVAKnock_bivariate.example$G_window
#'Zscore_pheno1_window=LAVAKnock_bivariate.example$Zscore_pheno1_window
#'Zscore_pheno2_window=LAVAKnock_bivariate.example$Zscore_pheno2_window
#'chr=window$chr
#'window_start=window$start
#'window_end=window$end
#'Results_LAVAKnock_bivariate_window=LAVAKnock_bivariate(Zscore_pheno1_window=Zscore_pheno1_window,
#'Zscore_pheno2_window=Zscore_pheno2_window,
#'G_window=G_window,chr=chr,window_start,window_end,n=20000,prune.thresh=99)
#'Results_LAVAKnock_bivariate_window
#'
#' @import Matrix
#' @import SKAT
#' @import MASS
#' @import SPAtest
#' @import CompQuadForm
#' @import irlba
#' @import matrixsampling
#' @export
LAVAKnock_bivariate=function(Zscore_pheno1_window=Zscore_pheno1_window,Zscore_pheno2_window=Zscore_pheno2_window,G_window=G_window,chr=chr,window_start,window_end,n=20000,prune.thresh=99){
  
  #missing genotype imputation
  G_window[G_window<0 | G_window>2]<-0 #NA
  
  #sparse matrix operation
  MAF<-colMeans(G_window)/2;MAC<-colSums(G_window)
  s<-colMeans(G_window^2)-colMeans(G_window)^2
  SNP.index<-which(MAF>0 & MAC>=25 & s!=0 & !is.na(MAF))# & MAC>10
  
  G_window<-G_window[,SNP.index,drop=F]
  
  #get positions and reorder G_window
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_window)))
  G_window<-G_window[,order(pos),drop=F]
  
  G_window<-Matrix(G_window,sparse=T)
  #get positions
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_window)))
  
  AF<-colMeans(G_window)/2
  X_window=(G_window-matrix(2*AF, nrow=n, ncol=length(AF), byrow=TRUE))/matrix(sqrt(2*AF*(1-AF)), nrow=n, ncol=length(AF), byrow=TRUE)
  dim(X_window)  #20000    10
  
  #conduct LAVA bivariate test on org Z-score and knockoff Z-score
  set.seed(12345+1)
  univ.bivariate_rg_window_org=univ.bivariate_rg_GhostKnockoff(X=X_window,Zscore=cbind(Zscore_pheno1_window$org,Zscore_pheno2_window$org),n=20000)
  univ.bivariate_rg_window_knockoff1=univ.bivariate_rg_GhostKnockoff(X=X_window,Zscore=cbind(Zscore_pheno1_window$knock1,Zscore_pheno2_window$knock1),n=20000)
  univ.bivariate_rg_window_knockoff2=univ.bivariate_rg_GhostKnockoff(X=X_window,Zscore=cbind(Zscore_pheno1_window$knock2,Zscore_pheno2_window$knock2),n=20000)
  univ.bivariate_rg_window_knockoff3=univ.bivariate_rg_GhostKnockoff(X=X_window,Zscore=cbind(Zscore_pheno1_window$knock3,Zscore_pheno2_window$knock3),n=20000)
  univ.bivariate_rg_window_knockoff4=univ.bivariate_rg_GhostKnockoff(X=X_window,Zscore=cbind(Zscore_pheno1_window$knock4,Zscore_pheno2_window$knock4),n=20000)
  univ.bivariate_rg_window_knockoff5=univ.bivariate_rg_GhostKnockoff(X=X_window,Zscore=cbind(Zscore_pheno1_window$knock5,Zscore_pheno2_window$knock5),n=20000)
  
  Results_bivariate_locus<-data.frame(chr,window_start,window_end,
                                      univ.bivariate_rg_window_org$rg,univ.bivariate_rg_window_org$bivar_p,
                                      univ.bivariate_rg_window_knockoff1$rg,univ.bivariate_rg_window_knockoff1$bivar_p,
                                      univ.bivariate_rg_window_knockoff2$rg,univ.bivariate_rg_window_knockoff2$bivar_p,
                                      univ.bivariate_rg_window_knockoff3$rg,univ.bivariate_rg_window_knockoff3$bivar_p,
                                      univ.bivariate_rg_window_knockoff4$rg,univ.bivariate_rg_window_knockoff4$bivar_p,
                                      univ.bivariate_rg_window_knockoff5$rg,univ.bivariate_rg_window_knockoff5$bivar_p)
  
  colnames(Results_bivariate_locus)=c('chr','window.start','window.end',
                                      'rg.orginal','pval.orginal',
                                      'rg.knockoff1','pval.knockoff1',
                                      'rg.knockoff2','pval.knockoff2',
                                      'rg.knockoff3','pval.knockoff3',
                                      'rg.knockoff4','pval.knockoff4',
                                      'rg.knockoff5','pval.knockoff5')
  
  return(Results_bivariate_locus)
}

#' Data example of LAVA-Knock knockoff filter for multiple testing.
#'
#' This example dataset contains the LAVAKnock bivariateresults on 1,468 windows. Each window has an original p-values and M=5 knockoff p-values
#' 
#' @name Example.LAVAKnock
#' @docType data
#' @keywords data
#' @usage data("Results_LAVAKnock.example")
#' @examples
#' data("Results_LAVAKnock.example")
#'dim(Results_LAVAKnock.example) #1468    7
"Results_LAVAKnock.example"

#' Conduct LAVA-Knock knockoff filter for multiple windows.
#'
#' This function conducts the knockoff filter of LAVA-Knock on multiple windows and select significant windows.
#'
#' @param M number of knockoffs. The recommended number is 5.
#' @param p0 p-values of the original Z-scores of multiple windows, which is a w*1 vector for w windows.
#' @param p_ko matrix of the p-values for M knockoffs, which is a w*M vector for w windows.
#' @param fdr target FDR level. The recommended level is 0.1.
#' @param window_id id of the windows considered for multiple tesing.
#' @param Rej.Bound calculate ratios for top Rej.Bound tau values. The recommended level is 20000.
#' @return \item{W}{knockoff statistics for each window.}
#' @return \item{Qvaluel}{Qvalue for each window.}
#' @return \item{W.threshold}{threshold of the W statistics with target fdr level.}
#' @return \item{window_sign}{Significant windows with q-values less then the fdr threshold.}
#' 
#' @examples
#'Results_LAVAKnock=LAVAKnock(M=5,p0=Results_LAVAKnock.example$pval.orginal,
#'p_ko=cbind(Results_LAVAKnock.example$pval.knockoff1,
#'           Results_LAVAKnock.example$pval.knockoff2,
#'           Results_LAVAKnock.example$pval.knockoff3,
#'           Results_LAVAKnock.example$pval.knockoff4,
#'           Results_LAVAKnock.example$pval.knockoff5),
#'fdr = 0.1,window_id=Results_LAVAKnock.example$window,Rej.Bound=20000)
#'Results_LAVAKnock$W.threshold#1.071494
#'sum(Results_LAVAKnock$W>=Results_LAVAKnock$W.threshold) #278
#'sum(Results_LAVAKnock$Qvalue<=0.1) #278
#'print(paste0(length(Results_LAVAKnock$window_sign), ' detected windows under M=5 knockoffs'))
#'#278 detected windows under M=5 knockoffs
#'
#' @export
LAVAKnock<-function(M=5,p0=Results_LAVAKnock.example$pval.orginal,
                    p_ko=cbind(Results_LAVAKnock.example$pval.knockoff1,
                               Results_LAVAKnock.example$pval.knockoff2,
                               Results_LAVAKnock.example$pval.knockoff3,
                               Results_LAVAKnock.example$pval.knockoff4,
                               Results_LAVAKnock.example$pval.knockoff5),
                    fdr = 0.1,window_id=Results_LAVAKnock.example$window,Rej.Bound=20000){
  
  T0=-log10(p0)
  T_ko=-log10(p_ko)
  T=cbind(T0,T_ko)
  T[is.na(T)]<-0
  
  if(M!=1){
    W=(T[,1]-apply(T[,2:(M+1)],1,median))*(T[,1]>=apply(T[,2:(M+1)],1,max))
    #we can try max-second largest in knockoff
    #W=(T[,1]-apply(T[,2:(M+1)],1,max))*(T[,1]>=apply(T[,2:(M+1)],1,max))
  }else{
    W=T[,1]-T[,2]
  }
  
  #let Inf-Inf=0
  W[is.na(W)]=0
  
  kappa=apply(T,1,which.max)-1 #max T is from original data (0) or knockoff data (1 to 5)
  tau=apply(T,1,max)-apply(T,1,function(x)median(x[-which.max(x)]))
  #tau=apply(T,1,max)-apply(T,1,function(x)max(x[-which.max(x)]))
  
  #b=order(tau,decreasing=TRUE) 
  b=order(-tau) #tau[b] #tau from largest to smallest
  c_0=kappa[b]==0  #only calculate q-value for kappa=0
  
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  
  #calculate q value for each gene/window
  qvalue=rep(NA,length(tau))
  for(i in 1:length(b)){
    qvalue[b[i]]=min(ratio[i:min(length(b),Rej.Bound)])*c_0[i]+1-c_0[i] #only calculate q-value for kappa=0, q-value for kappa!=0 is 1
    if(i>Rej.Bound){break}
  }
  
  #qvalue==MK.q.byStat(kappa,tau,M,Rej.Bound=20000)
  
  #W statistics threshold, which can not be 0
  W.threshold=MK.threshold.byStat(kappa,tau,M=M,fdr=fdr,Rej.Bound=Rej.Bound)
  
  #window_sign=as.character(window_id[which(round(qvalue,digits=6)<=fdr)])
  window_sign=window_id[which(W>=W.threshold)]
  return(list(W=W,W.threshold=W.threshold,Qvalue=qvalue,window_sign=window_sign)) #kappa=kappa,tau=tau,
}


#' Data example of LAVA-Knock Z-scores with Ghostknockoff generation.
#'
#' This example dataset contains the Z-scores for two phenotypes and the genotype matrix of a unit: chr1:58387484-70992905. Due to computational and memory size constaints, we use 100 snps as an example to generation knockoff Z-scores using Ghostknockoff.
#' 
#' @name Example.Ghostknockoff-generation
#' @docType data
#' @keywords data
#' @usage data("Ghostknockoff_generation.example")
#' @examples
#' data(Ghostknockoff_generation.example)
#' chr=Ghostknockoff_generation.example$chr
#' unit_start=Ghostknockoff_generation.example$unit_start
#' unit_end=Ghostknockoff_generation.example$unit_end
#' print(paste0('chr',chr,':',unit_start,'-',unit_end))
#' #"chr1:58387484-70992905"
#' Zscore_0_unit=Ghostknockoff_generation.example$Zscore_0_unit
#' dim(Zscore_0_unit)
#' G_unit=Ghostknockoff_generation.example$G_unit
#' dim(G_unit)
"Ghostknockoff_generation.example"

#' Conduct LAVA-Knock Ghostknockoff generation for a unit.
#'
#' This function generate the knockoff Z-scores for a unit using Ghostknockoff with (shrinkage) empirical LD matrix.
#'
#' @param Zscore_0_unit Zscores of a unit, a p*2 matrix for two phenotypes
#' @param G_unit The genotype matrix or the reference genotype for a unit, a n*p matrix
#' @param LD.threshold LD.threshold for single-linkage hierarchical clustering to filter out highly correlated variants. The recommended level is 0.75.
#' @param n sample size.
#' @return \item{Zscore_pheno1}{5 knockoff Zscores for phenotype 1.}
#' @return \item{Zscore_pheno2}{5 knockoff Zscores for phenotype 2.}
#' 
#' @examples
#'data(Ghostknockoff_generation.example)
#'Zscore_0_unit=Ghostknockoff_generation.example$Zscore_0_unit
#'dim(Zscore_0_unit)
#'G_unit=Ghostknockoff_generation.example$G_unit
#'dim(G_unit)
#'
#'Results_Ghostknockoff_generation=Ghostknockoff_generation(G_unit,
#'Zscore_0_unit,LD.threshold=0.75,n=20000)
#'dim(Results_Ghostknockoff_generation$Zscore_pheno1)
#'dim(Results_Ghostknockoff_generation$Zscore_pheno2)
#'
#'@import corpcor
#'@import Matrix
#'@import GhostKnockoff
#'@export
Ghostknockoff_generation<-function(G_unit,Zscore_0_unit,LD.threshold=0.75,n=20000){
  
  #missing genotype imputation
  G_unit[G_unit<0 | G_unit>2]<-0 #NA
  
  #sparse matrix operation
  MAF<-colMeans(G_unit)/2;MAC<-colSums(G_unit)
  s<-colMeans(G_unit^2)-colMeans(G_unit)^2
  SNP.index<-which(MAF>0 & MAC>=25 & s!=0 & !is.na(MAF))# & MAC>10
  G_unit<-G_unit[,SNP.index,drop=F]
  
  #get positions and reorder G_unit
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_unit)))
  G_unit<-G_unit[,order(pos),drop=F]
  
  G_unit<-Matrix(G_unit,sparse=T)
  #get positions
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_unit)))
  length(pos) #96
  
  #remove tightly linked variants
  cor.X<-sparse.cor(Matrix(G_unit))$cor
  Sigma.distance = as.dist(1 - abs(cor.X))
  fit = hclust(Sigma.distance, method="single")
  corr_max = LD.threshold
  clusters = cutree(fit, h=1-corr_max)
  cluster.index<-match(unique(clusters),clusters)
  G_unit<-G_unit[,cluster.index];MAF<-MAF[cluster.index];MAC<-MAC[cluster.index];pos<-pos[cluster.index]
  N.SNP<-ncol(G_unit)
  N.SNP #96
  
  G<-G_unit
  MAF<-apply(G,2,mean)/2
  MAC<-apply(G,2,sum)
  N.SNP<-ncol(G)
  G<-Matrix(G,sparse=T)
  
  #summary stats -- always use shrinkage
  cor.G<-matrix(as.numeric(corpcor::cor.shrink(G,verbose=F)), nrow=ncol(G))
  range(c(cor.G)[round(c(cor.G),digits = 2)!=1.00]) #-0.5351557  0.7372171 shrinkage LD for n=20k
  
  ###generate Z-score using Ghostknockoff
  print('generating knockoff Z-scores of locus using GhostKnockoff')
  
  dim(Zscore_0_unit) #96  2
  
  #fit null model
  set.seed(12345+1)
  max.size=dim(cor.G)[1]
  max.size #96
  #equivalent to sdp
  fit.prelim<-GhostKnockoff.prelim(cor.G,M=5,method='asdp',max.size=dim(cor.G)[1]) 
  
  n.study<-n
  #generate knockoff Z-scores for two phenotypes
  GK.stat_pheno1_unit<-GhostKnockoff.fit(as.matrix(Zscore_0_unit[,1]),n.study=n,
                                         fit.prelim=fit.prelim,gamma=1,weight.study=NULL)
  #knockoff Z-score
  Zscore_knockoff_pheno1=GK.stat_pheno1_unit$GK.Zscore_k
  dim(Zscore_knockoff_pheno1)
  range(Zscore_knockoff_pheno1[,1])
  
  GK.stat_pheno2_unit<-GhostKnockoff.fit(as.matrix(Zscore_0_unit[,2]),n.study=n,
                                         fit.prelim=fit.prelim,gamma=1,weight.study=NULL)
  #knockoff Z-score
  Zscore_knockoff_pheno2=GK.stat_pheno2_unit$GK.Zscore_k
  dim(Zscore_knockoff_pheno2)
  range(Zscore_knockoff_pheno2[,1])
  
  #save Z-score and knockoff Z-score of Ghostknockoff
  Zscore_pheno1=cbind(Zscore_0_unit[,1],Zscore_knockoff_pheno1)
  colnames(Zscore_pheno1)=c('org','knock1','knock2','knock3','knock4','knock5')
  
  Zscore_pheno2=cbind(Zscore_0_unit[,2],Zscore_knockoff_pheno2)
  colnames(Zscore_pheno2)=c('org','knock1','knock2','knock3','knock4','knock5')
  
  return(list(Zscore_pheno1=Zscore_pheno1,Zscore_pheno2=Zscore_pheno2))
  
}


######### Other functions #########
sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat)) 
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

#knockoff filter
MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=20000){
  #b<-order(tau,decreasing=T)
  b<-order(-tau)
  c_0<-kappa[b]==0
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  ok<-which(ratio<=fdr)
  if(length(ok)>0){
    #ok<-ok[which(ok-ok[1]:(ok[1]+length(ok)-1)<=0)]
    return(tau[b][ok[length(ok)]])  #tau cannot be 0
  }else{return(Inf)}
}

univ.bivariate_rg_GhostKnockoff=function(X=X_window,Zscore=cbind(Zscore_pheno1_window$org,Zscore_pheno2_window$org),n=20000,prune.thresh=99){
  
  svd_window = svd(X/sqrt(n-1))
  Q_window = svd_window$v 
  lambda_window = svd_window$d * svd_window$d 
  
  cum.perc_window = cumsum(lambda_window / sum(lambda_window) * 100)
  keep_window = 1:min(which(cum.perc_window >= prune.thresh))
  K_window=length(keep_window)
  Q_pruned_window=Q_window[,keep_window]
  lambda_pruned_window=lambda_window[keep_window] 
  
  r_window =Zscore/ sqrt(Zscore^2 + n - 2)	# for continuous phenos, convert Z to r
  alpha_window = Q_pruned_window%*% diag(1/lambda_pruned_window) %*% t(Q_pruned_window) %*% r_window #inverse of LD * beta from GWAS
  
  delta_window = diag(c(sqrt(lambda_pruned_window))) %*% t(Q_pruned_window) %*% alpha_window	# using keep to filter out dropped PCs (since this one is updated every time a PC is dropped)
  
  R2_window = diag(t(r_window) %*% alpha_window) 
  eta2_window = (n-1) / (n-K_window-1) * (1-R2_window) #residual variance
  sigma2_window=eta2_window/(n-1) #sampling variance
  h2_window=1-eta2_window  # get h2
  
  T_univariate =diag(t(delta_window)%*%delta_window)/sigma2_window/K_window
  p_univariate=pf(T_univariate, K_window, n-K_window-1, lower.tail=F)
  
  omega_window = t(delta_window)%*%delta_window/K_window - diag(sigma2_window)
  cov=omega_window[1,2]
  rg=omega_window[1,2]/sqrt(omega_window[1,1]*omega_window[2,2])
  if(sum(is.na(rg)|h2_window<0)>0|sum(diag(omega_window)<0)>=1){  
    rg=NA;bivar_p=NA
  }else{
    bivar_p=signif(integral.p(bivariate.integral, K = K_window, omega = omega_window, sigma =diag(c(sigma2_window)), adap.thresh=c(1e-4, 1e-6)), 6)
  }
  return(list(p_univariate=p_univariate,h2_window=h2_window,K_window=K_window,cov=cov,rg=rg,bivar_p=bivar_p)) 
}

##functions in LAVA package to compute bivariate p-values
# the integral.func argument should just be one of bivariate.integral, multivariate.integral or pcov.integral
# for pcov.integral, omega and sigma should be formatted such that the first two pheno are X and Y, and the remainder is Z
integral.p = function(integral.func, K, omega, sigma, min.iter=10000, adap.thresh=c(1e-4, 1e-6)) {
  tot.iter = min.iter * 10^(0:length(adap.thresh))
  adap.thresh = c(adap.thresh,0)  # adding dummy 0 at the end to simplify loop code
  
  p = 1; curr.iter = 0
  for (i in 1:length(tot.iter)) {
    add.iter = tot.iter[i] - curr.iter
    add.p = integral.func(K, omega, sigma, n.iter=add.iter)
    p = (curr.iter*p + add.iter*add.p) / tot.iter[i]
    curr.iter = tot.iter[i]
    if (all(is.na(p)) || all(p[!is.na(p)] >= adap.thresh[i])) break
  }
  return(p)
}

bivariate.integral = function(K, omega, sigma, n.iter=1000, add.reverse=T) {
  if (!add.reverse) {
    omega.null = diag(diag(omega))
    sig.use = matrix(0,3,3); sig.use[1,1] = sigma[1,1]
    theta = matrix(0,3,3); theta[-1,-1] = omega.null*K
    
    sig.xy = sigma[1,2]
    sig.xys = sig.xy/sigma[1,1]
    var.y = sigma[2,2] - (sigma[1,2]^2)/sigma[1,1]
    
    params = apply(matrixsampling::rwishart(n.iter, K, Sigma=sig.use, Theta=theta), 3, bivar.cond.stats, K=K, sig.xy, sig.xys, var.y)   # first row is means, second is SDs
    return(conditional.norm(omega[1,2], params[1,], params[2,]))
    
  } else {
    p1 = bivariate.integral(K, omega, sigma, n.iter/2, add.reverse=F)
    p2 = bivariate.integral(K, omega[2:1,2:1], sigma[2:1,2:1], n.iter/2, add.reverse=F)    
    return((p1+p2)/2)
  }
}

# this is an internal function for the apply in integral.p(), defined here for clarity
# draw will be the 3x3 matrix drawn from the wishart
bivar.cond.stats = function(draw, K, sig.xy, sig.xys, var.y) {
  m = draw[2,3] + draw[1,3] + sig.xys*(draw[1,2] + draw[1,1])
  m = m/K - sig.xy #mean
  
  v = var.y * (draw[2,2] + 2*draw[1,2] + draw[1,1])    
  v = v / K^2
  v = ifelse(v <= 0, NA, sqrt(v)) #sd
  return(c(m,v))  
}

conditional.norm = function(obs, means, sds) {
  obs = abs(obs)
  prob = suppressWarnings(pnorm(obs, mean=means, sd=sds, lower.tail=F))
  prob = prob + suppressWarnings(pnorm(-obs, mean=means, sd=sds, lower.tail=T))
  return(mean(prob, na.rm=T))
}

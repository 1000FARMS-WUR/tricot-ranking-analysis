
################################################################################
# Ranking data analysis of on-farme trials Model Functions
# 
# Authority: Wageningen university
# Author: Joost van Heerwaarden and Hugo Dorado
# Institution/Project:Wageningen university
# Date:24/10/2025
# Description: 
#   This script implements functions for fitting and evaluating Thurstonian 
#   model for tricot trials. Includes log-likelihood computation, optimization, 
#   ranking conversion, and estimation of variance components.
################################################################################

# Load required packages
require(tidyverse)
require(PlackettLuce)
require(gosset)
require(qvcalc)

#----------------------------------------------------------
# Function: contrast.diff.mat
# Purpose: Create a contrast matrix for the Thurstonian model
#----------------------------------------------------------
contrast.diff.mat <- function(len){ 
  if(len == 2){return(t(matrix(c(1,-1))))}
  
  c0 <- matrix(c(1, rep(0, len-2)), ncol=1)
  dM <- cbind(c0, matrix(0, ncol=len-1, nrow=len-1, byrow=TRUE))
  
  for(i in 2:(len-1)){
    cs <- dM[, i-1]
    cs[cs == -1] <- 0 
    cs <- cs * -1
    cs[i] <- 1
    dM[, i] <- cs
  }
  
  cf <- matrix(c(rep(0, len-2), -1), ncol=1)
  dM[, ncol(dM)] <- cf
  dM
}

#----------------------------------------------------------
# Function: thurstonian_loglikehood
# Purpose: Compute log-likelihood for the Thurstonian model
#----------------------------------------------------------
thurstonian_loglikehood <- function(alpha_bet){
  require(mvtnorm)
  ln <- length(alpha_bet)
  C <- contrast.diff.mat(ln)
  mu <- as.matrix(alpha_bet)
  
  Sm <- C %*% diag(length(alpha_bet)) %*% t(C)
  D <- if(nrow(Sm) == 1){Sm^(-1/2)} else {diag(diag(Sm)^(-1/2))}
  
  Smean <- D %*% C %*% mu
  Smf <- D %*% Sm %*% D
  pmv <- pmvnorm(lower=rep(-Inf, ln-1), upper=Smean[,1],
                 mean=rep(0, ln-1), sigma=Smf)
  
  log(pmv[1])
}

#----------------------------------------------------------
# Function: tricotNLLMU0
# Purpose: Compute log-likelihood for tricot data
#----------------------------------------------------------
tricotNLLMU0 <- function(data, par, model="t", ref){
  library(gtools)
  par <- c(ref=0, par)
  names(par)[1] <- paste0('mu_', ref)
  
  rankdata <- data$ranking
  matrixTricot <- as.data.frame(as.matrix(rankdata))
  colnames(matrixTricot) <- paste0('mu_', colnames(matrixTricot))
  
  muMat <- data.frame(nam.mu=names(par), par=par)
  matrixTricot <- cbind(Ind=row.names(matrixTricot), matrixTricot) %>%
    pivot_longer(cols = !Ind, names_to = "Item", values_to = "Ranking") %>%
    filter(Ranking != 0) %>%
    left_join(muMat, by=c('Item'='nam.mu')) %>%
    select(!Item) %>%
    pivot_wider(names_from = Ranking, values_from = par)
  
  matrixTricot <- matrixTricot[,-1]
  matrixTricot <- matrixTricot[, mixedorder(colnames(matrixTricot))]
  matrixTricot <- as.matrix(matrixTricot)
  
  loglk <- if(model == 't'){
    apply(matrixTricot, 1, thurstonian_loglikehood)
  } else {
    logworth <- exp(matrixTricot)
    apply(logworth, 1, pl_likehood)
  }
  
  -1 * sum(loglk * data$freq, na.rm=TRUE) # -1 to be minimized
}

#----------------------------------------------------------
# Function: thurstonianMod
# Purpose: Fit Thurstonian model (Case I) via optimization
#----------------------------------------------------------
thurstonianMod <- function(rnks, ref=NULL){
  mu.names <- dimnames(rnks)[[2]]
  ref <- if(is.null(ref)){mu.names[1]} else {ref}
  
  cp <- list(maxit=10000, temp=50, trace=TRUE, REPORT=5000)
  matMu <- runif(ncol(rnks), -5, 5)
  names(matMu) <- paste0('mu_', mu.names)
  
  rnks.freq <- aggregate(rnks)
  matMu <- matMu[-which(names(matMu) == paste0('mu_', ref))]
  
  opt <- optim(par=matMu, tricotNLLMU0, model="t", data=rnks.freq,
               method="BFGS", control=cp, hessian=TRUE, ref=ref)
  
  est_mu <- c(0, opt$par)
  mu_ref_nm <- paste0('mu_', ref)
  names(est_mu)[1] <- mu_ref_nm
  
  hess_mu <- cbind(muG1=0, rbind(muG1=0, opt$hessian))
  row.names(hess_mu)[1] <- mu_ref_nm
  colnames(hess_mu)[1] <- mu_ref_nm
  
  list(est_mu=est_mu, hess_mu=hess_mu, lkhd=opt)
}

#----------------------------------------------------------
# Function: convert_rankings
# Purpose: Convert tricot quantitative data into ranking format
#----------------------------------------------------------
convert_rankings <- function(data_simulated, names.genotype,
                             env='env.vec', farm='farm.vec', 
                             genotype='genotype', trait.value.rank='trait.value.rank'){
  require(PlackettLuce)
  
  data_simulated <- do.call(rbind, split(data_simulated[c(env,farm,genotype,trait.value.rank)],
                                         data_simulated[[farm]]))
  
  data_simulated_wide <- reshape(data_simulated, direction='wide',
                                 idvar=c(env, farm), timevar='genotype')
  
  names(data_simulated_wide) <- gsub(paste0(trait.value.rank, "\\."), '', names(data_simulated_wide))
  data_simulated_wide <- data_simulated_wide[, c(env, farm, names.genotype)]
  data_simulated_wide[is.na(data_simulated_wide)] <- 0
  
  ranks <- as.rankings(data_simulated_wide[,-c(1:2)])
  rownames(ranks) <- paste0(data_simulated_wide[[env]], ".", data_simulated_wide[[farm]])
  ranks
}

#----------------------------------------------------------
# Function: th_means_se
# Purpose: Compute Thurstonian model means and standard errors
#----------------------------------------------------------
th_means_se <- function(th_mod){
  require(qvcalc)
  hess_mat <- th_mod$lkhd$hessian
  inv_hess <- solve(hess_mat)
  means <- th_mod[[1]]
  names(means)[1] <- 'ref'
  
  vcov_mat <- cbind(ref=0, rbind(ref=0, inv_hess))
  qvcalc::qvcalc(vcov_mat, estimates=means)
}

#----------------------------------------------------------
# Function: qH2
# Purpose: Estimate broad-sense heritability (H²) from a mixed model
#----------------------------------------------------------
qH2 <- function(model){
  varcomp <- summary(model)$varcomp
  BLUPl <- summary(model, coef=TRUE)$coef.random
  BLUPl <- as.data.frame(BLUPl)
  BLUPl <- BLUPl[grep('genotype_', rownames(BLUPl)),]
  BLUPl$PEV <- BLUPl$std.error^2
  PEV1 <- BLUPl$PEV
  1 - mean(PEV1) / varcomp['genotype','component']
}


#----------------------------------------------------------
# Function: genotypic_means_table
# Purpose: Obtaining the genotypic means by enviroment using the thurstonian model
#----------------------------------------------------------


genotypic_means_table <- function(  
                                    data  ,                            # Dataset
                                    env   = 'environment',             # Column name for environment variable
                                    id    = 'registration_surveyid',   # Column name for farm or experimental unit ID
                                    items = 'genotype',                # Column name for genotype identifiers
                                    input = 'yield_ranked',            # Column name for ranked variable
                                    ref   = NULL                       # Reference genotype name
                                  
){
  if(is.null(ref)){
    stop("Error: 'ref' cannot be NULL. Please provide a valid reference object.")
  }
  dataset_ls <- split(data,data[env])
  
  dataset_rnks_ls <- lapply(dataset_ls,function(w){
    
    # Apply the Thurstonian model to each environment subset
    ranks=rank_numeric(data = w,id = id,items = items,input = input)
    th_mod_out <- thurstonianMod(rnks = ranks,ref = ref)
    th_means <- th_means_se(th_mod_out)$qvframe
    row.names(th_means)[1] <- ref
    th_means <- data.frame(genotype=row.names(th_means),th_means)
    row.names(th_means) <- NULL
    cbind(th_means,environment=unique(w[,env]),ref_gen=ref)
    
  })
  fnds <- do.call(rbind,dataset_rnks_ls)
  row.names(fnds) <- NULL
  fnds[,items] <- sub("^mu_", "", fnds[,items])
  fnds
}





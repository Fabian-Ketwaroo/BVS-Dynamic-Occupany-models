
# Rscript demonstrating BVS via PG works for DOM model

rm(list= ls())
library(nimble)
library(coda)
library(pgdraw)
library(data.table)
source("BVS_DOM_samplers.R") # load MCMC samplers

#  Simulated  data 
simdata <- function(M, J, T, ncov, true_betas_psi, true_betas_phi, true_betas_eta, true_betas_p){
  
  S = M*(T-1)
  N =  M*T*J
  
  # Design matrix for psi
  X = matrix(rnorm(M*(ncov-1) ), nrow = M, ncol = ncov-1)  # first 4 continouous and last one categorical with 3 levels
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  A_original = data.frame(x5 = as.factor(sample(x=c(1,2, 3), size=M, replace=TRUE, prob=rep(1/3, 3))) )
  A <- model.matrix(~ ., A_original) # creates a design (or model) matrix
  originalX =cbind(originalX, A_original)
  X_psi = cbind(X,  A[,-1] )
  #head(X_psi)
  #head(originalX)
  
  psi <-  ilogit(X_psi %*% true_betas_psi) # initial occupancy probability 
  mean(psi)
  
  # Design matrix for survival (phi) 
  X = matrix(rnorm(S*(ncov-1) ), nrow = S, ncol = ncov-1)  # first 4 continouous and last one categorical with 3 levels
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  A_original = data.frame(x5 = as.factor(sample(x=c(1,2, 3), size=S, replace=TRUE, prob=rep(1/3, 3))) )
  A <- model.matrix(~ ., A_original) # creates a design (or model) matrix
  head(originalX)
  originalX =cbind(originalX, A_original)
  X_phi = cbind(X,  A[,-1] )
  
  # Survival -here the 3rd and 5th beta are non-significant
  
  
  phi <-  matrix(ilogit(X_phi %*% true_betas_phi),M,T-1) # occupancy probability 
  mean(apply(phi, 2, mean))
  
  
  # Design matrix for colonization (eta) 
  X_eta = X_phi
  eta <-  matrix(ilogit(X_eta %*% true_betas_eta),M,T-1) # occupancy probability 
  mean(apply(eta, 2, mean))
  
  
  # Design matrix for p
  X = matrix(rnorm(N*(ncov-1) ), nrow = N, ncol = ncov-1)  # first 4 continouous and last one categorical with 3 levels
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  A_original = data.frame(x5 = as.factor(sample(x=c(1,2, 3), size=N, replace=TRUE, prob=rep(1/3, 3))) )
  A <- model.matrix(~ ., A_original) # creates a design (or model) matrix
  head(originalX)
  originalX =cbind(originalX, A_original)
  X_p = cbind(X,  A[,-1] )
  
  
  
  p <-  array(ilogit(X_p %*% true_betas_p), dim= c(M,J,T)) # occupancy probability
  mean(apply(p, 2, mean))
  
  z = matrix(0, M, T) # site specific occupancy
  Y =  array(0, dim = c(M,J,T))
  z[, 1] = rbinom(M, 1, psi)
  
  
  for (t in 2:T) {
    z[,t] = rbinom(M, 1, z[,t-1]*phi[,t-1] + (1-z[,t-1])*eta[, t-1])
  }
  
  for (t in 1:T) {
    Y[1:M, 1:J,t] = matrix(rbinom(M*J, 1, z[1:M,t]*p[1:M,1:J,t]), M,J)
    
  }
  
  y = c(Y) # data in scalar form
  
  return(list(y=y,Y= Y,X_p=X_p, X_phi= X_phi, X_eta= X_eta, X_psi = X_psi, originalX = originalX ))
  
}

M = 500 # sites 
J = 5 # sampling occasions
T = 10 # years

ncov = 5 # number of explanatory variables in the model excluding the intercept
# Here, the 1st 4 effects are continuous and the last effect is categorical with 3 levels
true_betas_psi <- c(0.75, 0.3, 0, 0.5, 0, 0.6, -0.4) # includes intercept
true_betas_phi <- c(1, 0.3, 0, 0.5, 0, -0.2, 0.9) # includes intercept
true_betas_eta <- c(-0.5, 0.3, 0, 0.5, 0, 0.35, -0.7) # includes intercept
true_betas_p <- c(1, 0.3, 0, 0.5, 0, -1, 0.33) # includes intercept

set.seed(1)
data = simdata(M, J, T, ncov, true_betas_psi, true_betas_phi, true_betas_eta, true_betas_p )
y=data$y; X_p=data$X_p; X_phi= data$X_phi; X_eta= data$X_eta; X_psi = data$X_psi; originalX = data$originalX

# Prerequisites needed to preform BVS
{
  column_covariate = 1:ncov
  
  # true or false if the covaraites is categorical or numerical 
  classCovariates <- c(rep(TRUE, ncov-1), FALSE) # first 4 continuous and last one categorical with 3 levels
  classCovariates
  
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1 # fix the intercept
  
  if(any(column_covariate != 0)){
    
    indexes_covariates <- c()
    indexes_covariates[1] <- 1
    k <- 2
    for (i in 1:ncov) {
      if(classCovariates[i]){
        indexes_covariates[k] <- i + 1
        k <- k + 1
      } else {
        num_levels <- length(unique(originalX[,i]))
        indexes_covariates[k + 0:(num_levels-2)] <- i + 1
        k <- k + num_levels - 1
      }
    }
    
  }
  
  indexes_covariates 
  
  
  # compute prior correlated (C) matrix 
  C <- matrix(0, nrow = ncol(X_p) - 1, ncol = ncol(X_p) - 1)
  l <- 0 # l loop through the covariates
  for (i in 1:ncov) {
    if(classCovariates[i]){ # if it's a numerical covariates
      C[l + 1, l + 1] <- 1
      l <- l + 1
    } else { # if it's a categorical covariates
      num_levels <- length(unique(originalX[,i]))
      C[l + 1:(num_levels-1), l + 1:(num_levels-1)] <- .5 * (diag(1, nrow = (num_levels-1)) +
                                                               matrix(1, nrow = (num_levels-1), ncol = (num_levels-1)))
      l <- l + num_levels - 1
    }
  }
  
  C
  Sigma <- cov(X_p[,-1, drop = F])
  sumCSigma <- sum(C * Sigma)
  
  mu_psi = 0.5
  phi_mu =4
  phi_beta = 1/4
  
  # Prior parameters for beta, beta ~ MVN(b, B)
  b_psi <- c(mu_psi, rep(0, ncol(X_p) - 1))
  sigma_beta_psi <- C * phi_beta
  b_psi #b
  
  B_psi <- matrix(0, nrow = ncol(X_p), ncol = ncol(X_p))
  B_psi[1, 1] <- phi_mu
  B_psi[2:(ncol(X_p)), 2:(ncol(X_p))] <- sigma_beta_psi
  B_psi # B
  
  numVars = dim(B_psi)[1] #number of covariates (including the intercept)
  
  nsamples = length(y)
  nsamples
  
  # indicator variable for z such that observation process can be scalar in NIMBLE model
  zindex = numeric(nsamples)
  zindex[1:(M*J)] = rep(1:M,J)
  
  for (t in 2:T) {
    zindex[((((t-1) *M*J)+1)): (M*J*t) ] = rep( ((t-1)*M+1):(t*M),J)
  }
  
  # shell distributions for PG
  dPG <- nimbleFunction(
    run = function(x = double(1), log=integer(0, default=1)){
      returnType(double(0))
      ans <- 0
      if(log) return(ans)
      else return(exp(ans))
    }
  )
  
  rPG <- nimbleFunction(
    run = function(n = integer(0)) {
      stop('rPG should never be run')
      returnType(double(1))
      return(0)
    })
  
  
  compute_predictor = nimbleFunction(
    run = function(X = double(2), beta= double(1),  indexes_covariates= double(1), gamma = double(1)){
      returnType(double(2))
      
      L <- length(indexes_covariates)
      index_present <- logical(L)
      
      for (i in 1:L) {
        if(any(indexes_covariates[i] ==which(gamma==1))){
          index_present[i] <- TRUE
        } else index_present[i] <- FALSE
        
      }
      
      X_gamma <- X[,index_present,drop = FALSE]
      beta_gamma <- beta[index_present]
      S <- dim(X_gamma)[1]
      linpred <- X_gamma %*% beta_gamma#[1:S,1]  
      
      return(linpred)
      
    }
  )
  
  
}

win.data <- list(S= M*T, D=M*(T-1), nsamples = length(y), nsites = M, numVars = numVars, b = b_psi, B = B_psi, zindex = zindex, ncov = ncov+1, indexes_covariates= indexes_covariates )
str(win.data)


nimOcc <- nimbleCode({
  
  # Priors
  
  # BVS on psi
  Omega_psi[1:nsites] ~ dPG()
  gamma_psi[1:ncov]  ~ dmnorm(b[1:ncov], B[1:ncov, 1:ncov])  
  beta_psi[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
  linpred_psi[1:nsites] <- compute_predictor(X_psi[1:nsites, 1:numVars], beta_psi[1:numVars], indexes_covariates[1:numVars],gamma_psi[1:ncov])[1:nsites,1]   
  
  
  # BVS on phi
  Omega_phi[1:D] ~ dPG()
  gamma_phi[1:ncov]  ~ dmnorm(b[1:ncov], B[1:ncov, 1:ncov])  
  beta_phi[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
  linpred_phi[1:D] <- compute_predictor(X_phi[1:D, 1:numVars], beta_phi[1:numVars], indexes_covariates[1:numVars],gamma_phi[1:ncov])[1:D,1]   
  
  
  # BVS on colonization
  Omega_eta[1:D] ~ dPG()
  gamma_eta[1:ncov]  ~ dmnorm(b[1:ncov], B[1:ncov, 1:ncov])  
  beta_eta[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
  linpred_eta[1:D] <- compute_predictor(X_eta[1:D, 1:numVars], beta_eta[1:numVars], indexes_covariates[1:numVars],gamma_eta[1:ncov])[1:D,1]   
  
  
  # BVS on p
  Omega_p[1:nsamples] ~ dPG()
  gamma_p[1:ncov]  ~ dmnorm(b[1:ncov], B[1:ncov, 1:ncov])  
  beta_p[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
  linpred_p[1:nsamples] <- compute_predictor(X_p[1:nsamples, 1:numVars], beta_p[1:numVars], indexes_covariates[1:numVars],gamma_p[1:ncov])[1:nsamples,1]   
  
  
  
  for (s in 1:nsites) {
    z[s] ~ dbern(psi[s])
    logit(psi[s]) <- linpred_psi[s]
  }
  
  for (d in 1:D) {
    logit(phi[d]) <- linpred_phi[d]
    logit(eta[d]) <- linpred_eta[d]
  }
  
  for (s in (nsites+1):S) {
    z[s] ~  dbern(z[s-nsites]*phi[s-nsites]  + (1-z[s-nsites])*eta[s-nsites]) #scalar representation on survival and colonization
  }
  
  for (n in 1:nsamples) {
    y[n] ~ dbern(z[zindex[n]]*p[n])
    logit(p[n]) <- linpred_p[n]
  }
  
  
})

# Inital value for z
zints = c(apply(data$Y, c(1,3), max))
# Initial value for PG parameter (Omega)
up = matrix(1,M*T*J,1)
Omega_p = pgdraw(rep(1, each = M*T*J),0.5*up)
up = matrix(1,win.data$D,1)
Omega_phi = pgdraw(rep(1, each = win.data$D),0.5*up)
up = matrix(1,win.data$D,1)
Omega_eta = pgdraw(rep(1, each = win.data$D),0.5*up)
u = matrix(1,M,1)
Omega_psi = pgdraw(rep(1, each = M),0.5*u)


ints = list(z = zints, beta_phi = rep(0, numVars), gamma_phi = rep(1,ncov+1), Omega_phi = Omega_phi ,beta_eta = rep(0, numVars), gamma_eta = rep(1,ncov+1), Omega_eta= Omega_eta ,beta_p = rep(0, numVars), gamma_p = rep(1,ncov+1), Omega_p= Omega_p, beta_psi = rep(0, numVars), gamma_psi = rep(1,ncov+1), Omega_psi= Omega_psi )
Occmodel = nimbleModel(nimOcc, constants = win.data,data = list(y=y, X_p=X_p, X_phi= X_phi, X_eta= X_eta, X_psi = X_psi), inits = ints)
Occmodel$calculate()


cOccmodel <- compileNimble(Occmodel )
Occconf = configureMCMC(Occmodel)

# Removing NIMBLE's default samplers so we can add our samplers to preform BVS
Occconf$removeSampler(c("Omega_p", "gamma_p", "beta_p", "Omega_phi", "gamma_phi", "beta_phi","Omega_eta", "gamma_eta", "beta_eta", "Omega_psi", "gamma_psi", "beta_psi" ))
Occconf$addSampler(target = "Omega_p", type = "PG_sampler_p", control = list(indexes_covariates = indexes_covariates, zindex= zindex ) )
Occconf$addSampler(target = "gamma_p", type = "gamma_sampler_p",control = list(ncov= ncov, indexes_covariates = indexes_covariates, zindex= zindex ))
Occconf$addSampler(target = "beta_p", type = "beta_sampler_p", control = list(indexes_covariates = indexes_covariates, zindex= zindex) )
Occconf$addSampler(target = "Omega_phi", type = "PG_sampler_phi", control = list(indexes_covariates = indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
Occconf$addSampler(target = "gamma_phi", type = "gamma_sampler_phi",control = list(ncov= ncov, indexes_covariates = indexes_covariates,   nsites = win.data$nsites, nyears = T))
Occconf$addSampler(target = "beta_phi", type = "beta_sampler_phi", control = list(indexes_covariates = indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
Occconf$addSampler(target = "Omega_eta", type = "PG_sampler_eta", control = list(indexes_covariates = indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
Occconf$addSampler(target = "gamma_eta", type = "gamma_sampler_eta",control = list(ncov= ncov, indexes_covariates = indexes_covariates,   nsites = win.data$nsites, nyears = T))
Occconf$addSampler(target = "beta_eta", type = "beta_sampler_eta", control = list(indexes_covariates = indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
Occconf$addSampler(target = "Omega_psi", type = "PG_sampler_psi", control = list(indexes_covariates = indexes_covariates ) )
Occconf$addSampler(target = "gamma_psi", type = "gamma_sampler_psi",control = list(ncov= ncov, indexes_covariates = indexes_covariates ))
Occconf$addSampler(target = "beta_psi", type = "beta_sampler_psi", control = list(indexes_covariates = indexes_covariates) )

Occconf$resetMonitors()
Occconf$addMonitors(c("beta_p",  "gamma_p", "beta_psi", "gamma_psi", "beta_phi", "gamma_phi", "beta_eta", "gamma_eta" ))


Occmcmc = buildMCMC(Occconf)
Occmod = compileNimble(Occmcmc, project = Occmodel, resetFunctions = TRUE)
Occmcmc.out <- runMCMC(Occmod, niter = 25000, nburnin = 10000, thin= 5, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
Occmcmc.out$summary$all.chains


MCMCvis::MCMCsummary(Occmcmc.out$samples) # Quick check of convergence and ESS

#=====================================================================================================================================================================

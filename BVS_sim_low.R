
# Rscript simulation code when M=500, J=5 and high levels on initial occupancy, persistence and colonization

rm(list= ls())
library(nimble)
library(coda)
library(pgdraw)
library(data.table)
source("BVS_DOM_samplers.R") # load MCMC samplers

M = 500 # sites 
J = 5 # sampling occasions
T = 10 # years
ncov = 10
nsim = 10
seed_numbers = 1:nsim

true_betas_psi <- c(-0.25, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
true_betas_phi <- c(-1, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
true_betas_eta <- c(-1.5, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
true_betas_p <- c(-0.25, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept

# Arrays to store MCMC summaries
beta_p = beta_psi = beta_phi= beta_eta= array(0, dim = c(15, 5, nsim))
gamma_p = gamma_psi= gamma_phi= gamma_eta= array(0, dim = c(ncov+1, 5, nsim))

simData <- function(M, J, T,true_betas_psi ,true_betas_phi, true_betas_eta, true_betas_p ){
  
  
  S = M*(T-1)
  N =  M*T*J
  
  ncov = 10 # number of explanatory variables in the model excluding the intercept
  # first 5 continuous and last 5 categorical 
  
  
  # Design matrix for psi
  X = matrix(0, nrow = M, ncol = 5)
  X[,1] = runif(M, -2,2)*(sqrt(3)/2)
  X[,2] = runif(M, -3,3)*(2/(3*sqrt(3)))
  X[,3:5] = rnorm(M*3, 0,1)
  
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  X6 = rbinom(M, 1, 0.8)
  X7 = data.frame(x7 = as.factor(sample(x=c(1,2, 3,4), size=M, replace=TRUE, prob=c(0.1,0.1,0.4,0.4) ) ))
  X8 = data.frame(x8 = as.factor(sample(x=c(1,2, 3), size=M, replace=TRUE, prob=c(0.3,0.3,0.4)  ) ) )
  X9 = rbinom(M,1, 0.4)
  X10 =  data.frame(x10 = as.factor(sample(x=c(1,2, 3), size=M, replace=TRUE, prob=c(0.1,0.3,0.6)  ) ) )
  
  originalX =cbind(originalX, X6, X7, X8, X9, X10)
  #head(originalX)
  
  X_psi <- model.matrix(~ ., originalX) # creates a design (or model) matrix
  #head(X_psi)
  
  psi <-  ilogit(X_psi %*% true_betas_psi) # initial occupancy probability 
  #mean(psi)
  
  # Design matrix for survival (phi) 
  X = matrix(0, nrow = S, ncol = 5)
  X[,1] = runif(S, -2,2)*(sqrt(3)/2)
  X[,2] = runif(S, -3,3)*(2/(3*sqrt(3)))
  X[,3:5] = rnorm(S*3, 0,1)
  
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  X6 = rbinom(M, 1, 0.8)
  X7 = data.frame(x7 = as.factor(sample(x=c(1,2, 3,4), size=S, replace=TRUE, prob=c(0.1,0.1,0.4,0.4) ) ))
  X8 = data.frame(x8 = as.factor(sample(x=c(1,2, 3), size=S, replace=TRUE, prob=c(0.3,0.3,0.4)  ) ) )
  X9 = rbinom(S,1, 0.4)
  X10 =  data.frame(x10 = as.factor(sample(x=c(1,2, 3), size=S, replace=TRUE, prob=c(0.1,0.3,0.6)  ) ) )
  
  originalX =cbind(originalX, X6, X7, X8, X9, X10)
  #head(originalX)
  
  X_phi <- model.matrix(~ ., originalX) # creates a design (or model) matrix
  #head(X_phi)
  
  phi <-  matrix(ilogit(X_phi %*% true_betas_phi),M,T-1) # occupancy probability 
  #mean(apply(phi, 2, mean))
  
  
  # Design matrix for colonization (eta) 
  X_eta = X_phi
  eta <-  matrix(ilogit(X_eta %*% true_betas_eta),M,T-1) # occupancy probability 
  #mean(apply(eta, 2, mean))
  
  
  # Design matrix for p
  X = matrix(0, nrow = N, ncol = 5)
  X[,1] = runif(N, -2,2)*(sqrt(3)/2)
  X[,2] = runif(N, -3,3)*(2/(3*sqrt(3)))
  X[,3:5] = rnorm(N*3, 0,1)
  
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  X6 = rbinom(M, 1, 0.8)
  X7 = data.frame(x7 = as.factor(sample(x=c(1,2, 3,4), size=N, replace=TRUE, prob=c(0.1,0.1,0.4,0.4) ) ))
  X8 = data.frame(x8 = as.factor(sample(x=c(1,2, 3), size=N, replace=TRUE, prob=c(0.3,0.3,0.4)  ) ) )
  X9 = rbinom(N,1, 0.4)
  X10 =  data.frame(x10 = as.factor(sample(x=c(1,2, 3), size=N, replace=TRUE, prob=c(0.1,0.3,0.6)  ) ) )
  
  originalX =cbind(originalX, X6, X7, X8, X9, X10)
  #head(originalX)
  
  X_p <- model.matrix(~ ., originalX) # creates a design (or model) matrix
  #head(X_p)
  
  
  p <-  array(ilogit(X_p %*% true_betas_p), dim= c(M,J,T)) # occupancy probability
  #mean(apply(p, 2, mean))
  
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
  
  
  column_covariate = 1:ncov
  
  # true or false if the covaraites is categorical or numerical 
  classCovariates <- c(rep(TRUE, 5), rep(FALSE, 5)) # first 4 continouous and last one categorical with 3 levels
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
  
  #indexes_covariates 
  
  
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
  
  #C
  Sigma <- cov(X_p[,-1, drop = F])
  sumCSigma <- sum(C * Sigma)
  
  mu_psi = 0.5
  phi_mu =4
  phi_beta = 1/4
  
  # Prior parameters for beta, beta ~ MVN(b, B)
  b_psi <- c(mu_psi, rep(0, ncol(X_p) - 1))
  sigma_beta_psi <- C * phi_beta
  #b_psi #b
  
  B_psi <- matrix(0, nrow = ncol(X_p), ncol = ncol(X_p))
  B_psi[1, 1] <- phi_mu
  B_psi[2:(ncol(X_p)), 2:(ncol(X_p))] <- sigma_beta_psi
  #B_psi # B
  
  numVars = dim(B_psi)[1] #number of covariates (including the intercept)
  
  
  return(list(y=y, Y=Y, X_psi=X_psi, X_phi=X_phi, X_eta=X_eta, X_p=X_p, indexes_covariates= indexes_covariates, B_psi=B_psi, b_psi=b_psi, numVars= numVars))
  
}

# shell distribution for PG
dPG <- nimbleFunction(
  run = function(x = double(1), log=integer(0, default=1)){
    returnType(double(0))
    ans <- 0
    if(log) return(ans)
    else return(exp(ans))
  }
)

# some non-significant dist

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


for (i in 1:nsim) {
  
  
  #  Simulated  data 
  set.seed(seed_numbers[i])
  
  data = simData(M,J, T, true_betas_psi, true_betas_phi, true_betas_eta,true_betas_p)
  
  # indicator variable for z such that observation process can be scalar
  nsamples = length(data$y)
  zindex = numeric(nsamples)
  zindex[1:(M*J)] = rep(1:M,J)
  
  for (t in 2:T) {
    zindex[((((t-1) *M*J)+1)): (M*J*t) ] = rep( ((t-1)*M+1):(t*M),J)
  }
  
  win.data <- list(S= M*T, D= M*(T-1), nsamples = nsamples, nsites = M, numVars = data$numVars, b = data$b_psi, B = data$B_psi, zindex = zindex, ncov = ncov+1, indexes_covariates= data$indexes_covariates )
  str(win.data)
  
  # Inital value for z
  zints = c(apply(data$Y, c(1,3), max))
  # Initial value for PG parameter (Omega)
  up = matrix(1,nsamples,1)
  Omega_p = pgdraw(rep(1, each = nsamples),0.5*up)
  up = matrix(1,win.data$D,1)
  Omega_phi = pgdraw(rep(1, each = win.data$D),0.5*up)
  up = matrix(1,win.data$D,1)
  Omega_eta = pgdraw(rep(1, each = win.data$D),0.5*up)
  u = matrix(1,M,1)
  Omega_psi = pgdraw(rep(1, each = M),0.5*u)
  
  
  
  
  
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
  
  
  
  ints = list(z = zints, beta_phi = rep(0, data$numVars), gamma_phi = rep(1,ncov+1), Omega_phi = Omega_phi ,beta_eta = rep(0, data$numVars), gamma_eta = rep(1,ncov+1), Omega_eta= Omega_eta ,beta_p = rep(0, data$numVars), gamma_p = rep(1,ncov+1), Omega_p= Omega_p, beta_psi = rep(0, data$numVars), gamma_psi = rep(1,ncov+1), Omega_psi= Omega_psi )
  Occmodel = nimbleModel(nimOcc, constants = win.data,data = list(y=data$y, X_p=data$X_p, X_phi= data$X_phi, X_eta= data$X_eta, X_psi = data$X_psi), inits = ints)
  Occmodel$calculate()
  
  
  cOccmodel <- compileNimble(Occmodel )
  Occconf = configureMCMC(Occmodel)
  
  # Removing NIMBLE's default samplers so we can add the samplers to preform BVS
  Occconf$removeSampler(c("Omega_p", "gamma_p", "beta_p", "Omega_phi", "gamma_phi", "beta_phi","Omega_eta", "gamma_eta", "beta_eta", "Omega_psi", "gamma_psi", "beta_psi" ))
  Occconf$addSampler(target = "Omega_p", type = "PG_sampler_p", control = list(indexes_covariates = data$indexes_covariates, zindex= zindex ) )
  Occconf$addSampler(target = "gamma_p", type = "gamma_sampler_p",control = list(ncov= ncov, indexes_covariates = data$indexes_covariates, zindex= zindex ))
  Occconf$addSampler(target = "beta_p", type = "beta_sampler_p", control = list(indexes_covariates = data$indexes_covariates, zindex= zindex) )
  Occconf$addSampler(target = "Omega_phi", type = "PG_sampler_phi", control = list(indexes_covariates = data$indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
  Occconf$addSampler(target = "gamma_phi", type = "gamma_sampler_phi",control = list(ncov= ncov, indexes_covariates = data$indexes_covariates,   nsites = win.data$nsites, nyears = T))
  Occconf$addSampler(target = "beta_phi", type = "beta_sampler_phi", control = list(indexes_covariates = data$indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
  Occconf$addSampler(target = "Omega_eta", type = "PG_sampler_eta", control = list(indexes_covariates = data$indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
  Occconf$addSampler(target = "gamma_eta", type = "gamma_sampler_eta",control = list(ncov= ncov, indexes_covariates = data$indexes_covariates,   nsites = win.data$nsites, nyears = T))
  Occconf$addSampler(target = "beta_eta", type = "beta_sampler_eta", control = list(indexes_covariates = data$indexes_covariates,   nsites =  win.data$nsites, nyears = T) )
  Occconf$addSampler(target = "Omega_psi", type = "PG_sampler_psi", control = list(indexes_covariates = data$indexes_covariates ) )
  Occconf$addSampler(target = "gamma_psi", type = "gamma_sampler_psi",control = list(ncov= ncov, indexes_covariates = data$indexes_covariates ))
  Occconf$addSampler(target = "beta_psi", type = "beta_sampler_psi", control = list(indexes_covariates = data$indexes_covariates) )
  Occconf$resetMonitors()
  Occconf$addMonitors(c("beta_p",  "gamma_p", "beta_psi", "gamma_psi", "beta_phi", "gamma_phi", "beta_eta", "gamma_eta" ))
  
  
  Occmcmc = buildMCMC(Occconf)
  Occmod = compileNimble(Occmcmc, project = Occmodel, resetFunctions = TRUE)
  #Occmcmc.out <- runMCMC(Occmod, niter = 25000, nburnin = 10000, thin= 5, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
  Occmcmc.out <- runMCMC(Occmod, niter = 45000, nburnin = 20000, thin= 5, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
  beta_eta[,,i]= Occmcmc.out$summary$all.chains[1:15,]
  beta_p[,,i]= Occmcmc.out$summary$all.chains[16:30,]
  beta_phi[,,i]= Occmcmc.out$summary$all.chains[31:45,]
  beta_psi[,,i]= Occmcmc.out$summary$all.chains[46:60,]
  gamma_eta[,,i]= Occmcmc.out$summary$all.chains[61:71,]
  gamma_p[,,i]= Occmcmc.out$summary$all.chains[72:82,]
  gamma_phi[,,i]= Occmcmc.out$summary$all.chains[83:93,]
  gamma_psi[,,i]= Occmcmc.out$summary$all.chains[94:104,]
  
  
  save(beta_eta, file = "BVS_DOM_1_beta_eta.Rdata")
  save(beta_p, file = "BVS_DOM_1_beta_p.Rdata")
  save(beta_phi, file = "BVS_DOM_1_beta_phi.Rdata")
  save(beta_psi, file = "BVS_DOM_1_beta_psi.Rdata")
  save(gamma_eta, file = "BVS_DOM_1_gamma_eta.Rdata")
  save(gamma_p, file = "BVS_DOM_1_gamma_p.Rdata")
  save(gamma_phi, file = "BVS_DOM_1_gamma_phi.Rdata")
  save(gamma_psi, file = "BVS_DOM_1_gamma_psi.Rdata")
  
  
}

#=====================================================================================================================================================================

# Results 

true_betas_psi <- c(1.8, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
true_betas_eta <- c(1, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
true_betas_p <- c(1.8, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept


library(ggmcmc)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/M500J5_all_high/BVS_DOM_1_gamma_p.Rdata")

head(gamma_p)
pip_p = t(gamma_p[-1,1,])
pip_p # looks good
true_betas_p <- c(1.8, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept

colnames(pip_p)= 1:10

ggs_caterpillar(ggs(as.mcmc(pip_p), keep_original_order = T), horizontal = FALSE, sort = FALSE, thin_ci = c(0.5,0.5))  + xlab("PIP") + ylab("Variable") + ggtitle(" ") + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_vline(xintercept=0.5,linetype=2)


load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/M500J5_all_high/BVS_DOM_1_gamma_phi.Rdata")

head(gamma_phi)
pip_phi = t(gamma_phi[-1,1,]) # take the mean
pip_phi # looks good
true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 

colnames(pip_phi)= 1:10

ggs_caterpillar(ggs(as.mcmc(pip_phi), keep_original_order = T), horizontal = FALSE, sort = FALSE, thin_ci = c(0.5,0.5)) + xlab("PIP") + ylab("Variable") + ggtitle(" ") + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_vline(xintercept=0.5,linetype=2)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/M500J5_all_high/BVS_DOM_1_gamma_psi.Rdata")

head(gamma_psi)
pip_psi = t(gamma_psi[-1,1,]) # take the mean
pip_psi # looks good
true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 

colnames(pip_psi)= 1:10

ggs_caterpillar(ggs(as.mcmc(pip_psi), keep_original_order = T), horizontal = FALSE, sort = FALSE, thin_ci = c(0.5,0.5)) + xlab("PIP") + ylab("Variable") + ggtitle(" ") + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_vline(xintercept=0.5,linetype=2)

load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/M500J5_all_high/BVS_DOM_1_gamma_eta.Rdata")

head(gamma_eta)
pip_eta = t(gamma_eta[-1,1,]) # take the mean
pip_eta # looks good
true_betas_eta <- c(1, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept

colnames(pip_eta)= 1:10

ggs_caterpillar(ggs(as.mcmc(pip_eta), keep_original_order = T), horizontal = FALSE, sort = FALSE, thin_ci = c(0.5,0.5)) + xlab("PIP") + ylab("Variable") + ggtitle(" ") + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_vline(xintercept=0.5,linetype=2)

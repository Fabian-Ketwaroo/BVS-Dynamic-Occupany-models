#=====================================================================================================================================
# SOM model code with BVS for Northern owls

rm(list= ls())
library(nimble)
library(coda)
library(pgdraw)
library(data.table)
source("BVS_Northern_owls_samplers.R") # load MCMC samplers

load("D:/BVS Dyanmic Occupancy/BVS DOM/Northern spotted owls.Rdata")
head(northern_owls)
M = dim(northern_owls)[1]
J= 6

Y = as.matrix(northern_owls[, 1:6], nrow = M, ncol=6)
head(Y)
y_na = c(Y) # data in scalar form 
y =y_na[ which(is.na(y_na)==0)] 
(nsamples = length(y))

Jindex =c()
for (i in 1:M) {
  Jindex[i] <- sum(!is.na(Y[i,]))  
}

zindex_array =  array(NA, dim = c(M,J))

for (i in 1:M) {
  zindex_array[i, 1:Jindex[i]] = rep(i,Jindex[i])
}


zindex_na = c(zindex_array)
zindex =zindex_na[ which(is.na(zindex_na)==0)] 
summary(zindex)
length(zindex) == nsamples

#=======================================================================================================================================
# Covariates for p
#=======================================================================================================================================

head(northern_owls)
# Time of sampling- day(0) or night(1- 
tim = as.matrix(northern_owls[,9:14 ], nrow = M, ncol=6)
head(tim)
tim_na = c(tim) # data in scalar form 
day_night =tim_na[ which(is.na(tim_na)==0)] 
length(day_night) == length(y)
day_night = as.factor(day_night+1)


# BAOW
Ba = cbind(northern_owls$BAOW, replicate(5,northern_owls$BAOW))
head(Ba)
dim(Ba)
Ba_na = c(Ba) # data in scalar form 
BAOW = Ba[which(is.na(y_na)==0)]
length(BAOW) == length(y)
BAOW = as.factor(BAOW+1)

# Edge
northern_owls$EDGE_scale = scale(northern_owls$EDGE)
Ed = cbind(northern_owls$EDGE_scale[,1], replicate(5,northern_owls$EDGE_scale[,1]))
head(Ed)
Ed_na = c(Ed) # data in scalar form 
EDGE = Ed[which(is.na(y_na)==0)]
length(EDGE) == length(y)



originalX = as.data.frame(cbind(BAOW, EDGE, day_night))
originalX$BAOW = as.factor(originalX$BAOW)
originalX$day_night = as.factor(originalX$day_night)
head(originalX)

X = originalX
head(X)
#X = as.data.frame(X)
X <- model.matrix(~ ., X) # creates a design (or model) matrix
X_p = X
head(X_p)

ncov_p = 3
column_covariate = 1:ncov_p

# true or false if the covaraites is categorical or numerical 
classCovariates <- c(FALSE, TRUE, FALSE) 
classCovariates


indexes_covariates <- c()
indexes_covariates[1] <- 1 # fix the intercept

if(any(column_covariate != 0)){
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  k <- 2
  for (i in 1:ncov_p) {
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
indexes_covariates_p = indexes_covariates  

# compute prior correlated (C) matrix 
C <- matrix(0, nrow = ncol(X_p) - 1, ncol = ncol(X_p) - 1)
l <- 0 # l loop through the covariates
for (i in 1:ncov_p) {
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

d_bar_p <- ncov_p # prior beliefs about the number of important covariates
d_bar_p <- ifelse(d_bar_p <= ncov_p, d_bar_p, 2)

mu_psi = 0.5
phi_mu =4
phi_beta = 1/4

# Prior parameters for beta, beta ~ MVN(b, B)
b_p <- c(mu_psi, rep(0, ncol(X_p) - 1))
sigma_beta_psi <- C * phi_beta
b_p #b

B_p <- matrix(0, nrow = ncol(X_p), ncol = ncol(X_p))
B_p[1, 1] <- phi_mu
B_p[2:(ncol(X_p)), 2:(ncol(X_p))] <- sigma_beta_psi
B_p # B

numVars_p = dim(B_p)[1] #number of covariates (including the intercept)

#==========================================================================================================================================
# Covariates for psi
#=======================================================================================================================================

head(northern_owls)
X = northern_owls[,7:8] 
head(X)
X[, "BAOW"] = as.factor(X[, "BAOW"]+1)
originalX = X
X = as.data.frame(X)
X <- model.matrix(~ ., X) # creates a design (or model) matrix
head(X)
X[,3] = scale(X[,3])
X_psi = X
head(X_psi)

ncov_psi = 2
column_covariate = 1:ncov_psi

# true or false if the covaraites is categorical or numerical 
classCovariates <- c(FALSE, TRUE) 
classCovariates


indexes_covariates <- c()
indexes_covariates[1] <- 1 # fix the intercept

if(any(column_covariate != 0)){
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  k <- 2
  for (i in 1:ncov_psi) {
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
indexes_covariates_psi = indexes_covariates  

# compute prior correlated (C) matrix 
C <- matrix(0, nrow = ncol(X_psi) - 1, ncol = ncol(X_psi) - 1)
l <- 0 # l loop through the covariates
for (i in 1:ncov_psi) {
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
Sigma <- cov(X_psi[,-1, drop = F])
sumCSigma <- sum(C * Sigma)

d_bar <- ncov_psi # prior beliefs about the number of important covariates
d_bar <- ifelse(d_bar <= ncov_psi, d_bar, 2)

mu_psi = 0.5
phi_mu =4
phi_beta = 1/4

# Prior parameters for beta, beta ~ MVN(b, B)
b_psi <- c(mu_psi, rep(0, ncol(X_psi) - 1))
sigma_beta_psi <- C * phi_beta
b_psi #b

B_psi <- matrix(0, nrow = ncol(X_psi), ncol = ncol(X_psi))
B_psi[1, 1] <- phi_mu
B_psi[2:(ncol(X_psi)), 2:(ncol(X_psi))] <- sigma_beta_psi
B_psi # B

numVars_psi = dim(B_psi)[1] #number of covariates (including the intercept)

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


win.data <- list( nsamples = length(y), nsites = M, numVars_p = numVars_p, numVars_psi = numVars_psi,b_p = b_p, B_p = B_p,b_psi = b_psi, B_psi = B_psi, zindex = zindex,ncov_psi = ncov_psi+1, ncov_p = ncov_p+1, indexes_covariates_p= indexes_covariates_p, indexes_covariates_psi= indexes_covariates_psi )
str(win.data)


nimOcc <- nimbleCode({
  
  # Priors
  
  # BVS on psi
  Omega_psi[1:nsites] ~ dPG()
  gamma_psi[1:ncov_psi]  ~ dmnorm(b_psi[1:ncov_psi], B_psi[1:ncov_psi, 1:ncov_psi])  
  beta_psi[1:numVars_psi] ~ dmnorm(b_psi[1:numVars_psi], B_psi[1:numVars_psi, 1:numVars_psi]) 
  linpred_psi[1:nsites] <- compute_predictor(X_psi[1:nsites, 1:numVars_psi], beta_psi[1:numVars_psi], indexes_covariates_psi[1:numVars_psi],gamma_psi[1:ncov_psi])[1:nsites,1]   
  
  
  # BVS on p
  Omega_p[1:nsamples] ~ dPG()
  gamma_p[1:ncov_p]  ~ dmnorm(b_p[1:ncov_p], B_p[1:ncov_p, 1:ncov_p])  
  beta_p[1:numVars_p] ~ dmnorm(b_p[1:numVars_p], B_p[1:numVars_p, 1:numVars_p]) 
  linpred_p[1:nsamples] <- compute_predictor(X_p[1:nsamples, 1:numVars_p], beta_p[1:numVars_p], indexes_covariates_p[1:numVars_p],gamma_p[1:ncov_p])[1:nsamples,1]   
  
  
  
  for (s in 1:nsites) {
    z[s] ~ dbern(psi[s])
    logit(psi[s]) <- linpred_psi[s]
  }
  
  
  for (n in 1:nsamples) {
    y[n] ~ dbern(z[zindex[n]]*p[n])
    logit(p[n]) <- linpred_p[n]
  }
  
  
})

# Inital value for z
zints =  rep(1, M)#c(apply(Y, c(1,3), max)) 
# Initial value for PG parameter (Omega)
up = matrix(1,nsamples,1)
Omega_p = pgdraw(rep(1, each = nsamples),0.5*up)
u = matrix(1,M,1)
Omega_psi = pgdraw(rep(1, each = M),0.5*u)

#ints = list(z = zints, beta_phi = rep(0, numVars), gamma_phi = rep(1,ncov+1), Omega_phi = Omega_phi ,beta_eta = rep(0, numVars), gamma_eta = rep(1,ncov+1), Omega_eta= Omega_eta ,beta_p = rep(0, numVars), gamma_p = rep(1,ncov+1), Omega_p= Omega_p, beta_psi = rep(0, numVars), gamma_psi = rep(1,ncov+1), Omega_psi= Omega_psi )
ints = list(z = zints, beta_p = rep(0, numVars_p), gamma_p = rep(1,ncov_p+1), Omega_p= Omega_p, beta_psi = rep(0, numVars_psi), gamma_psi = rep(1,ncov_psi+1), Omega_psi= Omega_psi )
Occmodel = nimbleModel(nimOcc, constants = win.data,data = list(y=y, X_p=X_p, X_psi = X_psi), inits = ints)
Occmodel$calculate()

cOccmodel <- compileNimble(Occmodel )
Occconf = configureMCMC(Occmodel)

# Removing NIMBLE's default samplers so we can add the samplers to preform BVS
Occconf$removeSampler(c("Omega_p", "gamma_p", "beta_p", "Omega_psi", "gamma_psi", "beta_psi" ))
Occconf$addSampler(target = "Omega_p", type = "PG_sampler_p", control = list(indexes_covariates = indexes_covariates_p, zindex= zindex ) )
Occconf$addSampler(target = "gamma_p", type = "gamma_sampler_p",control = list(ncov= ncov_p, indexes_covariates = indexes_covariates_p, zindex= zindex ))
Occconf$addSampler(target = "beta_p", type = "beta_sampler_p", control = list(indexes_covariates = indexes_covariates_p, zindex= zindex) )
Occconf$addSampler(target = "Omega_psi", type = "PG_sampler_psi", control = list(indexes_covariates = indexes_covariates_psi ) )
Occconf$addSampler(target = "gamma_psi", type = "gamma_sampler_psi",control = list(ncov= ncov_psi, indexes_covariates = indexes_covariates_psi ))
Occconf$addSampler(target = "beta_psi", type = "beta_sampler_psi", control = list(indexes_covariates = indexes_covariates_psi) )
Occconf$resetMonitors()
Occconf$addMonitors(c("beta_p",  "gamma_p", "beta_psi", "gamma_psi"))


Occmcmc = buildMCMC(Occconf)
Occmod = compileNimble(Occmcmc, project = Occmodel, resetFunctions = TRUE)
Occmcmc.out <- runMCMC(Occmod, niter = 25000, nburnin = 10000, thin= 5, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
Occmcmc.out$summary$all.chains

MCMCvis::MCMCsummary(Occmcmc.out$samples) # Quick check of convergence and ESS




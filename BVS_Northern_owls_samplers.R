#PG BVS samplers for occupany models

diagMatrixProd <- nimbleFunction(
  run = function( X = double(2), D= double(1) ){
    returnType(double(2))
    
    result <- matrix(0, dim(X)[1], length(D))
    
    for (i in 1:dim(result)[1]) {
      for (j in 1:dim(result)[2]) {
        result[i, j] = X[i,j] * D[j]
      }
      
    }
    
    return(result)
  }
)


# Gamma log likelihood
dLgamma <- nimbleFunction(
  run= function( x = double(1) , X= double(2),  indexes_covariates = double(1), b = double(1),  B= double(2),  Omega= double(1),  k= double(1), log = integer(0, default= 1) ){
    returnType(double(0))
    
    gamma <- x
    m <- length(indexes_covariates)
    index_present <-  numeric(m)
    l = 0
    
    #gamma = rep(1, numVars) # data to test likelihood works correctly
    #  gamma = c(1,1,1,0,0,0)
    
    for (i in 1:m ) {
      if(gamma[indexes_covariates[i]] == 1 ){
        index_present[l+1] = i
        l = l+1
      }
    }
    
    
    X_gamma <- matrix(0, dim(X)[1], l)
    b_gamma <- numeric(l)
    B_gamma <- matrix(0, l, l)
    
    
    for (i in 1:l) {
      X_gamma[,i] = X[, index_present[i]]
      b_gamma[i] = b[index_present[i]]
      for (j in 1:l) {
        B_gamma[i,j] <- B[index_present[i], index_present[j]]
        
      }
    }
    
    tX = t(X_gamma)
    tXOmega = diagMatrixProd(tX, Omega)
    cholXgOmX = chol(tXOmega  %*% X_gamma + inverse(B_gamma) )
    
    firstTerm = (.5) * logdet(inverse(B_gamma)) - logdet(cholXgOmX)
    tXKbplusBb = t(X_gamma) %*% k + inverse(B_gamma) %*% b_gamma  
    v = solve(t(cholXgOmX), tXKbplusBb) # problem with trimatl and trimatl and so I've taken it out 
    vtv = t(v) %*% v
    
    secondTerm = - .5 * ( (t(b_gamma) %*% inverse(B_gamma) %*% b_gamma) - vtv)
    
    loglikelihood <- firstTerm + secondTerm[1,1]
    
    if(log) return(loglikelihood)
    else return(exp(loglikelihood))
    
  }
)



# Sampler for Polya Gamma variables (Omega)

aterm <- nimbleFunction(
  run = function(n = integer(), x= double(0), t = double(0)){
    returnType(double(0))
    f=0
    if(x <= t) {
      f = log(pi) + log(n + 0.5) + 1.5*(-log(pi/2) - log(x)) - 2*(n + 0.5)*(n + 0.5)/x
    }
    else {
      f = log(pi) + log(n + 0.5) - x * pi^2/2 * (n + 0.5)*(n + 0.5)
    }    
    return(exp(f))
    
  }
)

exprnd <- nimbleFunction(
  run = function(mu = double(0)){
    returnType(double(0))
    return(-mu*log(1.0 - runif(1,0.0,1.0) ))
  }
)


truncgamma <- nimbleFunction(
  run=function( ){
    returnType(double(0))
    
    c <- pi/2
    done <- FALSE
    while(!done){
      X = exprnd(1.0) * 2.0 + c
      gX = sqrt(pi/2) / sqrt(X)
      
      if(runif(1,0,1) <= gX ) done = TRUE
    }
    
    return(X)
  }
)



randinvg <- nimbleFunction(
  run=function(mu = double(0)){
    returnType(double(0))
    u = rnorm(1,0,1)
    V = u*u
    out = mu + 0.5*mu *( mu*V - sqrt(4.0*mu*V + mu*mu * V*V) )
    
    if(runif(1,0,1) > mu /(mu+out) )  out = mu*mu / out
    
    return(out)
  }
)

tinvgauss <- nimbleFunction(
  run= function(z= double(0), t = double(0)) {
    returnType(double(0))
    
    mu <- 1.0/z
    
    done <- FALSE
    
    # Pick sampler
    if(mu > t) {
      # Sampler based on truncated gamma 
      # Algorithm 3 in the Windle (2013) PhD thesis, page 128
      while(!done) {
        u = runif(1,0.0, 1.0)
        X = 1.0 / truncgamma()
        
        if ( log(u) < (-z*z*0.5*X) ) {
          #nimbreak()
          done = TRUE
          #break
        }
      }
    }  
    else {
      # Rejection sampler
      X = t + 1.0
      while(X >= t) {
        X = randinvg(mu)
      }
    } 
    
    return(X)
    
  }
)


samplepg <- nimbleFunction(
  run= function( z= double(0) ){
    returnType(double(0))
    
    # PG(b, z) = 0.25 * J*(b, z/2)
    z = abs(z)*0.5
    # Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
    MATH_2_PI <- 0.636619772367581343075535053490057448137838582961825794990
    t = MATH_2_PI
    
    done <- breakWhile <- FALSE
    
    # Compute p, q and the ratio q / (q + p)
    # (derived from scratch; derivation is not in the original paper)
    K = z*z/2.0 + pi^2/8.0
    logA = log(4.0) - log(pi) - z
    logK = log(K)
    Kt = K * t
    w = sqrt(pi/2)
    
    logf1 = logA + pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt
    logf2 = logA + 2*z + pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt
    p_over_q = exp(logf1) + exp(logf2)
    ratio = 1.0 / (1.0 + p_over_q)
    
    #breakWhile <- FALSE
    
    # Main sampling loop; page 130 of the Windle PhD thesis
    while(!breakWhile) {
      # Step 1: Sample X ? g(x|z)
      u = runif(1,0.0,1.0)
      if(u < ratio) {
        # truncated exponential
        X = t + exprnd(1.0)/K
      }
      else{
        # truncated Inverse Gaussian
        X = tinvgauss(z, t)
      }
      
      # Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
      i = 1
      Sn = aterm(0, X, t)
      U = runif(1,0.0,1.0) * Sn
      asgn = -1
      even = FALSE
      
      while(!done) {
        
        Sn = Sn + asgn * aterm(i, X, t)
        
        # Accept if n is odd
        if(!even & (U <= Sn)) {
          X = X * 0.25
          return(X)
        }
        
        # Return to step 1 if n is even
        if(even & (U > Sn) ) {
          # nimbreak()
          breakWhile = done = TRUE 
          # break # this is where the problem is 
        }
        
        even = !even
        asgn = -asgn
        i = i+1
        
      }
      
    }
    return(X) 
  }
)



pgr <- nimbleFunction(
  run = function( n = integer(), z = double(0) ){
    returnType(double(0))
    
    x = 0
    for (i in 1:n) {
      x =  x + samplepg(z)
    }
    
    return(x)
    
  }
)




sample_Omega = nimbleFunction(
  run = function( X = double(2), beta = double(1), n = double(1) ){
    returnType(double(1))
    
    nsize = length(n) # n is the number of samples for each observation
    Omega_vec = numeric(nsize)
    m <- dim(X)[2]
    
    for (i in 1:nsize) {
      b = X[i, 1:m ] %*%  beta[1:m] 
      Omega_vec[i] = pgr(n[i], b[1,1]) 
    }
    
    return(Omega_vec)
    
  }
)



dMatrixProd <- nimbleFunction(
  run = function( X = double(2), D= double(1) ){
    returnType(double(2))
    
    result <- matrix(0, dim(X)[1], length(D))
    
    for (i in 1:dim(result)[1]) {
      for (j in 1:dim(result)[2]) {
        result[i, j] = X[i,j] * D[j]
      }
      
    }
    
    return(result)
  }
)

mvrnormQuick <- nimbleFunction(
  run = function(mu = double(1), cholsigma = double(2)){
    returnType(double(1))
    ncols = dim(cholsigma)[2] 
    Y =  rnorm(ncols, 0,1)  
    return (c(mu + cholsigma %*% Y)) # shouldn' this be a matrix?, no it returns a vector in c++
    #return ( Y) # shouldn' this be a matrix? 
  }
)


sample_beta <- nimbleFunction(
  run = function(X= double(2), B= double(2), b = double(1), Omega = double(1) ,k= double(1) ){
    returnType(double(1))
    
    tX <- t(X)
    tXOmega = dMatrixProd(tX, Omega)
    
    L = t(chol(tXOmega %*% X + inverse(B)) )
    
    tmp = solve(L, tX %*% k + inverse(B) %*% b)
    alpha = solve(t(L),tmp) # up to here compiles, # here alpha is a matrix 
    
    result = mvrnormQuick(c(alpha), t(inverse(L))) # here is the problem as alpha is take to be a vector 
    
    
    return(result) 
    
  }
  
)

# Gamma sampler 
gamma_sampler_p <- nimbleFunction(
  name = 'gamma_sampler_p',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    # calcNodes <- model$getDependencies(target)
    X1 <- model$origData$X_p
    #N <- dim(X)[1]; M <- dim(X)[2]
    #k <- model$origData$k[1:N]
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    b <- model$b_p
    B <- model$B_p
    fixedIndexes <- 1
    ncov <-  if(!is.null(control$ncov)) control$ncov else 7
    d_bar <- 2 #ncov  
    d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
    # indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:ncov
    zindex <- control$zindex
    # J <- if(!is.null(control$J)) control$J else N
  },
  run = function() {
    
    
    k <-  model$y[model$z[zindex]==1] -0.5 
    X <- X1[model$z[zindex]==1,, drop=FALSE]
    N <- dim(X)[1]; M <- dim(X)[2]
    
    
    #  k <- model$y - (21/2) # works with K fixed
    Omega <- model$Omega_p[1:N]
    gamma <- model$gamma_p
    gamma_star <- gamma
    h_ratio <- 1
    M1 <- ncov+1
    
    if( runif(1,0, 1) < 0.33333) { # add
      
      if( (sum(gamma[1:M1]) - fixedIndexes) != ncov ){
        
        
        # Find zero covariates
        numOfZeroCov = ncov - ( sum(gamma[1:M1])  - fixedIndexes ) # correct to here
        zeroCov = nimNumeric(numOfZeroCov)
        
        i <- 1
        for (l in (fixedIndexes+1):length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariate_to_update = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))]
        
        
        # covariate_to_update =  sample_int(zeroCov[1:numOfZeroCov])
        
        #
        gamma_star[covariate_to_update] = 1
        h_ratio = (ncov -(sum(gamma[1:M1]) - fixedIndexes)) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) - 1 + ((ncov - d_bar) / (d_bar)) )
      }
      
    } else if( runif(1,0, 1) < .5){ # delete
      #
      if( ( sum(gamma[1:M1]) - fixedIndexes) != 0){
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        
        i = 1
        for (l in (fixedIndexes+1):length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariate_to_update = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariate_to_update] = 0
        
        h_ratio = (ncov - (sum(gamma[1:M1]) - fixedIndexes) + ((ncov - d_bar) / (d_bar)) ) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) + 1)
      }
      #
    }  else { # swap
      # 
      if( (sum(gamma[1:M1]) - fixedIndexes) != 0 & (sum(gamma[1:M1]) - fixedIndexes) != ncov ) {
        
        # Find zero covariates
        numOfZeroCov = ncov - (sum(gamma[1:M1]) - fixedIndexes)
        zeroCov = integer(numOfZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariates2_to_swap = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))] #sample_int(zeroCov)
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariates1_to_swap = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariates1_to_swap] = 0
        gamma_star[covariates2_to_swap] = 1
        
        h_ratio = 1
        
      }
      # 
    }
    
    model[[target]][1:M1] <<- gamma_star[1:M1] # update proposal
    
    # Compute log likelihood of proposed and current value
    L_gamma_star = dLgamma(gamma_star[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    L_gamma = dLgamma(gamma[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    h_ratio = h_ratio * exp(L_gamma_star - L_gamma) # MH ratio
    
    if( runif(1,0, 1) < h_ratio ){
      jump <- TRUE
    }  else jump <- FALSE
    
    # keep the model and mvSaved objects consistent
    if(jump) copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    else copy(from = mvSaved, to = model, row = 1,  nodes = target, logProb = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)


# sampler for detection probability (p)
PG_sampler_p <- nimbleFunction(
  name = 'PG_sampler_p',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X1 = model$origData$X_p
    # nsize = dim(X)[1] # assume that all observations have the same number of samples
    #  m <- dim(X)[2]
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:m
    zindex <- control$zindex
    # J <- if(!is.null(control$J)) control$J else nsize
  },
  run = function() {
    
    #  k <- model$y[model$z==1] -0.5 
    X = X1[model$z[zindex]==1,, drop=FALSE]
    nsize <- dim(X)[1]; m <- dim(X)[2]
    
    n <- rep(1, nsize)  #rep(model$K, N)
    # n <- rep(1, length(model$origData$k))
    #Kst <- model$K 
    #n <-rep(Kst,nsize)
    #  n <- model$K
    gamma <- model$gamma_p
    beta <- model$beta_p
    
    # resize beta and X 
    l <- 0
    index_present <- integer(length(gamma))
    
    for (i in 1:length(gamma)) {
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] <- i
        l <- l + 1
      }
    }
    
    X_gamma <- matrix(0, nsize,l)
    beta_gamma <- numeric(l)
    
    for (i in 1:l) {
      X_gamma[,i] <- X[, index_present[i]]
      beta_gamma[i] <- beta[index_present[i]]
    }
    
    model[[target]][1:nsize] <<- sample_Omega(X_gamma[1:nsize, 1:l], beta_gamma[1:l], n[1:nsize]) # this matches the result obtain from Alex's function
    
    # for (i in 1:nsize) {
    #   b = X[i, 1:m ] %*%  beta[1:m] 
    # model[[target]][i] <<-  pgr(n[i], b[1,1]) 
    # }
    
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)


beta_sampler_p <- nimbleFunction(
  name = 'beta_sampler_p',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    #calcNodes <- model$getDependencies(target) # all the betas
    X1 <- model$origData$X_p
    #N <- dim(X)[1]; M <- dim(X)[2]
    #k <- model$origData$k[1:N] # data needs to be in setup code  
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:M
    zindex <- control$zindex
    # J <- if(!is.null(control$J)) control$J else N
    b <- model$b_p
    B <- model$B_p
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    ##calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesDepStage <- model$isStoch(calcNodesDepStage)   ## should be made faster
    calcNodesDepStageDeterm <- calcNodesDepStage[!isStochCalcNodesDepStage]
    calcNodesDepStageStoch <- calcNodesDepStage[isStochCalcNodesDepStage]
    
  },
  run = function() {
    
    #Kst <- rep(1, N)  #rep(model$K, N)
    #Kst <- rep(model$K, N)
    #k <- model$y[1:N] - (Kst[1:N]/2) # works with K fixed
    # k <- model$origData$k[1:N] #model$y[model$z==1] - 0.5 # works with K fixed
    #  k <- model$y - (21/2) # works with K fixed
    
    k <-  model$y[model$z[zindex]==1] -0.5 
    X <- X1[model$z[zindex]==1,, drop=FALSE]
    N <- dim(X)[1]; M <- dim(X)[2]
    
    gamma <- model$gamma_p
    Omega <- model$Omega_p[1:N]
    
    # resize X, b and B
    l = 0
    index_present = integer(length(indexes_covariates))
    
    # gamma <- c(1,1,1,0,0)
    # indexes_covariates <- 1:5
    for (i in 1:length(indexes_covariates)) {
      
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] = i
        l = l + 1
      }
    }
    
    X_gamma2 = matrix(0, dim(X)[1], l)
    b_gamma = numeric(l)
    B_gamma = matrix(0, l,l)
    
    for (i in 1:l) {
      X_gamma2[,i] = X[, index_present[i]]
      b_gamma[i] = b[index_present[i]]
      for (j in 1:l) {
        B_gamma[i,j] = B[index_present[i], index_present[j]]
      }
    }
    
    # sample beta
    beta_new  <- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N])
    #model[[target]][1:l] <<- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N]) # when l not equal to length of gamma, I'm updating the wrong gamma
    
    for (i in 1:l) {
      # beta[index_present[i]] = beta_new[i]
      model[[target]][index_present[i]] <<- beta_new[i]  # need to test here if betas that are not updated remain the same or what should the new value of beta be
    }
    
    update_other_nodes <- model$calculateDiff(calcNodesDepStage)
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageStoch, logProbOnly = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)


# Samplers for persistence probability (phi)

# Gamma sampler 
gamma_sampler_phi <- nimbleFunction(
  name = 'gamma_sampler_phi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    # calcNodes <- model$getDependencies(target)
    X1 <- model$origData$X_phi
    #N <- dim(X)[1]; M <- dim(X)[2]
    #k <- model$origData$k[1:N]
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    b <- model$b
    B <- model$B
    fixedIndexes <- 1
    ncov <-  if(!is.null(control$ncov)) control$ncov else 7
    d_bar <- 2 #ncov  
    d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
    # indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:ncov
    #zindex <- control$zindex
    nsites <- control$nsites
    nyears <- control$nyears
    # J <- if(!is.null(control$J)) control$J else N
  },
  run = function() {
    
    
    d <- model$z[(nsites+1):(nsites*nyears)]
    c <- model$z[1:(nsites*(nyears-1))]
    k <- d[c==1] -0.5
    X <- X1[c==1,, drop=FALSE]
    N <- dim(X)[1]; M <- dim(X)[2]  
    
    # k <- model$z[(nsites+1):(nsites*nyears)][model$z[1:(nsites*(nyears-1))]==1] -0.5
    # X <- X1[model$z[1:(nsites*(nyears-1))]==1,, drop=FALSE]
    # N <- dim(X)[1]; M <- dim(X)[2]  
    # 
    #  k <- model$y - (21/2) # works with K fixed
    Omega <- model$Omega_phi[1:N]
    gamma <- model$gamma_phi
    gamma_star <- gamma
    h_ratio <- 1
    M1 <- ncov+1
    
    if( runif(1,0, 1) < 0.33333) { # add
      
      if( (sum(gamma[1:M1]) - fixedIndexes) != ncov ){
        
        
        # Find zero covariates
        numOfZeroCov = ncov - ( sum(gamma[1:M1])  - fixedIndexes ) # correct to here
        zeroCov = nimNumeric(numOfZeroCov)
        
        i <- 1
        for (l in (fixedIndexes+1):length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariate_to_update = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))]
        
        
        # covariate_to_update =  sample_int(zeroCov[1:numOfZeroCov])
        
        #
        gamma_star[covariate_to_update] = 1
        h_ratio = (ncov -(sum(gamma[1:M1]) - fixedIndexes)) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) - 1 + ((ncov - d_bar) / (d_bar)) )
      }
      
    } else if( runif(1,0, 1) < .5){ # delete
      #
      if( ( sum(gamma[1:M1]) - fixedIndexes) != 0){
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        
        i = 1
        for (l in (fixedIndexes+1):length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariate_to_update = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariate_to_update] = 0
        
        h_ratio = (ncov - (sum(gamma[1:M1]) - fixedIndexes) + ((ncov - d_bar) / (d_bar)) ) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) + 1)
      }
      #
    }  else { # swap
      # 
      if( (sum(gamma[1:M1]) - fixedIndexes) != 0 & (sum(gamma[1:M1]) - fixedIndexes) != ncov ) {
        
        # Find zero covariates
        numOfZeroCov = ncov - (sum(gamma[1:M1]) - fixedIndexes)
        zeroCov = integer(numOfZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariates2_to_swap = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))] #sample_int(zeroCov)
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariates1_to_swap = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariates1_to_swap] = 0
        gamma_star[covariates2_to_swap] = 1
        
        h_ratio = 1
        
      }
      # 
    }
    
    model[[target]][1:M1] <<- gamma_star[1:M1] # update proposal
    
    # Compute log likelihood of proposed and current value
    L_gamma_star = dLgamma(gamma_star[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    L_gamma = dLgamma(gamma[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    h_ratio = h_ratio * exp(L_gamma_star - L_gamma) # MH ratio
    
    if( runif(1,0, 1) < h_ratio ){
      jump <- TRUE
    }  else jump <- FALSE
    
    # keep the model and mvSaved objects consistent
    if(jump) copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    else copy(from = mvSaved, to = model, row = 1,  nodes = target, logProb = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)


PG_sampler_phi <- nimbleFunction(
  name = 'PG_sampler_phi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X1 = model$origData$X_phi
    # nsize = dim(X)[1] # assume that all observations have the same number of samples
    #  m <- dim(X)[2]
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:m
    nsites <- control$nsites
    nyears <- control$nyears
    # zindex <- control$zindex
    # J <- if(!is.null(control$J)) control$J else nsize
  },
  run = function() {
    
    #  k <- model$y[model$z==1] -0.5 
    
    c <- model$z[1:(nsites*(nyears-1))]
    X <- X1[c==1,, drop=FALSE]
    nsize <- dim(X)[1]
    
    n <- rep(1, nsize)  #rep(model$K, N)
    # n <- rep(1, length(model$origData$k))
    #Kst <- model$K 
    #n <-rep(Kst,nsize)
    #  n <- model$K
    gamma <- model$gamma_phi
    beta <- model$beta_phi
    
    # resize beta and X 
    l <- 0
    index_present <- integer(length(gamma))
    
    for (i in 1:length(gamma)) {
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] <- i
        l <- l + 1
      }
    }
    
    X_gamma <- matrix(0, nsize,l)
    beta_gamma <- numeric(l)
    
    for (i in 1:l) {
      X_gamma[,i] <- X[, index_present[i]]
      beta_gamma[i] <- beta[index_present[i]]
    }
    
    model[[target]][1:nsize] <<- sample_Omega(X_gamma[1:nsize, 1:l], beta_gamma[1:l], n[1:nsize]) # this matches the result obtain from Alex's function
    
    # for (i in 1:nsize) {
    #   b = X[i, 1:m ] %*%  beta[1:m] 
    # model[[target]][i] <<-  pgr(n[i], b[1,1]) 
    # }
    
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)

beta_sampler_phi <- nimbleFunction(
  name = 'beta_sampler_phi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    #calcNodes <- model$getDependencies(target) # all the betas
    X1 <- model$origData$X_phi
    #N <- dim(X)[1]; M <- dim(X)[2]
    #k <- model$origData$k[1:N] # data needs to be in setup code  
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:M
    nsites <- control$nsites
    nyears <- control$nyears
    #zindex <- control$zindex
    # J <- if(!is.null(control$J)) control$J else N
    b <- model$b
    B <- model$B
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    ##calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesDepStage <- model$isStoch(calcNodesDepStage)   ## should be made faster
    calcNodesDepStageDeterm <- calcNodesDepStage[!isStochCalcNodesDepStage]
    calcNodesDepStageStoch <- calcNodesDepStage[isStochCalcNodesDepStage]
    
  },
  run = function() {
    
    
    d <- model$z[(nsites+1):(nsites*nyears)]
    c <- model$z[1:(nsites*(nyears-1))]
    k <- d[c==1] -0.5
    X <- X1[c==1,, drop=FALSE]
    N <- dim(X)[1]; M <- dim(X)[2]  
    
    gamma <- model$gamma_phi
    Omega <- model$Omega_phi[1:N]
    
    # resize X, b and B
    l = 0
    index_present = integer(length(indexes_covariates))
    
    # gamma <- c(1,1,1,0,0)
    # indexes_covariates <- 1:5
    for (i in 1:length(indexes_covariates)) {
      
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] = i
        l = l + 1
      }
    }
    
    X_gamma2 = matrix(0, dim(X)[1], l)
    b_gamma = numeric(l)
    B_gamma = matrix(0, l,l)
    
    for (i in 1:l) {
      X_gamma2[,i] = X[, index_present[i]]
      b_gamma[i] = b[index_present[i]]
      for (j in 1:l) {
        B_gamma[i,j] = B[index_present[i], index_present[j]]
      }
    }
    
    # sample beta
    beta_new  <- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N])
    #model[[target]][1:l] <<- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N]) # when l not equal to length of gamma, I'm updating the wrong gamma
    
    for (i in 1:l) {
      # beta[index_present[i]] = beta_new[i]
      model[[target]][index_present[i]] <<- beta_new[i]  # need to test here if betas that are not updated remain the same or what should the new value of beta be
    }
    
    update_other_nodes <- model$calculateDiff(calcNodesDepStage)
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageStoch, logProbOnly = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)

# Colonization

# Gamma sampler 
gamma_sampler_eta <- nimbleFunction(
  name = 'gamma_sampler_eta',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    # calcNodes <- model$getDependencies(target)
    X1 <- model$origData$X_eta
    #N <- dim(X)[1]; M <- dim(X)[2]
    #k <- model$origData$k[1:N]
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    b <- model$b
    B <- model$B
    fixedIndexes <- 1
    ncov <-  if(!is.null(control$ncov)) control$ncov else 7
    d_bar <- 2 #ncov  
    d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
    # indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:ncov
    #zindex <- control$zindex
    nsites <- control$nsites
    nyears <- control$nyears
    # J <- if(!is.null(control$J)) control$J else N
  },
  run = function() {
    
    
    d <- model$z[(nsites+1):(nsites*nyears)]
    c <- model$z[1:(nsites*(nyears-1))]
    k <- d[c==0] -0.5
    X <- X1[c==0,, drop=FALSE]
    N <- dim(X)[1]; M <- dim(X)[2]  
    
    # k <- model$z[(nsites+1):(nsites*nyears)][model$z[1:(nsites*(nyears-1))]==1] -0.5
    # X <- X1[model$z[1:(nsites*(nyears-1))]==1,, drop=FALSE]
    # N <- dim(X)[1]; M <- dim(X)[2]  
    # 
    #  k <- model$y - (21/2) # works with K fixed
    Omega <- model$Omega_eta[1:N]
    gamma <- model$gamma_eta
    gamma_star <- gamma
    h_ratio <- 1
    M1 <- ncov+1
    
    if( runif(1,0, 1) < 0.33333) { # add
      
      if( (sum(gamma[1:M1]) - fixedIndexes) != ncov ){
        
        
        # Find zero covariates
        numOfZeroCov = ncov - ( sum(gamma[1:M1])  - fixedIndexes ) # correct to here
        zeroCov = nimNumeric(numOfZeroCov)
        
        i <- 1
        for (l in (fixedIndexes+1):length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariate_to_update = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))]
        
        
        # covariate_to_update =  sample_int(zeroCov[1:numOfZeroCov])
        
        #
        gamma_star[covariate_to_update] = 1
        h_ratio = (ncov -(sum(gamma[1:M1]) - fixedIndexes)) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) - 1 + ((ncov - d_bar) / (d_bar)) )
      }
      
    } else if( runif(1,0, 1) < .5){ # delete
      #
      if( ( sum(gamma[1:M1]) - fixedIndexes) != 0){
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        
        i = 1
        for (l in (fixedIndexes+1):length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariate_to_update = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariate_to_update] = 0
        
        h_ratio = (ncov - (sum(gamma[1:M1]) - fixedIndexes) + ((ncov - d_bar) / (d_bar)) ) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) + 1)
      }
      #
    }  else { # swap
      # 
      if( (sum(gamma[1:M1]) - fixedIndexes) != 0 & (sum(gamma[1:M1]) - fixedIndexes) != ncov ) {
        
        # Find zero covariates
        numOfZeroCov = ncov - (sum(gamma[1:M1]) - fixedIndexes)
        zeroCov = integer(numOfZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariates2_to_swap = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))] #sample_int(zeroCov)
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariates1_to_swap = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariates1_to_swap] = 0
        gamma_star[covariates2_to_swap] = 1
        
        h_ratio = 1
        
      }
      # 
    }
    
    model[[target]][1:M1] <<- gamma_star[1:M1] # update proposal
    
    # Compute log likelihood of proposed and current value
    L_gamma_star = dLgamma(gamma_star[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    L_gamma = dLgamma(gamma[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    h_ratio = h_ratio * exp(L_gamma_star - L_gamma) # MH ratio
    
    if( runif(1,0, 1) < h_ratio ){
      jump <- TRUE
    }  else jump <- FALSE
    
    # keep the model and mvSaved objects consistent
    if(jump) copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    else copy(from = mvSaved, to = model, row = 1,  nodes = target, logProb = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)


PG_sampler_eta <- nimbleFunction(
  name = 'PG_sampler_eta',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X1 = model$origData$X_eta
    # nsize = dim(X)[1] # assume that all observations have the same number of samples
    #  m <- dim(X)[2]
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:m
    nsites <- control$nsites
    nyears <- control$nyears
    # zindex <- control$zindex
    # J <- if(!is.null(control$J)) control$J else nsize
  },
  run = function() {
    
    #  k <- model$y[model$z==1] -0.5 
    
    c <- model$z[1:(nsites*(nyears-1))]
    X <- X1[c==0,, drop=FALSE]
    nsize <- dim(X)[1]
    
    n <- rep(1, nsize)  #rep(model$K, N)
    # n <- rep(1, length(model$origData$k))
    #Kst <- model$K 
    #n <-rep(Kst,nsize)
    #  n <- model$K
    gamma <- model$gamma_eta
    beta <- model$beta_eta
    
    # resize beta and X 
    l <- 0
    index_present <- integer(length(gamma))
    
    for (i in 1:length(gamma)) {
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] <- i
        l <- l + 1
      }
    }
    
    X_gamma <- matrix(0, nsize,l)
    beta_gamma <- numeric(l)
    
    for (i in 1:l) {
      X_gamma[,i] <- X[, index_present[i]]
      beta_gamma[i] <- beta[index_present[i]]
    }
    
    model[[target]][1:nsize] <<- sample_Omega(X_gamma[1:nsize, 1:l], beta_gamma[1:l], n[1:nsize]) # this matches the result obtain from Alex's function
    
    # for (i in 1:nsize) {
    #   b = X[i, 1:m ] %*%  beta[1:m] 
    # model[[target]][i] <<-  pgr(n[i], b[1,1]) 
    # }
    
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)


beta_sampler_eta <- nimbleFunction(
  name = 'beta_sampler_eta',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    #calcNodes <- model$getDependencies(target) # all the betas
    X1 <- model$origData$X_eta
    #N <- dim(X)[1]; M <- dim(X)[2]
    #k <- model$origData$k[1:N] # data needs to be in setup code  
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:M
    nsites <- control$nsites
    nyears <- control$nyears
    #zindex <- control$zindex
    # J <- if(!is.null(control$J)) control$J else N
    b <- model$b
    B <- model$B
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    ##calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesDepStage <- model$isStoch(calcNodesDepStage)   ## should be made faster
    calcNodesDepStageDeterm <- calcNodesDepStage[!isStochCalcNodesDepStage]
    calcNodesDepStageStoch <- calcNodesDepStage[isStochCalcNodesDepStage]
    
  },
  run = function() {
    
    
    d <- model$z[(nsites+1):(nsites*nyears)]
    c <- model$z[1:(nsites*(nyears-1))]
    k <- d[c==0] -0.5
    X <- X1[c==0,, drop=FALSE]
    N <- dim(X)[1]; M <- dim(X)[2]  
    
    
    # nsites <- 100
    # nyears <- 10
    # k <- model$z[(nsites+1):(nsites*nyears)][model$z[1:(nsites*(nyears-1))]==1] -0.5
    # X <- X1[model$z[1:(nsites*(nyears-1))]==1,, drop=FALSE]
    # N <- dim(X)[1]; M <- dim(X)[2]  
    # 
    
    gamma <- model$gamma_eta
    Omega <- model$Omega_eta[1:N]
    
    # resize X, b and B
    l = 0
    index_present = integer(length(indexes_covariates))
    
    # gamma <- c(1,1,1,0,0)
    # indexes_covariates <- 1:5
    for (i in 1:length(indexes_covariates)) {
      
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] = i
        l = l + 1
      }
    }
    
    X_gamma2 = matrix(0, dim(X)[1], l)
    b_gamma = numeric(l)
    B_gamma = matrix(0, l,l)
    
    for (i in 1:l) {
      X_gamma2[,i] = X[, index_present[i]]
      b_gamma[i] = b[index_present[i]]
      for (j in 1:l) {
        B_gamma[i,j] = B[index_present[i], index_present[j]]
      }
    }
    
    # sample beta
    beta_new  <- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N])
    #model[[target]][1:l] <<- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N]) # when l not equal to length of gamma, I'm updating the wrong gamma
    
    for (i in 1:l) {
      # beta[index_present[i]] = beta_new[i]
      model[[target]][index_present[i]] <<- beta_new[i]  # need to test here if betas that are not updated remain the same or what should the new value of beta be
    }
    
    update_other_nodes <- model$calculateDiff(calcNodesDepStage)
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageStoch, logProbOnly = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)


# Initial occupancy (psi)

# Gamma sampler 
gamma_sampler_psi <- nimbleFunction(
  name = 'gamma_sampler_psi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    # calcNodes <- model$getDependencies(target)
    X <- model$origData$X_psi
    N <- dim(X)[1]; M <- dim(X)[2]
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    b <- model$b_psi
    B <- model$B_psi
    fixedIndexes <- 1
    ncov <-  if(!is.null(control$ncov)) control$ncov else 7
    d_bar <- 2 #ncov  
    d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
    # indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:ncov
    # J <- if(!is.null(control$J)) control$J else N
  },
  run = function() {
    
    Kst <- rep(1, N)  #rep(model$K, N)
    # Kst <- model$K
    k <- model$z[1:N] - (Kst[1:N]/2) 
    #  k <- model$y - (21/2) # works with K fixed
    Omega <- model$Omega_psi
    gamma <- model$gamma_psi
    gamma_star <- gamma
    h_ratio <- 1
    M1 <- ncov+1
    
    if( runif(1,0, 1) < 0.33333) { # add
      
      if( (sum(gamma[1:M1]) - fixedIndexes) != ncov ){
        
        
        # Find zero covariates
        numOfZeroCov = ncov - ( sum(gamma[1:M1])  - fixedIndexes ) # correct to here
        zeroCov = nimNumeric(numOfZeroCov)
        
        i <- 1
        for (l in (fixedIndexes+1):length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariate_to_update = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))]
        
        
        # covariate_to_update =  sample_int(zeroCov[1:numOfZeroCov])
        
        #
        gamma_star[covariate_to_update] = 1
        h_ratio = (ncov -(sum(gamma[1:M1]) - fixedIndexes)) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) - 1 + ((ncov - d_bar) / (d_bar)) )
      }
      
    } else if( runif(1,0, 1) < .5){ # delete
      #
      if( ( sum(gamma[1:M1]) - fixedIndexes) != 0){
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        
        i = 1
        for (l in (fixedIndexes+1):length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariate_to_update = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariate_to_update] = 0
        
        h_ratio = (ncov - (sum(gamma[1:M1]) - fixedIndexes) + ((ncov - d_bar) / (d_bar)) ) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) + 1)
      }
      #
    }  else { # swap
      # 
      if( (sum(gamma[1:M1]) - fixedIndexes) != 0 & (sum(gamma[1:M1]) - fixedIndexes) != ncov ) {
        
        # Find zero covariates
        numOfZeroCov = ncov - (sum(gamma[1:M1]) - fixedIndexes)
        zeroCov = integer(numOfZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariates2_to_swap = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))] #sample_int(zeroCov)
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariates1_to_swap = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariates1_to_swap] = 0
        gamma_star[covariates2_to_swap] = 1
        
        h_ratio = 1
        
      }
      # 
    }
    
    model[[target]][1:M1] <<- gamma_star[1:M1] # update proposal
    
    # Compute log likelihood of proposed and current value
    L_gamma_star = dLgamma(gamma_star[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    L_gamma = dLgamma(gamma[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    h_ratio = h_ratio * exp(L_gamma_star - L_gamma) # MH ratio
    
    if( runif(1,0, 1) < h_ratio ){
      jump <- TRUE
    }  else jump <- FALSE
    
    # keep the model and mvSaved objects consistent
    if(jump) copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    else copy(from = mvSaved, to = model, row = 1,  nodes = target, logProb = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)

PG_sampler_psi <- nimbleFunction(
  name = 'PG_sampler_psi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X = model$origData$X_psi
    nsize = dim(X)[1] # assume that all observations have the same number of samples
    m <- dim(X)[2]
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:m
    # J <- if(!is.null(control$J)) control$J else nsize
  },
  run = function() {
    
    n <- rep(1, nsize)  #rep(model$K, N)
    #Kst <- model$K 
    #n <-rep(Kst,nsize)
    #  n <- model$K
    gamma <- model$gamma_psi
    beta <- model$beta_psi
    
    # resize beta and X 
    l <- 0
    index_present <- integer(length(gamma))
    
    for (i in 1:length(gamma)) {
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] <- i
        l <- l + 1
      }
    }
    
    X_gamma <- matrix(0, nsize,l)
    beta_gamma <- numeric(l)
    
    for (i in 1:l) {
      X_gamma[,i] <- X[, index_present[i]]
      beta_gamma[i] <- beta[index_present[i]]
    }
    
    model[[target]][1:nsize] <<- sample_Omega(X_gamma[1:nsize, 1:l], beta_gamma[1:l], n[1:nsize]) # this matches the result obtain from Alex's function
    
    # for (i in 1:nsize) {
    #   b = X[i, 1:m ] %*%  beta[1:m] 
    # model[[target]][i] <<-  pgr(n[i], b[1,1]) 
    # }
    
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)


beta_sampler_psi <- nimbleFunction(
  name = 'beta_sampler_psi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    #calcNodes <- model$getDependencies(target) # all the betas
    X <- model$origData$X_psi
    N <- dim(X)[1]; M <- dim(X)[2]
    #indexes_covariates <- c(1,2,3,4,5,6,7,8)  #model$origData$indexes_covariates
    #J=  c(6 , 6 , 8 , 6 , 8,  8 , 8,  8,  8,  8,  8 , 8,  8 , 8,  8 , 8,  8,  8,  8,  8,  8,  8,  8,  6,  8,  8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:M
    # J <- if(!is.null(control$J)) control$J else N
    b <- model$b_psi
    B <- model$B_psi
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    ##calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesDepStage <- model$isStoch(calcNodesDepStage)   ## should be made faster
    calcNodesDepStageDeterm <- calcNodesDepStage[!isStochCalcNodesDepStage]
    calcNodesDepStageStoch <- calcNodesDepStage[isStochCalcNodesDepStage]
    
  },
  run = function() {
    
    Kst <- rep(1, N)  #rep(model$K, N)
    #Kst <- rep(model$K, N)
    # Kst <- model$K
    k <- model$z[1:N] - (Kst[1:N]/2) # works with K fixed
    #  k <- model$y - (21/2) # works with K fixed
    gamma <- model$gamma_psi
    Omega <- model$Omega_psi
    
    # resize X, b and B
    l = 0
    index_present = integer(length(indexes_covariates))
    
    # gamma <- c(1,1,1,0,0)
    # indexes_covariates <- 1:5
    for (i in 1:length(indexes_covariates)) {
      
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] = i
        l = l + 1
      }
    }
    
    X_gamma2 = matrix(0, dim(X)[1], l)
    b_gamma = numeric(l)
    B_gamma = matrix(0, l,l)
    
    for (i in 1:l) {
      X_gamma2[,i] = X[, index_present[i]]
      b_gamma[i] = b[index_present[i]]
      for (j in 1:l) {
        B_gamma[i,j] = B[index_present[i], index_present[j]]
      }
    }
    
    # sample beta
    beta_new  <- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N])
    #model[[target]][1:l] <<- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N]) # when l not equal to length of gamma, I'm updating the wrong gamma
    
    for (i in 1:l) {
      # beta[index_present[i]] = beta_new[i]
      model[[target]][index_present[i]] <<- beta_new[i]  # need to test here if betas that are not updated remain the same or what should the new value of beta be
    }
    
    update_other_nodes <- model$calculateDiff(calcNodesDepStage)
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageStoch, logProbOnly = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)


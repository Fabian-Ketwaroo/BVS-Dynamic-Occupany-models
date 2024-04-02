
rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)
library(tidyverse)

sumSim <- function(samples, true_value, PIP, nsim, indexes_covariates ){
  
  RB = CV= matrix(0, length(true_value), nsim)
  
  for (i in 1:nsim) {
    
    L <- length(indexes_covariates)
    index_present <- logical(L)
    
    for (l in 1:L) {
      if(any(indexes_covariates[l] == which(PIP[,1,i]>=0.5))){
        index_present[l] <- TRUE
      } else index_present[l] <- FALSE
      
    }
    
    beta_gamma <- samples[index_present,,i]
    true_beta_gamma = true_value[index_present]
    
    RB[index_present,i] =  (beta_gamma[,1] - true_beta_gamma)/true_beta_gamma
    CV[index_present,i] =   beta_gamma[,3]/ abs(beta_gamma[,1])  
  }
  
  
  
  
  
  out = list( RB = RB, CV = CV)
  return(out)
  
}


nsim = 10
ncov = 10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/indexes_covariates.Rdata")

#================================================================================================================================================
# p
#================================================================================================================================================

#######################################################################################################
# High
#######################################################################################################


# high case
# true_betas_psi <- c(1.8, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
# true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
# true_betas_eta <- c(1, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
true_betas_p <- c(1.8, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept

# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_beta_p.Rdata")


s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p1_M500_high_RB =  s$RB[which(s$RB!=0)]
p1_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_beta_p.Rdata")


s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p2_M500_high_RB =  s$RB[which(s$RB!=0)]
p2_M500_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p3_M500_high_RB =  s$RB[which(s$RB!=0)]
p3_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p4_M500_high_RB =  s$RB[which(s$RB!=0)]
p4_M500_high_CV =  s$CV[which(s$RB!=0)]

#========================================================================
# M= 200
#======================================================================


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p5_M200_high_RB =  s$RB[which(s$RB!=0)]
p5_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p6_M200_high_RB =  s$RB[which(s$RB!=0)]
p6_M200_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p7_M200_high_RB =  s$RB[which(s$RB!=0)]
p7_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p8_M200_high_RB =  s$RB[which(s$RB!=0)]
p8_M200_high_CV =  s$CV[which(s$RB!=0)]


#######################################################################################################
#Low
#######################################################################################################

# true_betas_psi <- c(-0.25, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
# true_betas_phi <- c(-1, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
# true_betas_eta <- c(-1.5, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
true_betas_p <- c(-0.25, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept



# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_beta_p.Rdata")


s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p1_M500_low_RB =  s$RB[which(s$RB!=0)]
p1_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_beta_p.Rdata")


s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p2_M500_low_RB =  s$RB[which(s$RB!=0)]
p2_M500_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p3_M500_low_RB =  s$RB[which(s$RB!=0)]
p3_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p4_M500_low_RB =  s$RB[which(s$RB!=0)]
p4_M500_low_CV =  s$CV[which(s$RB!=0)]




# M= 200


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_beta_p.Rdata")


s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p5_M200_low_RB =  s$RB[which(s$RB!=0)]
p5_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p6_M200_low_RB =  s$RB[which(s$RB!=0)]
p6_M200_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p7_M200_low_RB =  s$RB[which(s$RB!=0)]
p7_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_gamma_p.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_beta_p.Rdata")

s = sumSim(beta_p, true_betas_p, gamma_p, nsim, indexes_covariates)
p8_M200_low_RB =  s$RB[which(s$RB!=0)]
p8_M200_low_CV =  s$CV[which(s$RB!=0)]

p8_M200_low_CV = p8_M200_low_CV[which(p8_M200_low_RB==-Inf|p8_M200_low_RB==Inf)]
boxplot(p8_M200_low_CV)

############################################################################################################

# Violin Plot

############################################################################################################



RB = c(p1_M500_high_RB, p1_M500_low_RB, p2_M500_high_RB, p2_M500_low_RB, p3_M500_high_RB, p3_M500_low_RB, p4_M500_high_RB, p4_M500_low_RB, p5_M200_high_RB, p5_M200_low_RB, p6_M200_high_RB, p6_M200_low_RB, p7_M200_high_RB, p7_M200_low_RB, p8_M200_high_RB, p8_M200_low_RB)
CV = c(p1_M500_high_CV, p1_M500_low_CV, p2_M500_high_CV, p2_M500_low_CV, p3_M500_high_CV, p3_M500_low_CV, p4_M500_high_CV, p4_M500_low_CV, p5_M200_high_CV, p5_M200_low_CV, p6_M200_high_CV, p6_M200_low_CV, p7_M200_high_CV, p7_M200_low_CV, p8_M200_high_CV, p8_M200_low_CV)
classes = factor(c(rep(1, length(c(p1_M500_high_RB, p1_M500_low_RB))), rep(2, length(c(p2_M500_high_RB, p2_M500_low_RB))), rep(3, length(c(p3_M500_high_RB, p3_M500_low_RB))), rep(4, length(c(p4_M500_high_RB, p4_M500_low_RB))), rep(5, length(c(p5_M200_high_RB, p5_M200_low_RB))), rep(6, length(c(p6_M200_high_RB, p6_M200_low_RB))), rep(7, length(c(p7_M200_high_RB, p7_M200_low_RB))), rep(8, length(c(p8_M200_high_RB, p8_M200_low_RB)))), levels = 1:8)
Level = factor(c( rep("High", length(p1_M500_high_RB)), rep("Low", length(p1_M500_low_RB)), rep("High", length(p2_M500_high_RB)), rep("Low", length(p2_M500_low_RB)), rep("High", length(p3_M500_high_RB)), rep("Low", length(p3_M500_low_RB)), rep("High", length(p4_M500_high_RB)), rep("Low", length(p4_M500_low_RB)),   rep("High", length(p5_M200_high_RB)), rep("Low", length(p5_M200_low_RB)), rep("High", length(p6_M200_high_RB)), rep("Low", length(p6_M200_low_RB)), rep("High", length(p7_M200_high_RB)), rep("Low", length(p7_M200_low_RB)), rep("High", length(p8_M200_high_RB)), rep("Low", length(p8_M200_low_RB))  ),levels  =c("High", "Low") )
length(RB) == length(classes) 
length(classes) == length(Level)

results = data.frame(RB=RB, CV = CV, classes = classes, Levels = Level)
head(results)

ggplot(results, aes(x = classes, y = RB, fill = Levels)) + geom_violin(trim = TRUE)   + stat_summary(fun =mean,geom="point", shape=16, size=2,  position = position_dodge2(width = 0.9,  preserve = "single")) + theme_classic() + ylab("Relative Bias") + xlab("Cases") +  theme(legend.position = "right", legend.title = element_text(size =14), legend.text = element_text(size =14),  axis.text = element_text(size=14), axis.title = element_text(size=14)) 

ggplot(results, aes(x = classes, y = CV, fill = Levels)) + geom_violin(trim = TRUE)   + stat_summary(fun =mean,geom="point", shape=16, size=2,  position = position_dodge2(width = 0.9,  preserve = "single")) + theme_classic() + ylab("Coefficient of Variation") + xlab("Cases") +  theme(legend.position = "right", legend.title = element_text(size =14), legend.text = element_text(size =14),  axis.text = element_text(size=14), axis.title = element_text(size=14)) 

#=======================================================================================================================================
#phi
#======================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)
library(tidyverse)

sumSim <- function(samples, true_value, PIP, nsim, indexes_covariates ){
  
  RB = CV= matrix(0, length(true_value), nsim)
  
  for (i in 1:nsim) {
    
    L <- length(indexes_covariates)
    index_present <- logical(L)
    
    for (l in 1:L) {
      if(any(indexes_covariates[l] == which(PIP[,1,i]>=0.5))){
        index_present[l] <- TRUE
      } else index_present[l] <- FALSE
      
    }
    
    beta_gamma <- samples[index_present,,i]
    true_beta_gamma = true_value[index_present]
    
    RB[index_present,i] =  (beta_gamma[,1] - true_beta_gamma)/true_beta_gamma
    CV[index_present,i] =   beta_gamma[,3]/ abs(beta_gamma[,1])  
  }
  
  
  
  
  
  out = list( RB = RB, CV = CV)
  return(out)
  
}


nsim = 10
ncov = 10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/indexes_covariates.Rdata")


#######################################################################################################
# High
#######################################################################################################


# high case
# true_betas_psi <- c(1.8, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
#true_betas_p= true_betas_phi
# true_betas_eta <- c(1, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
#true_betas_p <- c(1.8, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept

# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_beta_phi.Rdata")


s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p1_M500_high_RB =  s$RB[which(s$RB!=0)]
p1_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_beta_phi.Rdata")


s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p2_M500_high_RB =  s$RB[which(s$RB!=0)]
p2_M500_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p3_M500_high_RB =  s$RB[which(s$RB!=0)]
p3_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p4_M500_high_RB =  s$RB[which(s$RB!=0)]
p4_M500_high_CV =  s$CV[which(s$RB!=0)]

#========================================================================
# M= 200
#======================================================================


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p5_M200_high_RB =  s$RB[which(s$RB!=0)]
p5_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p6_M200_high_RB =  s$RB[which(s$RB!=0)]
p6_M200_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p7_M200_high_RB =  s$RB[which(s$RB!=0)]
p7_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p8_M200_high_RB =  s$RB[which(s$RB!=0)]
p8_M200_high_CV =  s$CV[which(s$RB!=0)]


#######################################################################################################
#Low
#######################################################################################################

# true_betas_psi <- c(-0.25, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
 true_betas_phi <- c(-1, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
# true_betas_eta <- c(-1.5, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
#true_betas_p <- c(-0.25, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept



# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_beta_phi.Rdata")


s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p1_M500_low_RB =  s$RB[which(s$RB!=0)]
p1_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_beta_phi.Rdata")


s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p2_M500_low_RB =  s$RB[which(s$RB!=0)]
p2_M500_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p3_M500_low_RB =  s$RB[which(s$RB!=0)]
p3_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p4_M500_low_RB =  s$RB[which(s$RB!=0)]
p4_M500_low_CV =  s$CV[which(s$RB!=0)]




# M= 200


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_beta_phi.Rdata")


s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p5_M200_low_RB =  s$RB[which(s$RB!=0)]
p5_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p6_M200_low_RB =  s$RB[which(s$RB!=0)]
p6_M200_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p7_M200_low_RB =  s$RB[which(s$RB!=0)]
p7_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_gamma_phi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_beta_phi.Rdata")

s = sumSim(beta_phi, true_betas_phi, gamma_phi, nsim, indexes_covariates)
p8_M200_low_RB =  s$RB[which(s$RB!=0)]
p8_M200_low_CV =  s$CV[which(s$RB!=0)]


############################################################################################################

# Violin Plot

############################################################################################################



RB = c(p1_M500_high_RB, p1_M500_low_RB, p2_M500_high_RB, p2_M500_low_RB, p3_M500_high_RB, p3_M500_low_RB, p4_M500_high_RB, p4_M500_low_RB, p5_M200_high_RB, p5_M200_low_RB, p6_M200_high_RB, p6_M200_low_RB, p7_M200_high_RB, p7_M200_low_RB, p8_M200_high_RB, p8_M200_low_RB)
classes = factor(c(rep(1, length(c(p1_M500_high_RB, p1_M500_low_RB))), rep(2, length(c(p2_M500_high_RB, p2_M500_low_RB))), rep(3, length(c(p3_M500_high_RB, p3_M500_low_RB))), rep(4, length(c(p4_M500_high_RB, p4_M500_low_RB))), rep(5, length(c(p5_M200_high_RB, p5_M200_low_RB))), rep(6, length(c(p6_M200_high_RB, p6_M200_low_RB))), rep(7, length(c(p7_M200_high_RB, p7_M200_low_RB))), rep(8, length(c(p8_M200_high_RB, p8_M200_low_RB)))), levels = 1:8)
Level = factor(c( rep("High", length(p1_M500_high_RB)), rep("Low", length(p1_M500_low_RB)), rep("High", length(p2_M500_high_RB)), rep("Low", length(p2_M500_low_RB)), rep("High", length(p3_M500_high_RB)), rep("Low", length(p3_M500_low_RB)), rep("High", length(p4_M500_high_RB)), rep("Low", length(p4_M500_low_RB)),   rep("High", length(p5_M200_high_RB)), rep("Low", length(p5_M200_low_RB)), rep("High", length(p6_M200_high_RB)), rep("Low", length(p6_M200_low_RB)), rep("High", length(p7_M200_high_RB)), rep("Low", length(p7_M200_low_RB)), rep("High", length(p8_M200_high_RB)), rep("Low", length(p8_M200_low_RB))  ),levels  =c("High", "Low") )
length(RB) == length(classes) 
length(classes) == length(Level)

results = data.frame(RB=RB, classes = classes, Levels = Level)
head(results)

ggplot(results, aes(x = classes, y = RB, fill = Levels)) + geom_violin(trim = TRUE)   + stat_summary(fun =mean,geom="point", shape=16, size=2,  position = position_dodge2(width = 0.9,  preserve = "single")) + theme_classic() + ylab("Relative Bias") + xlab("Cases") +  theme(legend.position = "right", legend.title = element_text(size =14), legend.text = element_text(size =14),  axis.text = element_text(size=14), axis.title = element_text(size=14)) 


#=======================================================================================================================================
#psi
#======================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)
library(tidyverse)

# samples = beta_psi
# true_value = true_betas_psi
# PIP = gamma_psi

sumSim <- function(samples, true_value, PIP, nsim, indexes_covariates ){
  
  RB = CV= matrix(0, length(true_value), nsim)
  
  for (i in 1:nsim) {
    
    L <- length(indexes_covariates)
    index_present <- logical(L)
    
    for (l in 1:L) {
      if(any(indexes_covariates[l] == which(PIP[,1,i]>=0.5))){
        index_present[l] <- TRUE
      } else index_present[l] <- FALSE
      
    }
    
    beta_gamma <- samples[index_present,,i]
    true_beta_gamma = true_value[index_present]
    
    if(is.matrix(beta_gamma)){
      RB[index_present,i] =  (beta_gamma[,1] - true_beta_gamma)/true_beta_gamma
      CV[index_present,i] =   beta_gamma[,3]/ abs(beta_gamma[,1]) 
    } else {
       RB[index_present,i] =  (beta_gamma[1] - true_beta_gamma)/true_beta_gamma
      CV[index_present,i] =   beta_gamma[3]/ abs(beta_gamma[1]) 
    }
    
  }
  
  
  
  
  
  out = list( RB = RB, CV = CV)
  return(out)
  
}


nsim = 10
ncov = 10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/indexes_covariates.Rdata")


#######################################################################################################
# High
#######################################################################################################


# high case
 true_betas_psi <- c(1.8, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
#true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
#true_betas_p= true_betas_phi
# true_betas_eta <- c(1, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
#true_betas_p <- c(1.8, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept

# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_beta_psi.Rdata")


s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p1_M500_high_RB =  s$RB[which(s$RB!=0)]
p1_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_beta_psi.Rdata")


s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p2_M500_high_RB =  s$RB[which(s$RB!=0)]
p2_M500_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p3_M500_high_RB =  s$RB[which(s$RB!=0)]
p3_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p4_M500_high_RB =  s$RB[which(s$RB!=0)]
p4_M500_high_CV =  s$CV[which(s$RB!=0)]

#========================================================================
# M= 200
#======================================================================


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p5_M200_high_RB =  s$RB[which(s$RB!=0)]
p5_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p6_M200_high_RB =  s$RB[which(s$RB!=0)]
p6_M200_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p7_M200_high_RB =  s$RB[which(s$RB!=0)]
p7_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p8_M200_high_RB =  s$RB[which(s$RB!=0)]
p8_M200_high_CV =  s$CV[which(s$RB!=0)]


#######################################################################################################
#Low
#######################################################################################################

 true_betas_psi <- c(-0.25, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
# true_betas_phi <- c(-1, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
# true_betas_eta <- c(-1.5, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
#true_betas_p <- c(-0.25, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept



# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p1_M500_low_RB =  s$RB[which(s$RB!=0)]
p1_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_beta_psi.Rdata")


s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p2_M500_low_RB =  s$RB[which(s$RB!=0)]
p2_M500_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p3_M500_low_RB =  s$RB[which(s$RB!=0)]
p3_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p4_M500_low_RB =  s$RB[which(s$RB!=0)]
p4_M500_low_CV =  s$CV[which(s$RB!=0)]




# M= 200


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p5_M200_low_RB =  s$RB[which(s$RB!=0)]
p5_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p6_M200_low_RB =  s$RB[which(s$RB!=0)]
p6_M200_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p7_M200_low_RB =  s$RB[which(s$RB!=0)]
p7_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_gamma_psi.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_beta_psi.Rdata")

s = sumSim(beta_psi, true_betas_psi, gamma_psi, nsim, indexes_covariates)
p8_M200_low_RB =  s$RB[which(s$RB!=0)]
p8_M200_low_CV =  s$CV[which(s$RB!=0)]


############################################################################################################

# Violin Plot

############################################################################################################



RB = c(p1_M500_high_RB, p1_M500_low_RB, p2_M500_high_RB, p2_M500_low_RB, p3_M500_high_RB, p3_M500_low_RB, p4_M500_high_RB, p4_M500_low_RB, p5_M200_high_RB, p5_M200_low_RB, p6_M200_high_RB, p6_M200_low_RB, p7_M200_high_RB, p7_M200_low_RB, p8_M200_high_RB, p8_M200_low_RB)
classes = factor(c(rep(1, length(c(p1_M500_high_RB, p1_M500_low_RB))), rep(2, length(c(p2_M500_high_RB, p2_M500_low_RB))), rep(3, length(c(p3_M500_high_RB, p3_M500_low_RB))), rep(4, length(c(p4_M500_high_RB, p4_M500_low_RB))), rep(5, length(c(p5_M200_high_RB, p5_M200_low_RB))), rep(6, length(c(p6_M200_high_RB, p6_M200_low_RB))), rep(7, length(c(p7_M200_high_RB, p7_M200_low_RB))), rep(8, length(c(p8_M200_high_RB, p8_M200_low_RB)))), levels = 1:8)
Level = factor(c( rep("High", length(p1_M500_high_RB)), rep("Low", length(p1_M500_low_RB)), rep("High", length(p2_M500_high_RB)), rep("Low", length(p2_M500_low_RB)), rep("High", length(p3_M500_high_RB)), rep("Low", length(p3_M500_low_RB)), rep("High", length(p4_M500_high_RB)), rep("Low", length(p4_M500_low_RB)),   rep("High", length(p5_M200_high_RB)), rep("Low", length(p5_M200_low_RB)), rep("High", length(p6_M200_high_RB)), rep("Low", length(p6_M200_low_RB)), rep("High", length(p7_M200_high_RB)), rep("Low", length(p7_M200_low_RB)), rep("High", length(p8_M200_high_RB)), rep("Low", length(p8_M200_low_RB))  ),levels  =c("High", "Low") )
length(RB) == length(classes) 
length(classes) == length(Level)

results = data.frame(RB=RB, classes = classes, Levels = Level)
head(results)

ggplot(results, aes(x = classes, y = RB, fill = Levels)) + geom_violin(trim = TRUE)   + stat_summary(fun =mean,geom="point", shape=16, size=2,  position = position_dodge2(width = 0.9,  preserve = "single")) + theme_classic() + ylab("Relative Bias") + xlab("Cases") +  theme(legend.position = "right", legend.title = element_text(size =14), legend.text = element_text(size =14),  axis.text = element_text(size=14), axis.title = element_text(size=14)) 


#=======================================================================================================================================
# eta
#======================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)
library(tidyverse)

# samples = beta_psi
# true_value = true_betas_psi
# PIP = gamma_psi

sumSim <- function(samples, true_value, PIP, nsim, indexes_covariates ){
  
  RB = CV= matrix(0, length(true_value), nsim)
  
  for (i in 1:nsim) {
    
    L <- length(indexes_covariates)
    index_present <- logical(L)
    
    for (l in 1:L) {
      if(any(indexes_covariates[l] == which(PIP[,1,i]>=0.5))){
        index_present[l] <- TRUE
      } else index_present[l] <- FALSE
      
    }
    
    beta_gamma <- samples[index_present,,i]
    true_beta_gamma = true_value[index_present]
    
    if(is.matrix(beta_gamma)){
      RB[index_present,i] =  (beta_gamma[,1] - true_beta_gamma)/true_beta_gamma
      CV[index_present,i] =   beta_gamma[,3]/ abs(beta_gamma[,1]) 
    } else {
      RB[index_present,i] =  (beta_gamma[1] - true_beta_gamma)/true_beta_gamma
      CV[index_present,i] =   beta_gamma[3]/ abs(beta_gamma[1]) 
    }
    
  }
  
  
  
  
  
  out = list( RB = RB, CV = CV)
  return(out)
  
}


nsim = 10
ncov = 10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/indexes_covariates.Rdata")


#######################################################################################################
# High
#######################################################################################################


# high case
#true_betas_psi <- c(1.8, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
#true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
#true_betas_p= true_betas_phi
 true_betas_eta <- c(1, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
#true_betas_p <- c(1.8, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept

# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_beta_eta.Rdata")


s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p1_M500_high_RB =  s$RB[which(s$RB!=0)]
p1_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_beta_eta.Rdata")


s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p2_M500_high_RB =  s$RB[which(s$RB!=0)]
p2_M500_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p3_M500_high_RB =  s$RB[which(s$RB!=0)]
p3_M500_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p4_M500_high_RB =  s$RB[which(s$RB!=0)]
p4_M500_high_CV =  s$CV[which(s$RB!=0)]

#========================================================================
# M= 200
#======================================================================


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p5_M200_high_RB =  s$RB[which(s$RB!=0)]
p5_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T10_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p6_M200_high_RB =  s$RB[which(s$RB!=0)]
p6_M200_high_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J5_T5_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p7_M200_high_RB =  s$RB[which(s$RB!=0)]
p7_M200_high_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M200_J2_T5_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p8_M200_high_RB =  s$RB[which(s$RB!=0)]
p8_M200_high_CV =  s$CV[which(s$RB!=0)]


#######################################################################################################
#Low
#######################################################################################################

#true_betas_psi <- c(-0.25, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
# true_betas_phi <- c(-1, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
 true_betas_eta <- c(-1.5, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
#true_betas_p <- c(-0.25, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept



# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p1_M500_low_RB =  s$RB[which(s$RB!=0)]
p1_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p2_M500_low_RB =  s$RB[which(s$RB!=0)]
p2_M500_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p3_M500_low_RB =  s$RB[which(s$RB!=0)]
p3_M500_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p4_M500_low_RB =  s$RB[which(s$RB!=0)]
p4_M500_low_CV =  s$CV[which(s$RB!=0)]




# M= 200


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T10_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p5_M200_low_RB =  s$RB[which(s$RB!=0)]
p5_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T10_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p6_M200_low_RB =  s$RB[which(s$RB!=0)]
p6_M200_low_CV =  s$CV[which(s$RB!=0)]


#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J5_T5_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p7_M200_low_RB =  s$RB[which(s$RB!=0)]
p7_M200_low_CV =  s$CV[which(s$RB!=0)]

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_gamma_eta.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M200_J2_T5_low_beta_eta.Rdata")

s = sumSim(beta_eta, true_betas_eta, gamma_eta, nsim, indexes_covariates)
p8_M200_low_RB =  s$RB[which(s$RB!=0)]
p8_M200_low_CV =  s$CV[which(s$RB!=0)]


############################################################################################################

# Violin Plot

############################################################################################################



RB = c(p1_M500_high_RB, p1_M500_low_RB, p2_M500_high_RB, p2_M500_low_RB, p3_M500_high_RB, p3_M500_low_RB, p4_M500_high_RB, p4_M500_low_RB, p5_M200_high_RB, p5_M200_low_RB, p6_M200_high_RB, p6_M200_low_RB, p7_M200_high_RB, p7_M200_low_RB, p8_M200_high_RB, p8_M200_low_RB)
classes = factor(c(rep(1, length(c(p1_M500_high_RB, p1_M500_low_RB))), rep(2, length(c(p2_M500_high_RB, p2_M500_low_RB))), rep(3, length(c(p3_M500_high_RB, p3_M500_low_RB))), rep(4, length(c(p4_M500_high_RB, p4_M500_low_RB))), rep(5, length(c(p5_M200_high_RB, p5_M200_low_RB))), rep(6, length(c(p6_M200_high_RB, p6_M200_low_RB))), rep(7, length(c(p7_M200_high_RB, p7_M200_low_RB))), rep(8, length(c(p8_M200_high_RB, p8_M200_low_RB)))), levels = 1:8)
Level = factor(c( rep("High", length(p1_M500_high_RB)), rep("Low", length(p1_M500_low_RB)), rep("High", length(p2_M500_high_RB)), rep("Low", length(p2_M500_low_RB)), rep("High", length(p3_M500_high_RB)), rep("Low", length(p3_M500_low_RB)), rep("High", length(p4_M500_high_RB)), rep("Low", length(p4_M500_low_RB)),   rep("High", length(p5_M200_high_RB)), rep("Low", length(p5_M200_low_RB)), rep("High", length(p6_M200_high_RB)), rep("Low", length(p6_M200_low_RB)), rep("High", length(p7_M200_high_RB)), rep("Low", length(p7_M200_low_RB)), rep("High", length(p8_M200_high_RB)), rep("Low", length(p8_M200_low_RB))  ),levels  =c("High", "Low") )
length(RB) == length(classes) 
length(classes) == length(Level)

results = data.frame(RB=RB, classes = classes, Levels = Level)
head(results)

ggplot(results, aes(x = classes, y = RB, fill = Levels)) + geom_violin(trim = TRUE)   + stat_summary(fun =mean,geom="point", shape=16, size=2,  position = position_dodge2(width = 0.9,  preserve = "single")) + theme_classic() + ylab("Relative Bias") + xlab("Cases") +  theme(legend.position = "right", legend.title = element_text(size =14), legend.text = element_text(size =14),  axis.text = element_text(size=14), axis.title = element_text(size=14)) 





# M= 500 PIP plots for all parameters

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

M= 500
J= 5;T=10

#==========================================================================================================================================
#p
#==========================================================================================================================================


# high case
true_betas_psi <- c(1.8, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
true_betas_phi <- c(1.3, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
true_betas_eta <- c(1, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
true_betas_p <- c(1.8, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept




# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_p.Rdata")
A = gamma_p

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_p.Rdata")
B = gamma_p

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_p.Rdata")
C = gamma_p

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_p.Rdata")
D = gamma_p


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

p_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(p_PIP)

# Change line types and point shapes
ggplot(p_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))

#=============================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

M= 500
J= 5;T=10

# Low case
true_betas_psi <- c(-0.25, 0.4, 0, 0.5, 0, -0.6, -0.7,0,0,0, 0.25,0.65,0,0, 0) # includes intercept
true_betas_phi <- c(-1, -0.4, 0.87, 0.5, 0, -0.6, 0.7,-0.75,0.1,0.2, 0,0,0,-0.33, 0.22) # includes intercept 
true_betas_eta <- c(-1.5, -0.7, 0, 0.65, -0.33, -0.2, 0.8,0,0,0, 0.71,-0.25,1,0, 0) #  includes intercept
true_betas_p <- c(-0.25, 0, 0.25, 0.75, 0, -0.1, -0.2,0.84,-0.05,0.12, 0,0,0,-0.75, 0.12) #  includes intercept

# M= 500

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_p.Rdata")
A = gamma_p

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_p.Rdata")
B = gamma_p

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_p.Rdata")
C = gamma_p

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_p.Rdata")
D = gamma_p


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

p_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(p_PIP)

# Change line types and point shapes
ggplot(p_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))


#==========================================================================================================================================
#phi
#==========================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

# M= 500 high case


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_phi.Rdata")
A = gamma_phi

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_phi.Rdata")
B = gamma_phi

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_phi.Rdata")
C = gamma_phi

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_phi.Rdata")
D = gamma_phi


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

phi_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(phi_PIP)

# Change line types and point shapes
ggplot(phi_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))

#=============================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

# M= 500 low 

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_phi.Rdata")
A = gamma_phi

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_phi.Rdata")
B = gamma_phi

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_phi.Rdata")
C = gamma_phi

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_phi.Rdata")
D = gamma_phi


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

p_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(p_PIP)

# Change line types and point shapes
ggplot(p_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))


#==========================================================================================================================================
#zeta(eta)
#==========================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

# M= 500 high case


#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_eta.Rdata")
A = gamma_eta

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_eta.Rdata")
B = gamma_eta

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_eta.Rdata")
C = gamma_eta

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_eta.Rdata")
D = gamma_eta


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

phi_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(phi_PIP)

# Change line types and point shapes
ggplot(phi_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))

#=============================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

# M= 500 low 

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_eta.Rdata")
A = gamma_eta

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_eta.Rdata")
B = gamma_eta

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_eta.Rdata")
C = gamma_eta

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_eta.Rdata")
D = gamma_eta


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

p_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(p_PIP)

# Change line types and point shapes
ggplot(p_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))

#==========================================================================================================================================
#psi
#==========================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

# M= 500 high case

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_1_gamma_psi.Rdata")
A = gamma_psi

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T10_gamma_psi.Rdata")
B = gamma_psi

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J5_T5_gamma_psi.Rdata")
C = gamma_psi

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/High_case/BVS_DOM_M500_J2_T5_gamma_psi.Rdata")
D = gamma_psi


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

phi_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(phi_PIP)

# Change line types and point shapes
ggplot(phi_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))

#=============================================================================================================================================

rm(list=ls())
library(tidyverse)
library(ggmcmc)
library(coda)

# M= 500 low 

#J=5;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T10_low_gamma_psi.Rdata")
A = gamma_psi

#J=2;T=10
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T10_low_gamma_psi.Rdata")
B = gamma_psi

#J=5;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J5_T5_low_gamma_psi.Rdata")
C = gamma_psi

#J=2;T=5
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/BVS Dyanmic Occupancy/BVS DOM/Simulation results/Low case/BVS_DOM_M500_J2_T5_low_gamma_psi.Rdata")
D = gamma_psi


mean(A[1,1,])
mean(A[2,1,])
APIP = apply(A[-1,1,], 1, mean)
APIP
BPIP = apply(B[-1,1,], 1, mean)
CPIP = apply(C[-1,1,], 1, mean)
DPIP = apply(D[-1,1,], 1, mean)

variables = as.factor(rep(1:10, 4)) 
classes = factor(rep(c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"), each = 10), levels = c("J=5,T=10", "J=2,T=10", "J=5,T=5", "J=2,T=5"))
#classes = as.factor(rep(c("A", "B", "C", "D"), each = 10))
class(classes)
PIP = c(APIP, BPIP,CPIP, DPIP)

p_PIP = data.frame(Classes = classes, Variable = variables, PIP= PIP)
head(p_PIP)

# Change line types and point shapes
ggplot(p_PIP, aes(x=Variable, y=PIP, group=Classes)) +geom_line(aes(linetype=Classes))+geom_point(aes(shape=Classes), size = 2.5) + theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +geom_hline(yintercept=0.5,linetype=8) #+ scale_group_discrete(labels=c('1', '2', '3','4'))

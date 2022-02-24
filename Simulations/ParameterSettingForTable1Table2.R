rm(list=ls())

########################################################
## Generate parameter setting for Table 1             ##
## for simulation w.r.t null and alternative          ##
## modalilty: A & B, Truth: normal & disease          ##
## Variance components the same for normal & disease  ##
## i.e. b=1                                           ##
## AUC_B  = c(0.6,0.8,0.9)                            ##
## AUC_A = AUC_B (tau_DA = 0) for null simulation     ## 
## AUC_A = AUC_B + dAUC for alternative simulation    ## 
## dAUC = 0.03 & 0.05 & 0.051 & 0.06                  ##
## mu_N = 0; tau_NB = tau_DB = 0; tau_NA = 0          ##
## need to calculate: mu_D and tau_DA                 ##
########################################################

# the structure of combinations of varC (H or L) and varR (H or L) #
#varStruc = rep(c("HL","LL","HH","LH"), each =3)
varStruc = c("HL1","HL2","HL3","LL1","LL2","LL3","HH1","HH2","HH3","LH1","LH2","LH3")

# Low var_R and var_tau.R for Table 1 #
varRL1 = matrix(c(0.0055,0.0055,0.0055),ncol=1); varRL1 = cbind(varRL1,varRL1);
# High var_R and var_tau.R for Table 1 #
varRH1 = matrix(c(0.011,0.030,0.056),ncol=1); varRH1 = cbind(varRH1,varRH1); 
## var_R and var_tau.R for Table 1 ##
varR1 = rbind(varRL1,varRL1, varRH1,varRH1)

# High var_C, var_tau.C, var_R.C and var_epsilon for Table 1 (the same for normal and disease) #
varCH1 = matrix(c(0.3,0.3,0.2,0.2),ncol=4); varCH1 = rbind(varCH1,varCH1,varCH1); 
# Low var_C, var_tau.C, var_R.C and var_epsilon for Table 1 (the same for normal and disease) #
varCL1 = matrix(c(0.1,0.1,0.2,0.6),ncol=4); varCL1 = rbind(varCL1,varCL1,varCL1); 
## var_C, var_tau.C, var_R.C and var_epsilon for Table 1 (the same for normal and disease) ##
varC1 = rbind(varCH1,varCL1,varCH1,varCL1)

## var_C, var_tau.C, var_R.C and var_epsilon & var_R and var_tau.R for Table 1 ##
var1 = cbind(varC1,varR1)
colnames(var1) = c("var_C-","var_tau.C-", "var_R.C-", "var_eps-", "var_R", "var_tau.R")


# AUC for modality B #
AUCB = rep(c(0.6,0.8,0.9),4) 
# mu_N = 0; tau_NA = tau_DA = 0; tauNB = 0 #
muN = tauNB = tauDB = tauNA = 0 
#### calculate muD for Table 1 w.r.t. alternative simulation ####
muD1 = sqrt(2*apply(var1,1,sum))*qnorm(AUCB) - tauDB

#### the parameter setting for Table 1 w.r.t alternative simulation for different dAUC ####
# differneces of AUCA and AUCB #
dAUC = c(0.03,0.05,0.051,0.06)
for(i in 1:length(dAUC)){
  # AUC for modality A #
  AUCA = AUCB + dAUC[i]
  #### calculate tauDAS for Table 1 w.r.t. alternative simulation ####
  tauDA1 = sqrt(2*apply(var1,1,sum))*qnorm(AUCA) - muD1
  #### the parameter setting for Table 1 w.r.t alternative simulation for different dAUC ####
  temp = data.frame(varStruc=varStruc,mu_1 = muD1, tau_1A = tauDA1,
                     AUC_B = AUCB, b=rep(1,length(varStruc)),var = var1)
  colnames(temp) = c("varStruc ","mu1","tauA1","AUC_B","b", "Vc0","Vtauc0", "Vrc0", "Veps0","Vr0", "Vtaur0")
  assign(paste("para.alter",as.character(dAUC[i]*1000),sep=""),temp)
}

#### the parameter setting for Table 1 w.r.t null simulation ####
para.null = eval(as.name(paste("para.alter",as.character(dAUC[1]*1000),sep="")))
para.null$tauA1 = 0

## save the  parameter setting for Table 1 ##
setwd("/home/zhuang/Desktop/MRMC/CodesAdaptive/FullyCross/Perceus/ParameterSetting")
fname = paste(getwd(),sprintf("/ParameterSettingForTable1.Rdata"),sep="")
#save(para.null,para.alter3,para.alter5, file =fname)
paraSettings= grep("^[para.]",ls(),value = TRUE); 
save(list = paraSettings,file=fname); rm(paraSettings)
#load(file=fname)



##############################################################################
## Generate parameter setting for Table 2                                   ##
## for simulation w.r.t null and alternative                                ##
## modalilty: A & B, Truth: normal & disease                                ##
## Variance components var_R and var_tau.R the same for normal & disease    ##
## diseae (var_C,var_tau.C, var_R.C, var_eps) equals b^-2 times normal ones ##
## where b = c(0.84566,0.71082,0.55140)                                     ##
## AUC_B  = c(0.6,0.8,0.9)                                                  ##
## AUC_A = AUC_B (tau_DA = 0) for null simulation                           ##
## AUC_A = AUC_B + dAUC for alternative simulation                          ## 
## dAUC = 0.03 & 0.05 & 0.051 & 0.06                                        ##
## mu_N = 0; tau_NB = tau_DB = 0; tau_NA = 0                                ##
## need to calculate: mu_D and tau_DA                                       ##
##############################################################################

# Low var_R and var_tau.R for Table 2 #
varRL2 = matrix(c(0.0066,0.0082,0.0118),ncol=1); varRL2 = cbind(varRL2,varRL2);
# High var_R and var_tau.R for Table 2 #
varRH2 = matrix(c(0.0132,0.0447,0.1201),ncol=1); varRH2 = cbind(varRH2,varRH2);
## var_R and var_tau.R for Table 2 ##
varR2 = rbind(varRL2,varRL2, varRH2,varRH2)

## var_C, var_tau.C, var_R.C and var_epsilon for Table 2 with normal cases ##
## the same as Table 1
varCN2 = varC1

## var_C, var_tau.C, var_R.C and var_epsilon for Table 2 with disease cases ##
## b = c(0.84566,0.71082,0.55140)
b = c(0.84566,0.71082,0.55140); fullb = rep(b,4)
varCD2 = varCN2*fullb^(-2);
#varCD2 = round(varCD2,2)

## var_C-, var_tau.C-, var_R.C- and var_epsilon- &
## var_C+, var_tau.C+, var_R.C+ and var_epsilon+ &
## var_R and var_tau.R for Table 2 ##
var2 = cbind(varCN2,varCD2, varR2)
colnames(var2) = c("var_C-","var_tau.C-", "var_R.C-", "var_eps-",
                   "var_C+","var_tau.C+", "var_R.C+", "var_eps+","var_R", "var_tau.R")


# AUC for modality B (The same as for Table 1) #
AUCB = rep(c(0.6,0.8,0.9),4)
# mu_N = 0; tau_NB = tau_DB = 0; tauNA = 0 (The same as for Table 1) #
muN = tauNB = tauDB = tauNA = 0

#### calculate muD and tauDA for Table 2 w.r.t alternative simulation ####
muD2 = sqrt(apply(var2,1,sum)+apply(varR2,1,sum))*qnorm(AUCB) - tauDB


#### the parameter setting for Table 2 w.r.t alternative simulation for different dAUC ####
# differneces of AUCA and AUCB #
for(i in 1:length(dAUC)){
  # AUC for modality A #
  AUCA = AUCB + dAUC[i]
  #### calculate tauDAS for Table 1 w.r.t. alternative simulation ####
  tauDA2 = sqrt(apply(var2,1,sum)+apply(varR2,1,sum))*qnorm(AUCA) - muD2
  #### the parameter setting for Table 2 w.r.t alternative simulation for different dAUC ####
  temp = data.frame(varStruc=varStruc,mu_1 = muD2, tau_1A = tauDA2,
                    AUC_B = AUCB, b = fullb, var = cbind(varCN2,varR2))
  colnames(temp) = c("varStruc","mu1","tauA1","AUC_B","b","Vc0","Vtauc0", "Vrc0", "Veps0","Vr0", "Vtaur0")
  assign(paste("para.alter",as.character(dAUC[i]*1000),sep=""),temp)
}

#### the parameter setting for Table 1 w.r.t null simulation ####
para.null = eval(as.name(paste("para.alter",as.character(dAUC[1]*1000),sep="")))
para.null$tauA1 = 0


## save the  parameter setting for Table 2 ##
setwd("/home/zhuang/Desktop/MRMC/CodesAdaptive/FullyCross/Perceus/ParameterSetting")
fname = paste(getwd(),sprintf("/ParameterSettingForTable2.Rdata"),sep="")
#save(para.null,para.alter3,para.alter5,file =fname)
paraSettings= grep("^[para.]",ls(),value = TRUE); 
save(list = paraSettings,file=fname); rm(paraSettings)
#load(file=fname)
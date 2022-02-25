#######################################################################################
############################# Setting when running on Cluster #########################
## pass the random seed from Rscript command line
args = commandArgs(trailingOnly = TRUE)
if(length(args)< 5){
  stop("At least five arguments must be suppied:
       1. adaptive method: Adaptive1 or Adaptive2
       2. null or alternative simulation: null or alter50
       3. inverse of ratio of diseased cases to normal cases: 1 or 2
       4. underpower/overpower for alter simulations or NA for null simulation with initial sample size: under or over or NA
       5. random seed: number 1-1000", call. = FALSE)
}
AdapMethod = args[1]; Hypo = args[2]; Ratio.Nc0toNc1 = as.numeric(args[3]); Power0 = args[4]; SetSeed = as.numeric(args[5])

library(AdaptiveMRMC)
## define function for pasting directory
"%+%" = function(x,y) paste(x,y,sep="")

## set the directory to  source R functions, input files and output files
indir = "/home/zhuang/MRMC/Rpackage/Testing"
workdir = "/raidb/zhuang/" %+% AdapMethod
outdir = ifelse(Hypo == "null", workdir %+% "/result." %+% Hypo %+% "." %+% Ratio.Nc0toNc1,
                workdir %+% "/result." %+% Hypo %+% "." %+% Ratio.Nc0toNc1 %+% "." %+% Power0)

## input the parameters setted in data generation
## ParaSet = c(mu1,tauB1,b,Vc0,Vtauc0,Vrc0,Veps0,Vr0,Vtaur0)
Para = "para."; Tbl = "Table1";  # or: Tbl = "Table2";
load(indir %+% "/ParameterSetting/ParameterSettingFor" %+% Tbl %+% ".Rdata")
ParaSet = eval(as.name(Para %+% Hypo))

## the output file name
fname <- outdir %+% sprintf("/results_seed_%05d", SetSeed) %+% ".Rdata"

setwd(workdir) ## set the work directory

#######################################################################################
########################## Setting when running on local computer #####################
# rm(list=ls())
# library(AdaptiveMRMC)
#
# ## define function for pasting directory
# "%+%" = function(x,y) paste(x,y,sep="")
#
# ## decide which adaptive method to evaulate (1. resize reader only; 2. resize both reader and case)
# AdapMethod = "Adaptive1"; AdapMethod = "Adaptive2"
#
# ## More settings
# Hypo = "null";  Hypo = "alter50";
# Ratio.Nc0toNc1 = 1; Ratio.Nc0toNc1 = 2
# Power0 = "under"; #Power0 = "over"; #Power0 = NA
# SetSeed = 2
#
# ## set the directories of files input, running simulations and files output
# indir = "/home/zhuang/Desktop/MRMC/Rpackage/Testing"
# #indir = "C:/Users/Zhipeng.Huang/Desktop/MRMC/Rpackage/Testing"
# workdir = indir %+% "/" %+% AdapMethod
# outdir = ifelse(Hypo == "null", workdir %+% "/result." %+% Hypo %+% "." %+% Ratio.Nc0toNc1,
#                 workdir %+% "/result." %+% Hypo %+% "." %+% Ratio.Nc0toNc1 %+% "." %+% Power0)
#
# ## input the parameters setted in data generation
# ## ParaSet = c(mu1,tauB1,b,Vc0,Vtauc0,Vrc0,Veps0,Vr0,Vtaur0)
# Para = "para."; Tbl = "Table1";  # or: Tbl = "Table2";
# load(indir %+% "/ParameterSetting/ParameterSettingFor" %+% Tbl %+% ".Rdata")
# ParaSet = eval(as.name(Para %+% Hypo))
#
# ## the output file name
# fname <- outdir %+% sprintf("/results_seed_%05d", SetSeed) %+% ".Rdata"
#
# setwd(workdir) ## set the work directory

#######################################################################################
################################### More Settings #####################################
paraID = c(2,5,9,11)
ParaSet = ParaSet[paraID,]

SampleSet = numeric(0); names.SampleSet = c("Nc.0","Nr.init","Nr.1","Nr.max","Nr.expect","Nc.0.max")
if(Hypo == "alter50" && Ratio.Nc0toNc1 == 1 && Power0 == "under"){
  SampleSet = matrix(c(c(120,12,6,40,999,300),
                       c(100,5,3,30,999,220),
                       c(120,8,4,30,999,250),
                       c(120,12,6,40,999,250)),byrow = TRUE,ncol = length(names.SampleSet))
}else if(Hypo == "alter50" && Ratio.Nc0toNc1 == 1 && Power0 == "over"){
  SampleSet = matrix(c(c(250,30,4,40,7,300),
                       c(150,20,3,30,6,220),
                       c(200,25,5,40,11,300),
                       c(220,30,6,40,17,300)),byrow = TRUE,ncol = length(names.SampleSet))
}else if(Hypo == "alter50" && Ratio.Nc0toNc1 == 2 && Power0 == "under"){
  SampleSet = matrix(c(c(180,12,6,40,999,400),
                       c(150,5,3,30,999,300),
                       c(180,8,4,30,999,330),
                       c(180,12,6,40,999,330)),byrow = TRUE,ncol = length(names.SampleSet))
}else if(Hypo == "alter50" && Ratio.Nc0toNc1 == 2 && Power0 == "over"){
  SampleSet = matrix(c(c(340,30,4,40,7,400),
                       c(220,20,3,30,6,300),
                       c(280,25,5,40,11,400),
                       c(320,30,6,40,17,400)),byrow = TRUE,ncol = length(names.SampleSet))
}else if(Hypo == "null" && Ratio.Nc0toNc1 == 1){
  SampleSet = matrix(rep(c(200,20,10,40,999,300),4),byrow = TRUE,ncol = length(names.SampleSet))
}else if(Hypo == "null" && Ratio.Nc0toNc1 == 2){
  SampleSet = matrix(rep(c(240,20,10,40,999,320),4),byrow = TRUE,ncol = length(names.SampleSet))
}

## add one more column if adaptive method II is adopted
if(AdapMethod == "Adaptive2"){
  SampleSet = cbind(SampleSet[,1], SampleSet)
  names.SampleSet = c("Nc.10", "Nc.20", names.SampleSet[-1])
}
SampleSet = as.data.frame(SampleSet)
dimnames(SampleSet) = list(paraID, names.SampleSet)
SampleSet = data.frame(SampleSet,R.Nc0toNc1 = Ratio.Nc0toNc1)

## array to instore the simulation results
n.sim = 10; var.struc = ParaSet$varStruc
names.outcome = c("Nr.adap", "cv.adap", paste("M", 1:8, sep = ""), "Vr.inter", "Vc.inter", "Vcr.Ratio.inter",
                  "dAUC.inter", "var.dAUC.inter", "z.inter", "CE.inter", "CP.inter", "CP.adap.low", "CP.adap.up", "CP.adap", "s0", "s1",
                  "dAUC.init", "var.dAUC.init", "z.init", "dAUC.adap", "var.dAUC.adap", "z.adap",
                  "reject.inter.z", "reject.inter.t", "reject.init.z", "reject.init.t",
                  "reject.adap.z", "reject.adap.t", "reject.adap.cv")
if(AdapMethod == "Adaptive2"){
  names.outcome = c("Nr.adap", "Nc.20.adap", "cv.adap", paste("M", 1:8, sep = ""), "Vr.1", "Vc.1", "Vr.2", "Vc.2", "Vr.2.adap", "Vc.2.adap",
                    "dAUC.inter", "sd.inter", "sd.cond2", "sd.total", "z.inter", "CE.inter", "CP.inter", "CC1", "CC2",
                    "CP.adap.low", "CP.adap.up", "CP.adap", "N.sim.na", "sd.cond2.adap", "sd.total.adap",
                    "dAUC.init", "sd.init", "z.init", "dAUC.adap", "sd.adap", "z.adap", "varUstat.adap.na", "varMLE.adap.na",
                    "reject.inter.z", "reject.inter.t", "reject.init.z", "reject.init.t",
                    "reject.adap.z", "reject.adap.t", "reject.adap.cv")
}
sim.result = array(0,c(n.sim, length(names.outcome), length(var.struc)), dimnames = list(NULL, names.outcome, var.struc))

########### run the simulation ###########
runfun = ifelse(AdapMethod == "Adaptive2", sim.adaptiveII, sim.adaptiveI)
time.sim = proc.time()
for(kk in 1:dim(SampleSet)[1]){
  for(ss in 1:n.sim) sim.result[ss,,kk] <- try(runfun(SampleSet[kk,], ParaSet[kk,], 100000*SetSeed+1000*(ss-1)+10*kk),TRUE)
}
time.sim = proc.time() - time.sim

## save the raw result to the output file
save(sim.result, time.sim, file = fname)
#load(fname)

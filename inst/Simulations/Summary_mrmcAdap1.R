rm(list=ls())

time.summ = proc.time()
## part of the simulation setting
AdapMethod = "Adaptive1"
Hypo = "alter50"; #Hypo = "null"
Ratio.Nc0toNc1 = 1; #Ratio.Nc0toNc1 = 2
Power0 = "under"; #Power0 = "over"; #Power0 = NA
alpha0 = 0.025; z.alpha = qnorm(1-alpha0); beta0 = 0.2; #z.beta = qnorm(1-beta0)

if(Hypo != "null") library(AdaptiveMRMC)

## set the directories of parameters input
indir = "/home/zhuang/Desktop/MRMC/Rpackage/Testing"
#indir = "C:/Users/Zhipeng.Huang/Desktop/MRMC/Rpackage/Testing"

## the parameter settings for the simulations
Para = "para."; Tbl = "Table1";  #Tbl = "Table2";
load(paste(indir, "/ParameterSetting/ParameterSettingFor", Tbl, ".Rdata",sep=""))
ParaSet = eval(as.name(paste(Para, Hypo,sep="")))
paraID = c(2,5,9,11)
ParaSet = ParaSet[paraID,]

## the sample settings for the simulations
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
## transfer to data.frame and name it
SampleSet = as.data.frame(SampleSet)
dimnames(SampleSet) = list(paraID, names.SampleSet)
SampleSet = data.frame(SampleSet,R.Nc0toNc1 = Ratio.Nc0toNc1)

#######################################################################################
## calculate the theoratical power under alternative
if(Hypo != "null"){
  setwd(indir); reject.init.theory = reject.min.theory = reject.max.theory = numeric(0)
  for(kk in 1:dim(ParaSet)[1]){
    nr = SampleSet[kk,"Nr.init"]; nr.min = SampleSet[kk,"Nr.1"] + 1; nr.max = SampleSet[kk,"Nr.max"];
    nc0 = SampleSet[kk,"Nc.0"]; nc1 = round(nc0/Ratio.Nc0toNc1);
    tempt = analytical.var.AUC(nr, nc0, nc1, ParaSet[kk,])
    tempt.min = analytical.var.AUC(nr.min, nc0, nc1, ParaSet[kk,])
    tempt.max = analytical.var.AUC(nr.max, nc0, nc1, ParaSet[kk,])
    tempt$pow; tempt.max$pow; tempt.min$pow
    reject.init.theory = c(reject.init.theory,tempt$pow)
    reject.min.theory = c(reject.min.theory,tempt.min$pow)
    reject.max.theory = c(reject.max.theory,tempt.max$pow)
  }
}

#######################################################################################
## set the directories of raw simulation result input and result summary output
rawdir = ifelse(Hypo == "null", paste(indir, "/", AdapMethod, "/result.", Hypo, ".", Ratio.Nc0toNc1, sep=""),
                paste(indir, "/", AdapMethod, "/result.", Hypo, ".", Ratio.Nc0toNc1, ".", Power0, sep=""))
sumdir = paste(indir, "/Results", sep="")
setwd(rawdir)

#######################################################################################
names.outcome = c("Nr.adap", "cv.adap", paste("M", 1:8, sep = ""), "Vr.inter", "Vc.inter", "Vcr.Ratio.inter",
                  "dAUC.inter", "var.dAUC.inter", "z.inter", "CE.inter", "CP.inter", "CP.adap.low", "CP.adap.up", "CP.adap", "s0", "s1",
                  "dAUC.init", "var.dAUC.init", "z.init", "dAUC.adap", "var.dAUC.adap", "z.adap",
                  "reject.inter.z", "reject.inter.t", "reject.init.z", "reject.init.t",
                  "reject.adap.z", "reject.adap.t", "reject.adap.cv")

## number of parallel computations in cluster, simulation replicates in each cluster node and variance structures
n.paral = 1000; n.sim = 10; n.outcome = length(names.outcome)
names.varstr = ParaSet$varStruc; n.varstr = length(names.varstr)

###################################################
## save results one by one from 1:n.paral and 1:n.sim
fname.v = Sys.glob("result*"); n.paral = length(fname.v); n.total = n.paral*n.sim
comb.sim.result = comb.sim.result0 = array(0,c(n.total,n.outcome,n.varstr),dimnames = list(NULL,names.outcome,names.varstr))
seed.v = id.na.mat.1 = numeric(0)
for(setseed in 1:n.paral){
  fname = fname.v[setseed]; load(fname)
  id.seed = as.numeric(sub(".Rdata", "", sub(".*_", "", fname)))
  seed.v = c(seed.v,id.seed)
  comb.sim.result0[(n.sim*(setseed-1)+1):(n.sim*setseed),,] = sim.result

  for(i in 1:n.sim){
    tempt = apply(sim.result[i,,],2,as.numeric)
    if(sum(is.na(tempt))>0) id.na.mat.1 = rbind(id.na.mat.1,c(id.seed,i,which(is.na(tempt),arr.ind = T)[1,2]))
  }
}
dimnames(id.na.mat.1) = list(NULL, c("setseed","ss","kk"))
for(i in 1:n.total) comb.sim.result[i,,] = apply(comb.sim.result0[i,,],2,as.numeric)
id.na.mat.2 = which(is.na(comb.sim.result),arr.ind = T); dimnames(id.na.mat.2)[[2]][3]="kk"
val.na.v = unique(comb.sim.result0[id.na.mat.2]); val.na.v
id.NaN = which(comb.sim.result0=="NaN",arr.ind = T); id.NaN
#seed.v.fail = setdiff(1:1000,seed.v)

## plot critical value v.s. conditional power
id.check = 2
v.check = comb.sim.result[,,id.check];
id.na.check = unique(which(is.na(v.check),arr.ind = T)[,1]); id.na.check
id.vc.ne.check = which(v.check[,"Vc.inter"]<=0); id.vc.ne.check

plot(v.check[,"CP.inter"],v.check[,"cv.adap"])
abline(h=z.alpha,col="red")
abline(v=1-beta0,col="red")

if(length(c(id.na.check,id.vc.ne.check))==0){
  plot(v.check[,"CP.inter"],v.check[,"cv.adap"])
}else{
  plot(v.check[-c(id.na.check,id.vc.ne.check),"CP.inter"],v.check[-c(id.na.check,id.vc.ne.check),"cv.adap"])
}
abline(h=z.alpha,col="red")
abline(v=1-beta0,col="red")

plot(v.check[,"z.inter"],v.check[,"cv.adap"])
plot(v.check[,"z.inter"],v.check[,"CP.inter"])
plot(v.check[,"z.inter"],v.check[,"Nr.adap"])
plot(v.check[,"CP.inter"],v.check[,"Nr.adap"])


## the mean measurments over the simulations
result.mean = n.effect = n.cplowG = n.cplowG.interRej = n.cplowG.adapNoRej = n.cpupL = n.cpupL.adapRej = n.neg = n.neg.initRej = numeric(0)
for(kk in 1:n.varstr){
  id.na = unique(which(is.na(comb.sim.result[,,kk]),arr.ind = T)[,1])
  id.vc.ne = which(comb.sim.result[,"Vc.inter",kk]<=0)
  id.effect = setdiff(1:(n.total), c(id.na,id.vc.ne)); n.effect = c(n.effect, length(id.effect))
  result.mean = cbind(result.mean,apply(comb.sim.result[id.effect,,kk],2,mean))
  
  n.cplowG = c(n.cplowG, sum(comb.sim.result[id.effect,"CP.adap.low",kk]>1-beta0&comb.sim.result[id.effect,"CP.adap.up",kk]>1-beta0))
  n.cplowG.interRej = c(n.cplowG.interRej, sum(comb.sim.result[id.effect,"CP.adap.low",kk]>1-beta0
                                                 &comb.sim.result[id.effect,"CP.adap.up",kk]>1-beta0
                                                 &comb.sim.result[id.effect,"reject.inter.z",kk]==1))
  n.cplowG.adapNoRej = c(n.cplowG.adapNoRej, sum(comb.sim.result[id.effect,"CP.adap.low",kk]>1-beta0
                                          &comb.sim.result[id.effect,"CP.adap.up",kk]>1-beta0
                                          &comb.sim.result[id.effect,"reject.adap.cv",kk]==0))
  n.cpupL = c(n.cpupL, sum(comb.sim.result[id.effect,"CP.adap.low",kk]<1-beta0&comb.sim.result[id.effect,"CP.adap.up",kk]<1-beta0))
  n.cpupL.adapRej = c(n.cpupL.adapRej, sum(comb.sim.result[id.effect,"CP.adap.low",kk]<1-beta0
                                    &comb.sim.result[id.effect,"CP.adap.up",kk]<1-beta0
                                    &comb.sim.result[id.effect,"reject.adap.cv",kk]==1))
  n.neg = c(n.neg, sum(comb.sim.result[id.effect,"z.inter",kk]<=0))
  n.neg.initRej = c(n.neg.initRej, sum(comb.sim.result[id.effect,"z.inter",kk]<=0&comb.sim.result[id.effect,"reject.adap.cv",kk]==1))
}
if(Hypo == "null"){
  result.mean = rbind(result.mean, n.effect, 
                      n.cplowG, n.cplowG.interRej, n.cplowG.adapNoRej, n.cpupL, n.cpupL.adapRej, n.neg, n.neg.initRej)
}else{
  result.mean = rbind(result.mean, reject.init.theory, reject.min.theory, reject.max.theory, n.effect, 
                      n.cplowG, n.cplowG.interRej, n.cplowG.adapNoRej, n.cpupL, n.cpupL.adapRej, n.neg, n.neg.initRej)
}
colnames(result.mean) = names.varstr

n.digit = 4;
## mean error/power difference for the 12 variance structures
if(Hypo == "null"){
  Mean.Error.Diff = result.mean["reject.adap.z",]-result.mean["reject.adap.cv",]
  names(Mean.Error.Diff) = names.varstr
  print(round(Mean.Error.Diff, n.digit))
}else{
  Mean.Power.Diff = result.mean["reject.adap.cv",]-result.mean["reject.init.z",]
  names(Mean.Power.Diff) = names.varstr
  print(round(Mean.Power.Diff, n.digit))
}

## show part of the results
if(Hypo == "null"){
  print(round(result.mean[c("reject.inter.z", "reject.inter.t", "reject.init.z", "reject.init.t",
                            "reject.adap.z", "reject.adap.t", "reject.adap.cv", "Nr.adap", "n.effect",
                            "n.cplowG", "n.cplowG.interRej", "n.cplowG.adapNoRej", "n.cpupL", "n.cpupL.adapRej", "n.neg", "n.neg.initRej"),], n.digit))
}else{
  print(round(rbind(result.mean[c("reject.inter.z", "reject.inter.t", "reject.init.z", "reject.init.t",
                                  "reject.init.theory", "reject.min.theory", "reject.max.theory", "reject.adap.cv"),],
                    t(SampleSet[,c("Nr.1","Nr.init","Nr.expect")]), result.mean[c("Nr.adap","n.effect", 
                                                                                  "n.cplowG", "n.cplowG.interRej", "n.cplowG.adapNoRej",
                                                                                  "n.cpupL", "n.cpupL.adapRej",
                                                                                  "n.neg", "n.neg.initRej"),]), n.digit))
}
time.summ = proc.time() - time.summ

#######################################################################################
## save the simulation summary
fname = ifelse(Hypo == "null", paste("summary_", AdapMethod, "_", Hypo, "_", Ratio.Nc0toNc1, ".Rdata", sep=""),
               paste("summary_", AdapMethod, "_", Hypo, "_", Ratio.Nc0toNc1, "_", Power0, ".Rdata", sep=""))
rm(setseed, kk, sim.result, time.sim, comb.sim.result0, comb.sim.result, id.check, v.check)
save.image(file = paste(sumdir,"/", fname,sep=""))
#load(fname)

#######################################################################################



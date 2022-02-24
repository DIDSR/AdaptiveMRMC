rm(list=ls())
library(AdaptiveMRMC)

time.summ = proc.time()
## part of the simulation setting
AdapMethod = "Adaptive2"
Hypo = "alter50"; #Hypo = "null"
Ratio.Nc0toNc1 = 1; #Ratio.Nc0toNc1 = 2
Power0 = "under"; #Power0 = "over"; #Power0 = NA
alpha0 = 0.025; z.alpha = qnorm(1-alpha0); beta0 = 0.2; #z.beta = qnorm(1-beta0)

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

SampleSet = as.data.frame(SampleSet)
dimnames(SampleSet) = list(paraID, names.SampleSet)
SampleSet = data.frame(SampleSet,R.Nc0toNc1 = Ratio.Nc0toNc1)

#######################################################################################
## calculate the theoratical power under alternative
if(Hypo != "null"){
  setwd(indir); reject.init.theory = reject.min.theory = reject.max.theory = numeric(0)
  for(kk in 1:dim(ParaSet)[1]){
    nr = SampleSet[kk,"Nr.init"]; nr.min = SampleSet[kk,"Nr.1"] + 1; nr.max = SampleSet[kk,"Nr.max"]
    nc0 = SampleSet[kk,"Nc.10"]; nc1 = round(nc0/Ratio.Nc0toNc1);
    nc0.max = SampleSet[kk,"Nc.0.max"]; nc1.max = round(nc0.max/Ratio.Nc0toNc1)
    tempt = analytical.var.AUC(nr, nc0, nc1, ParaSet[kk,])
    tempt.min = analytical.var.AUC(nr.min, nc0, nc1, ParaSet[kk,])
    tempt.max = analytical.var.AUC(nr.max, nc0.max, nc1.max, ParaSet[kk,])
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
names.outcome = c("Nr.adap", "Nc.20.adap", "cv.adap", paste("M", 1:8, sep = ""), "Vr.1", "Vc.1", "Vr.2", "Vc.2", "Vr.2.adap", "Vc.2.adap",
                  "dAUC.inter", "sd.inter", "sd.cond2", "sd.total", "z.inter", "CE.inter", "CP.inter", "CC1", "CC2",
                  "CP.adap.low", "CP.adap.up", "CP.adap", "N.sim.na", "sd.cond2.adap", "sd.total.adap",
                  "dAUC.init", "sd.init", "z.init", "dAUC.adap", "sd.adap", "z.adap", "varUstat.adap.na", "varMLE.adap.na",
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
id.NaN = which(comb.sim.result0=="NaN",arr.ind = T); id.NaN; comb.sim.result0[id.NaN]
#seed.v.fail = setdiff(1:1000,seed.v)

#################################################################################
## check the simulations with NA and the simulations with negative Vr.2.adap
ana = apply(comb.sim.result[,"N.sim.na",]!=0,2,sum)
anan = apply(comb.sim.result[,"N.sim.na",]!=0&comb.sim.result[,"z.inter",]<=0,2,sum)
anap = apply(comb.sim.result[,"N.sim.na",]!=0&comb.sim.result[,"z.inter",]>0,2,sum)
anapg = apply(comb.sim.result[,"N.sim.na",]!=0&comb.sim.result[,"z.inter",]>0&comb.sim.result[,"CP.adap.low",]>1-beta0,2,sum)
anapl = apply(comb.sim.result[,"N.sim.na",]!=0&comb.sim.result[,"z.inter",]>0&comb.sim.result[,"CP.adap.low",]<=1-beta0,2,sum)
anaps = apply(comb.sim.result[,"N.sim.na",]!=0&comb.sim.result[,"z.inter",]>0&comb.sim.result[,"CP.adap.up",]<=1-beta0,2,sum)
avr2adap = apply(comb.sim.result[,"Vr.2.adap",]<=0,2,sum)
ana; anan; anap; anapg; anapl; anaps; avr2adap

id.check = 4;
atna = which(comb.sim.result[,"N.sim.na",id.check]!=0);
atnan = which(comb.sim.result[,"N.sim.na",id.check]!=0&comb.sim.result[,"z.inter",id.check]<=0)
atnap = which(comb.sim.result[,"N.sim.na",id.check]!=0&comb.sim.result[,"z.inter",id.check]>0)
atnapg = which(comb.sim.result[,"N.sim.na",id.check]!=0&comb.sim.result[,"z.inter",id.check]>0&comb.sim.result[,"CP.adap.low",id.check]>1-beta0)
atnapl = which(comb.sim.result[,"N.sim.na",id.check]!=0&comb.sim.result[,"z.inter",id.check]>0&comb.sim.result[,"CP.adap.low",id.check]<=1-beta0)
atnaps = which(comb.sim.result[,"N.sim.na",id.check]!=0&comb.sim.result[,"z.inter",id.check]>0&comb.sim.result[,"CP.adap.up",id.check]<=1-beta0)
atnaps2l = setdiff(atnapl,atnaps)
atvr2adap = which(comb.sim.result[,"Vr.2.adap",id.check]<=0)

intersect(atvr2adap,atna)
intersect(atvr2adap,atnan)
setdiff(atnan,atvr2adap)## whose simulaitons with z.inter <=0 have no intersection with atvr2adap,
## and NA happens with low Nr.adap and high Nc.adap;
## but it does not matter since we keep the original sample sizes for those simulations

intersect(atvr2adap,atnap)
intersect(atvr2adap,atnapg)
setdiff(atnapg,atvr2adap) ## whose simulaitons with cp.adap.low greater than targeted power have no intersection with atvr2adap,
## and NA happens with low Nr.adap and high Nc.adap with (increasing) monotony partially broken around this area,
## thus those NA do not matter since we set Nr.adap=Nr.min and Nc.adap=Nc.init for those simulaitons.

intersect(atvr2adap,atnapl) ## atnaps=c(atnaps, atnaps2l) and atvr2adap have intersection.
intersect(atvr2adap,atnaps)
setdiff(atnaps, atvr2adap) ## it seems that atnaps \in atvr2adap.
intersect(atvr2adap,atnaps2l) ## atnaps2l and atvr2adap have intersection.
setdiff(atvr2adap,atnaps2l)
setdiff(atnaps2l,atvr2adap)
## so it is atnaps2l and atvr2adap have intersection actually.
## Typically, whose simulaitons in atvr2adap have monotony broken totally.
## while whose simulations in atnaps/atvr2adap have monotony partially broken which happen with low Nr.adap and high Nc.adap,
## and the monotony of the adaptive conditional power near targeted power is remained;
## so it may not be neccessary to delete those simulations.

#### In our simulation, we have:
####    a. the initial number of readers Nr_inint and cases Nc_20 = Nc_21, and maximum number of reader Nr_max and cases Nc_max;
####    b. the number of readers Nr_1 and cases Nc_10 = Nc_11 in the interim analysis where Nc_10 = Nc_20,
####       note that we set the number of normal and disease cases equally without loss of generalization;
####    c. after the interim analsyis, the number of readers Nr.adap can be resized within (Nr_r+1, Nr_max),
####       but the number of cases Nc.adap can only be resized within (Nc_20, Nc_max)--
####       --in the other words, we allow the number of cases to increase only (>=N_20).
#### 1. Based on the above investigation, those simulations with Vr.2.adap non-positive must be deleted,
####    since the estimation of Vr.2.adap is worng caused by the mis-estimation of U-statistics Moments,
####    and the monotony of grid value of conditional power is broken as a whole;
####    (even though some of them don't produce NA, NA can be produced with even higher Nc.adap with relatively low Nr.adap in those cases).
#### 1.1 Recall Vr_2 = c_21*(M_1-M_5)+c_22*(M_2-M_6)+c_23*(M_3-M_7)+c_24*(M_4-M_8) and the definition of c_2x.
####     Theoritically, M_4>M_8, which guarantee that Vr_2 be positive => adaptive conditional variance be positive => no NA happens.
####     Yet in estimation, it is possible that \hat{M_4}-\hat{M_8} is negative (may be due to the U-statistics method);
####     thus causes the estimation of Vr_2 be non-positive => adaptive conditional variance be non-positve => NA happes,
####     especially when the Nr.adap is relative low while the Nc.adap is higher.
#### 1.2 Conditional on \hat{M_4}-\hat{M_8} <=0, when the Nc.adap is higher,
####     c_21, c_22 and c_23 are negligible and the estimation of Vr_2 is donimated be c_24*(M_4-M_8) which is negative.
####     Recall that cond.var = Vr_2/Nr_2+Vc_2-Nr_1*Vc_3^2/(Vr_1+Nr_1*Vc_1) with Vc_3 = Nc_10*Nc_11/(Nc_20*Nc_21)*Vc_1<=Vc_1.
####     Vc_2-Nr_1*Vc_3^2/(Vr_1+Nr_1*Vc_1) is positive, which is guarantted by therory and proved by simulations;
####     however, given that the estimation of Vr_2 is negative and Nr_adap is low, it is possible that the estimation of cond.var <=0,
####     which causes NA grid conditional powers.
#### 2. divide the simualations with NA grid conditional powers into the followings:
####    (2.1 z.inter<0, 2.2 cp.adap.adap.low>1-beta0, 2.3 cp.adap.up<=1-beta0, 2.4 cp.adap.low<=1-beta0&cp.adap.up>1-beta0)
#### 2.1-2.2 Those simulations with NA but also with z.inter<=0, or cp.adap.low>1-beta0 have no intersection with those with Vr.2.adap non-positive;
####      and they need not be deleted according the simulation setting (Nr.adap=Nr.init/Nr.min, Nc.adap=Nc.init/Nc.min) or/and the partially maintained monotony (2.2).
#### 2.3 Those simulations with NA and cp.adap.up<=1-beta0 are typically contained in those with Vr.2.adap non-positive;
####      however, those simulations need not be deleted due to our simulation setting: Nr.adap = Nr.max, Nc.0.adap = Nc.0.max.
#### 2.4 Those simulations with NA and cp.adap.low<=1-beta0 but cp.adap.up> 1-beta0 have intersection with those with Vr.2.adap non-positive;
####      obviously, those intersections need to be deleted;
####      but those without Vr.2.adap non-positive may not need to be deleted since the monotony of the grid conditional power around the targeted power is maintained.
#### 3. NA typically happen when the adaptive conditional variance is non-positive along with low Nr.adap and high Nc.adap:
####    in one situation, they are due to the breaking of monotony becasue of mis-estimated Vr.2.adap;
####    in the other suitaitons, they are due to the partially breaking of monotony with high Nc.adap and low Nr.adap caused by......
#### 4. About the strange points with low cp.inter (<0.25 typically) and low cv.adap (<alpha0):
####    they may or may not accompany with NA grid conditional powers.
####    For those along with NA simulatins, according to our setting, those points should not occur with 2.1 (cv.adap=alpha0) and 2.1 (cp.inter>1-beta0).
####    The majority of grid conditional powers are low (including cp.inter);
####    while minority grid conditional powers increase suddenly and dramatically (and/or produce NA) and pass the targeted power with low Nr.adap and high Nc.adap,
####    this is a sign of monotony broken and low Nr.adap and high Nc.adap cause low cv.adap.
####    Typically, those strang points disappear after deleting simulations with Vr.2.adap non-positive,
####    so they may be due to totoally breaking monotony caused by the mis-estimation of Vr.2.adap
#### 5. Conclusion: after investigate NA, Vr.2.adap and the strange points with low cp.inter/cv.adap, and their relations,
####    it is better to delete only those simulations with Vr.2.adap non-positive except those with cp.adap.up<1-beta0.

#########################################
## plot critical value v.s. conditional power
v.check = comb.sim.result[,,id.check]
id.na.check = unique(which(is.na(v.check),arr.ind = T)[,1]); id.na.check

plot(v.check[,"CP.inter"],v.check[,"cv.adap"])
abline(h=z.alpha,col="red"); abline(v=1-beta0,col="red")
plot(v.check[setdiff(1:n.total, atna),"CP.inter"],v.check[setdiff(1:n.total, atna),"cv.adap"])
abline(h=z.alpha,col="red"); abline(v=1-beta0,col="red")
plot(v.check[setdiff(1:n.total, atvr2adap),"CP.inter"],v.check[setdiff(1:n.total, atvr2adap),"cv.adap"])
abline(h=z.alpha,col="red"); abline(v=1-beta0,col="red")

# plot(v.check[,"z.inter"],v.check[,"cv.adap"])
# plot(v.check[,"z.inter"],v.check[,"CP.inter"])
# plot(v.check[,"z.inter"],v.check[,"Nr.adap"])
# plot(v.check[,"z.inter"],v.check[,"Nc.20.adap"])
# plot(v.check[,"CP.inter"],v.check[,"Nr.adap"])
# plot(v.check[,"CP.inter"],v.check[,"Nc.20.adap"])
# plot(v.check[,"Nr.adap"],v.check[,"Nc.20.adap"])
#
# plot(v.check[v.check[,"z.inter"]<=0,"CP.inter"],v.check[v.check[,"z.inter"]<=0,"cv.adap"])
# plot(v.check[v.check[,"z.inter"]>0,"CP.inter"],v.check[v.check[,"z.inter"]>0,"cv.adap"])
# abline(h=z.alpha,col="red")
# abline(v=1-beta0,col="red")

## check the key function is positive
keyfun = function(nr1, vr1, vc1, vc3adap, sdadap){
  num = (vr1 + nr1*vc1 - nr1*vc3adap)/sdadap
}
kf = array(0, dim=dim(comb.sim.result)[c(1,3)])
for(kk in 1:n.varstr){
  tpt = comb.sim.result[,,kk]
  kf[,kk] = keyfun(SampleSet[kk,"Nr.1"],tpt[,"Vr.1"],tpt[,"Vc.1"],tpt[,"Vc.2.adap"],tpt[,"sd.cond2"])
  print(paste(names.varstr[kk], " where Vc.1 < Vc3.adap:", sep = "")); print(which(tpt[,"Vc.1"]<tpt[,"Vc.2.adap"]))
  print(paste(names.varstr[kk], " where key function <= 0:", sep = "")); print(which(kf[,kk]<=0))
}
#################################################################################

## the mean measurments over the simulations
result.mean = n.effect = n.cplowG = n.cplowG.interRej = n.cplowG.adapNoRej = n.cpupL = n.cpupL.adapRej = n.neg = n.neg.initRej = numeric(0)
for(kk in 1:n.varstr){
  np.tvr2adap = which(comb.sim.result[,"Vr.2.adap",kk]<=0)
  na.cps = which(comb.sim.result[,"N.sim.na",kk]!=0&comb.sim.result[,"z.inter",kk]>0&comb.sim.result[,"CP.adap.up",kk]<=1-beta0)
  na.reject.adap = which(is.na(comb.sim.result[,"reject.adap.cv",kk])) ## 6 NAs happen for Hypo = 30; kk=3; need to check!!
  #id.na = unique(which(is.na(comb.sim.result[,,kk]),arr.ind = T)[,1]) ## id.na == na.reject.adap
  id.effect = setdiff(1:(n.total), union(setdiff(np.tvr2adap,na.cps),na.reject.adap)); n.effect = c(n.effect, length(id.effect))
  result.mean = cbind(result.mean, apply(comb.sim.result[id.effect,,kk],2,mean))
  
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
                            "reject.adap.z", "reject.adap.t", "reject.adap.cv", "Nr.adap", "Nc.20.adap", "n.effect",
                            "n.cplowG", "n.cplowG.interRej", "n.cplowG.adapNoRej", "n.cpupL", "n.cpupL.adapRej", "n.neg", "n.neg.initRej"),], n.digit))
}else{
  print(round(rbind(result.mean[c("reject.inter.z", "reject.inter.t", "reject.init.z", "reject.init.t",
                                  "reject.init.theory", "reject.min.theory", "reject.max.theory", "reject.adap.cv"),],
                    t(SampleSet[,c("Nr.1","Nr.init","Nr.expect","Nc.20")]),result.mean[c("Nr.adap","Nc.20.adap","n.effect", 
                                                                                         "n.cplowG", "n.cplowG.interRej", "n.cplowG.adapNoRej",
                                                                                         "n.cpupL", "n.cpupL.adapRej",
                                                                                         "n.neg", "n.neg.initRej"),]), n.digit))
}

## calculate the simulations where iMRMC fail with negative var(AUCAminusAUCB)
apply(is.na(comb.sim.result[,"sd.adap",]),2,sum); #which(is.na(comb.sim.result[,"sd.adap",3]))
apply(comb.sim.result[,"varUstat.adap.na",],2,sum);
apply(comb.sim.result[,"varMLE.adap.na",],2,sum);

time.summ = proc.time() - time.summ

#######################################################################################
## save the simulation summary
fname = ifelse(Hypo == "null", paste("summary_", AdapMethod, "_", Hypo, "_", Ratio.Nc0toNc1, ".Rdata", sep=""),
               paste("summary_", AdapMethod, "_", Hypo, "_", Ratio.Nc0toNc1, "_", Power0, ".Rdata", sep=""))
rm(setseed, kk, sim.result, time.sim, comb.sim.result0, comb.sim.result, id.check, v.check, tpt)
save.image(file = paste(sumdir,"/", fname,sep=""))
#load(fname)

#######################################################################################



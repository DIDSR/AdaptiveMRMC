#########################################################################################################################
#' Demonstration of adaptive method 1 with simulation data
#' The method is presented in the following reference
#' Huang Z, Samuelson F, Tcheuko L, Chen W.
#' Adaptive designs in multi-reader multi-case clinical trials of imaging devices.
#' Statistical Methods in Medical Research. 2020;29(6):1592-1611. doi:10.1177/0962280219869370
#'
library(iMRMC)
"%+%" = function(x,y) paste(x,y,sep="")

print("This is a demonstration of adaptive method 1 with simulation data:re-sizing the readers after an interim analysis")
print("The parameters (e.g., for simulation) are chosen arbitrarily and this is only to demo how the method works.")
print("#########################################################################")
print("STEP1: Sample size setting")
print("Nr.max = 30: maximum number of readers the study can possibly afford")
print("Nc.0.max = 250: maximum number of non-diseased patients whose images are available for the study")
print("Nc.1.max = 250: maximum number of diseased patients whose images are available for the study")
print("Assuming initial sizing: Nr.init = 8, Nc.0 = 100, Nc.1 = 100")
print("In adaptive method 1, we only resize the number of readers")
print("Interim analysis will be conducted after collecting Nr.1 = 3 readers' data")
Nr.max = 30; Nc.0.max = 250; Nc.1.max = 250
Nr.init = 8; Nr.1 = 3; Nc.0 = 100; Nc.1 = 100;
readline(prompt="Press [enter] to continue")

#indir = "C:/Users/WXC4/OneDrive - FDA/Documents/Research/AdaptiveMRMC/ZhipengHuang/MRMC/Rpackage/Testing";
Para = "para."; Tbl = "Table1";
#load(indir %+% "/ParameterSetting/ParameterSettingFor" %+% Tbl %+% ".Rdata")
print("#######################################################################")
print("STEP2: Collection of data before the interim analysis. Here we use simulat ion to generate data for demo purpose.")
load("data/ParameterSettingForTable1.Rdata")
ParaSet = eval(as.name(Para %+% "alter50"))
set.seed(2022);
dFrame.imrmc = mrmcRMscoresFC(Nr.max, Nc.0.max, Nc.1.max, ParaSet[5,])

data.interim = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.0), paste0("posCase", 1:Nc.1)) &
                        readerID %in% c(paste0("reader", 1:Nr.1), "-1"))
samples <- list("Nr.init" = Nr.init, "Nr.1" = Nr.1, "Nc.0" = Nc.0, "Nr.max" = Nr.max, "R.Nc0toNc1" = 1)
readline(prompt="Press [enter] to continue")

print("######################################################################")
print("STEP3: Interim analysis.")
print("Please wait...")
interim.adap = adaptiveI(data.interim, samples, 1)
print("Interium analysis results:")
print(sprintf("The conditional power in the interim analysis was found to be %f", interim.adap$CP.inter))
print(sprintf("The reader sample size is adaptively re-sized to %i to achieve a power of %f", interim.adap$nr.adap,interim.adap$CP.adap))
print(sprintf("The critical value for z test in the final analysis is computed as %f", round(interim.adap$cv.adap, 4)))

readline(prompt="Press [enter] to continue")


print("######################################################################")
print("STEP4: Collection of data after interim analysis and Final analysis.")
print("Please wait...")
data.adap = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.0), paste0("posCase", 1:Nc.1)) &
                     readerID %in% c(paste0("reader", 1:interim.adap$nr.adap), "-1"))
randseed = 2
wkdir = paste(getwd(), "/", randseed+3, sep = ""); dir.create(wkdir)
BDG.0 = doIMRMC(data = data.adap, workDir = wkdir)
BDG.adap = list(Ustat = BDG.0$Ustat,
                moments = rbind(BDG.0$varDecomp$BDG$Ustat$comp$modalityA.modalityB,
                                BDG.0$varDecomp$BDG$Ustat$coeff$modalityA.modalityB[1,]),
                MLEstat = BDG.0$MLEstat)
row.names(BDG.adap$moments) = c("modalityA", "modalityB", "crossAB", "coeff")
rm(BDG.0); unlink(wkdir, recursive = TRUE)

dAUC.adap = BDG.adap$Ustat$AUCAminusAUCB[3];
var.dAUC.adap = BDG.adap$Ustat$varAUCAminusAUCB[3]
z.adap = dAUC.adap/sqrt(var.dAUC.adap)
reject.adap.cv = as.numeric(z.adap > abs(interim.adap$cv.adap))

print("The final analysis results")
print(paste0("Num of readers = ", interim.adap$nr.adap, ", num of non-diseased cases = ", Nc.0, ", num of diseased cases = ", Nc.1))
print(paste0("The difference of AUC (SD) = ", round(dAUC.adap, 4), "(", round(sqrt(var.dAUC.adap),4), ")"))
print(paste0("The z statistic = ", round(z.adap,5), ", and the critical value obtained in the interium analysis = ", round(abs(interim.adap$cv.adap), 4)))
print(paste0("The study conclusion is: ", ifelse(reject.adap.cv,
                                                "the difference of AUC is statistically significant.",
                                                "the difference of AUC is not statistically significant.")))

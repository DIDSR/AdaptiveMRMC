## create .swa file to run parallel using swarm
rm(list=ls())

## set the work directory
setwd("/home/zhuang/Desktop/MRMC/Rpackage/Testing")
#setwd("C:/Users/Zhipeng.Huang/Desktop/MRMC/Rpackage/Testing")

N.paral = 1000;
AdapMethod = "Adaptive1"; AdapMethod = "Adaptive2"
Hypo = "null"; Hypo = "alter50";
Ratio.Nc0toNc1 = 1; Ratio.Nc0toNc1 = 2
Power0 = NA; Power0 = "under"; #Power0 = "over";

name.Rscript = "runsimAdaptives.R"

filename = paste("run.", AdapMethod, ".", Hypo, ".", Ratio.Nc0toNc1, ".", Power0, ".swa", sep="")
if(Hypo == "null") filename = paste("run.", AdapMethod, ".", Hypo, ".", Ratio.Nc0toNc1, ".swa", sep="")

sink(file = filename, append = TRUE)

for(i in 1:N.paral){
  cat('export DISPLAY=\":0.0\";', "Rscript", name.Rscript, AdapMethod, Hypo, Ratio.Nc0toNc1, Power0, i, '\n')
}

sink()


#' Implement simulation evaluation of adaptive method II where only both reader sample and case samples are resized
#'
#' This is the main function to implement simulations to evaluate the adaptive method II where both the reader size and case sizes are adapted. Right-sided hypothesis test is considered
#'
#' @param samples [list] of sample sizes:
#'
#'               $Nr.init:    [num] initial total number of readers
#'
#'               $Nr.1:       [num] number of readers for the first part of study (before the interim analysis)
#'
#'               $Nr.max:     [num] maximum number of readers
#'
#'               $Nc.10:      [num] number of non-diseased cases for the first part of the study (before the interim analysis)
#'
#'               $Nc.20:      [num] initial number of non-diseased cases for the second part of the study (after the interim analysis, and we set $Nc.20=$Nc.10)
#'
#'               $Nc.0.max:   [num] maximum number of non-diseased cases (to generate the same data as adaptive method II and make comparison)
#'
#'               $R.Nc0toNc1: [num] the ratio of the number of non-diseased cases to that of diseased cases
#'
#'               ($Nr.expect: [num] expected number of adaptive readers to reach targeted power)
#' @param paras [list] of parameters:
#'
#'               $mu1: [num] grand average
#'
#'               $tauA1: [num] modality effect (for the signal-present class)
#'
#'               $b, $Vc0, $Vtauc0, $Vrc0, $Veps0, $Vr0, $Vtaur0: [num] variance component parameters (see references)
#'
#'               ($varStruc: [factor] variance structure; $AUC_B: [num] base line AUC)
#' @param randseed [num] the random seed to generate data
#' @param alpha0 [num] the nominal type I error
#' @param beta0  [num] 1-beta0 is the targeted power
#'
#' @return  result.interim [list] with the updated reader/cases sizes, critical value, test decisions and other results:
#'
#'              $nr.adap:       [num] updated totoal reader sample size
#'
#'              $nc.20.adap:    [num] updated non-diseased case sample size for the second part of the study (the updated diseased case sample size is nc.21.adap = $nc.20.adap/$R.Nc0toNc1)
#'
#'              $cv.adap:       [num] updated critical value
#'
#'              $reject.init.z  [num] test decision using wald-test w.r.t. initial sample without interim analysis
#'
#'              $reject.init.z  [num] test decision using t-test w.r.t. initial sample without interim analysis
#'
#'              $reject.adap.z  [num] test decision using wald-test w.r.t. adjusted reader size and conventional critical value
#'
#'              $reject.adap.t  [num] test decision using t-test w.r.t. adjusted reader size and conventional critical value
#'
#'              $reject.adap.cv [num] test decision using wald-test w.r.t. adjusted reader size and updated critical value
#'
#'              $...: other intermediate and final results
#'
#' @export
## @import iMRMC
#' @importFrom iMRMC doIMRMC

sim.adaptiveII = function(samples, paras, randseed=666, alpha0 = 0.025, beta0 = 0.2){
  ## set the quantile z.alpha and z.beta w.r.t. type I error and targeted power
  z.alpha = qnorm(1-alpha0); z.beta = qnorm(1-beta0)

  ## sample size setting
  Nr.init = samples$Nr.init; Nr.1 = samples$Nr.1; Nr.max = samples$Nr.max;
  Nc.10 = samples$Nc.10; Nc.11 = round(Nc.10/samples$R.Nc0toNc1);
  Nc.20 = samples$Nc.20; Nc.21 = round(Nc.20/samples$R.Nc0toNc1);
  Nc.0.max = samples$Nc.0.max; Nc.1.max = round(Nc.0.max/samples$R.Nc0toNc1)

  ########################### Adaptive Design Method I Key Component ###############################
  ## generate the complete data using RM model
  set.seed(randseed+1); dFrame.imrmc = mrmcRMscoresFC(Nr.max, Nc.0.max, Nc.1.max, paras)

  ## the data for interim analysis, which is a part of the full data
  data.interim = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.10), paste0("posCase", 1:Nc.11)) &
                    readerID %in% c(paste0("reader", 1:Nr.1), "-1"))

  ## run BDG Ustat model for interim analysis and adjuste the reader/cases sample sizes and critical value
  ## the updated diseased case sample size is nc.21.adap = $nc.20.adap/$R.Nc0toNc1
  interim.adap = adaptiveII(data.interim, samples, randseed)
  Nr.adap =  interim.adap$nr.adap; Nc.20.adap =  interim.adap$nc.20.adap;
  Nc.21.adap = round(Nc.20.adap/samples$R.Nc0toNc1); cv.adap = interim.adap$cv.adap
  z.1 = interim.adap$z.inter
  ##################################################################################################

  ######## the adaptive data, which is a part of the full data ########
  if(Nc.20.adap == Nc.10){
    data.t0.adap = NULL
  }else{
    data.t0.adap = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% paste0("negCase", (Nc.10+1):Nc.20.adap) &
                            readerID =="-1")
  }
  if(Nc.21.adap == Nc.11){
    data.t1.adap = NULL
  }else{
    data.t1.adap = subset(dFrame.imrmc$dFrame.imrmc, caseID %in%  paste0("posCase", (Nc.11+1):Nc.21.adap) &
                            readerID =="-1")
  }
  data.2.adap = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.20.adap), paste0("posCase", 1:Nc.21.adap)) &
                         readerID %in% paste0("reader", (Nr.1+1):Nr.adap))
  data.adap = rbind(data.interim, data.t0.adap, data.t1.adap, data.2.adap)

  ## run BDG Ustat model for adaptive analysis
  #BDG.adap = doIMRMCfast(data.adap, randseed+3)
  wkdir = paste(getwd(), "/", randseed+3, sep = ""); dir.create(wkdir)
  BDG.0 = doIMRMC(data = data.adap, workDir = wkdir)
  BDG.adap = list(Ustat = BDG.0$Ustat,
                  moments = rbind(BDG.0$varDecomp$BDG$Ustat$comp$modalityA.modalityB,
                                  BDG.0$varDecomp$BDG$Ustat$coeff$modalityA.modalityB[1,]),
                  MLEstat = BDG.0$MLEstat)
  row.names(BDG.adap$moments) = c("modalityA", "modalityB", "crossAB", "coeff")
  rm(BDG.0); unlink(wkdir, recursive = TRUE)

  dAUC.adap = BDG.adap$Ustat$AUCAminusAUCB[3]
  sd.adap = ifelse(BDG.adap$Ustat$varAUCAminusAUCB[3]<=0,
                   sqrt(BDG.adap$MLEstat$varAUCAminusAUCB[3]),
                   sqrt(BDG.adap$Ustat$varAUCAminusAUCB[3]))
  z.adap = dAUC.adap/sd.adap
  varUstat.adap.na = ifelse(BDG.adap$Ustat$varAUCAminusAUCB[3]<=0, 1, 0)
  varMLE.adap.na = ifelse(BDG.adap$MLEstat$varAUCAminusAUCB[3]<=0, 1, 0)

  ######## the data for initial analysis, which is a part of the full data ########
  if(z.1<=0){
    data.init = data.adap; BDG.init = BDG.adap
    dAUC.init = dAUC.adap; sd.init = sd.adap; z.init = z.adap
  }else{
    data.init = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.10), paste0("posCase", 1:Nc.11)) &
                         readerID %in% c(paste0("reader", 1:Nr.init), "-1"))
    ## run BDG Ustat model for initial analysis
    #BDG.init = doIMRMCfast(data.init, randseed+4)
    wkdir = paste(getwd(), "/", randseed+4, sep = ""); dir.create(wkdir)
    BDG.0 = doIMRMC(data = data.init, workDir = wkdir)
    BDG.init = list(Ustat = BDG.0$Ustat,
                    moments = rbind(BDG.0$varDecomp$BDG$Ustat$comp$modalityA.modalityB,
                                    BDG.0$varDecomp$BDG$Ustat$coeff$modalityA.modalityB[1,]),
                    MLEstat = BDG.0$MLEstat)
    row.names(BDG.init$moments) = c("modalityA", "modalityB", "crossAB", "coeff")
    rm(BDG.0); unlink(wkdir, recursive = TRUE)

    dAUC.init = BDG.init$Ustat$AUCAminusAUCB[3];
    sd.init = sqrt(BDG.init$Ustat$varAUCAminusAUCB[3])
    z.init = dAUC.init/sd.init
  }

  ######## test decisions based on normal/t distribution and differenct critical values ########
  reject.inter.z = interim.adap$reject.inter.z;  interim.adap$reject.inter.z = NULL
  reject.inter.t = interim.adap$reject.inter.t;  interim.adap$reject.inter.t = NULL

  reject.init.z = as.numeric(z.init > z.alpha)
  reject.init.t = as.numeric(BDG.init$Ustat$rejectBDG[3] && z.init > 0)

  reject.adap.z = as.numeric(z.adap > z.alpha)
  reject.adap.t = as.numeric(ifelse(BDG.adap$Ustat$varAUCAminusAUCB[3]<=0,
                                      BDG.adap$MLEstat$rejectBDG[3],BDG.adap$Ustat$rejectBDG[3]) && z.adap > 0)
  reject.adap.cv = as.numeric(z.adap > abs(cv.adap))


  return(c(as.numeric(interim.adap),
           dAUC.init, sd.init, z.init, dAUC.adap, sd.adap, z.adap,
           varUstat.adap.na, varMLE.adap.na,
           reject.inter.z, reject.inter.t, reject.init.z, reject.init.t,
           reject.adap.z, reject.adap.t, reject.adap.cv))
  # return(data.frame(interim.adap,
  #          dAUC.init, sd.init, z.init, dAUC.adap, sd.adap, z.adap,
  #          varUstat.adap.na, varMLE.adap.na,
  #          reject.inter.z, reject.inter.t, reject.init.z, reject.init.t,
  #          reject.adap.z, reject.adap.t, reject.adap.cv))
}

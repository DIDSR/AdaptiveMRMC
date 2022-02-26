#' Implement simulation evaluation of adaptive method I where only reader sample is resized
#'
#' This is the main function to implement simulations to evaluate the adaptive method I where only the reader size is adapted. The cases sizes remain unchanged. Right-sided hypothesis test is considered
#'
#' @param samples [list] of sample sizes:
#'
#'               $Nr.init:    [num] initial total number of readers
#'
#'               $Nr.1:       [num] number of readers for the first part of study (before the interim analysis)
#'
#'               $Nr.max:     [num] maximum number of readers
#'
#'               $Nc.0:       [num] fixed number of non-diseased cases
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
#' @return  result.interim [list] with the updated reader size, critical value, test decisions and other results:
#'
#'              $nr.adap:       [num] updated total reader sample size
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
#' @references Huang Z, Samuelson F, Tcheuko L, Chen W.Statistical Methods in Medical Research. 2020;29(6):1592-1611.
## @import iMRMC
#' @importFrom iMRMC doIMRMC

sim.adaptiveI = function(samples, paras, randseed=666, alpha0 = 0.025, beta0 = 0.2){
  ## set the quantile z.alpha and z.beta w.r.t. type I error and targeted power
  z.alpha = qnorm(1-alpha0); z.beta = qnorm(1-beta0)

  ## sample size setting
  Nr.init = samples$Nr.init; Nr.1 = samples$Nr.1; Nr.max = samples$Nr.max;
  Nc.0 = samples$Nc.0; Nc.1 = round(Nc.0/samples$R.Nc0toNc1);
  Nc.0.max = samples$Nc.0.max; Nc.1.max = round(Nc.0.max/samples$R.Nc0toNc1)

  ########################### Adaptive Design Method I Key Component ###############################
  ## generate the complete data using RM model
  set.seed(randseed+1); dFrame.imrmc = mrmcRMscoresFC(Nr.max, Nc.0.max, Nc.1.max, paras)

  ## the data for interim analysis, which is a part of the full data
  data.interim = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.0), paste0("posCase", 1:Nc.1)) &
                    readerID %in% c(paste0("reader", 1:Nr.1), "-1"))

  ## run BDG Ustat model for interim analysis and adjust the reader sample size and critical value
  interim.adap = adaptiveI(data.interim, samples, randseed)
  Nr.adap =  interim.adap$nr.adap; cv.adap = interim.adap$cv.adap
  z.1 = interim.adap$z.inter
  ##################################################################################################

  ## the adaptive data, which is a part of the full data
  ## and run BDG Ustat model for interim analysis
  data.adap = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.0), paste0("posCase", 1:Nc.1)) &
                       readerID %in% c(paste0("reader", 1:Nr.adap), "-1"))

  #BDG.adap = doIMRMCfast(data.adap, randseed+3)
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

  ## the data for initial analysis, which is a part of the full data
  if(z.1<=0){
    data.init = data.adap; BDG.init = BDG.adap;
    dAUC.init = dAUC.adap; var.dAUC.init = var.dAUC.adap; z.init = z.adap
  }else{
    data.init = subset(dFrame.imrmc$dFrame.imrmc, caseID %in% c(paste0("negCase", 1:Nc.0), paste0("posCase", 1:Nc.1)) &
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
    var.dAUC.init = BDG.init$Ustat$varAUCAminusAUCB[3]
    z.init = dAUC.init/sqrt(var.dAUC.init)
  }

  ## test decisions based on normal/t distribution and different critical values
  reject.inter.z = interim.adap$reject.inter.z;  interim.adap$reject.inter.z = NULL
  reject.inter.t = interim.adap$reject.inter.t;  interim.adap$reject.inter.t = NULL

  reject.init.z = as.numeric(z.init > z.alpha);
  reject.init.t= as.numeric(BDG.init$Ustat$rejectBDG[3] && z.init > 0)

  reject.adap.z = as.numeric(z.adap > z.alpha);
  reject.adap.t= as.numeric(BDG.adap$Ustat$rejectBDG[3] && z.adap > 0)
  reject.adap.cv = as.numeric(z.adap > abs(cv.adap));

  return(c(as.numeric(interim.adap),dAUC.init, var.dAUC.init, z.init, dAUC.adap, var.dAUC.adap, z.adap,
           reject.inter.z, reject.inter.t, reject.init.z, reject.init.t,
           reject.adap.z, reject.adap.t,reject.adap.cv))
  # return(data.frame(interim.adap, dAUC.init, var.dAUC.init, z.init, dAUC.adap, var.dAUC.adap, z.adap,
  #          reject.inter.z, reject.inter.t, reject.init.z, reject.init.t,
  #          reject.adap.z, reject.adap.t,reject.adap.cv))
}

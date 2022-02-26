#' Adaptive method I where only reader sample is resized
#'
#' The main function to implement the adaptive method I where only the reader size is adapted after the data for the interim analysis has been collected.
#' The critical value and the reader sample size are updated simultaneously. The data for interim analysis has been converted to doIMRMC formatted data frame. The cases sizes remain unchanged. Right-sided hypothesis test is considered.
#'
#'
#' @param data.interim [data.frame] with (nc0+nc1)+nr.1*(nc0+nc1)*2 rows and 4 variables including:
#'
#'              $readerID: [Factor] w/ nr.1+1 levels "-1", "reader1", "reader2", ...
#'
#'              $caseID: [Factor] w/ nc0+nc1 levels "negCase1", "negCase2", ..., posCase1", "posCase2", ...
#'
#'              $modalityID: [Factor] w/ 3 levels "truth", "modility A", "modility B"
#'
#'              $score: [num] reader score
#' @param samples  [list] of sample sizes:
#'
#'               $Nc.0:       [num] fixed number of non-diseased cases
#'
#'               ($Nc.1):     [num] fixed diseased cases size
#'
#'               $R.Nc0toNc1: [num] the ratio of the number of non-diseased cases to that of diseased cases (we assume $R.Nc0toNc1>=1)
#'
#'               $Nr.init:    [num] initial total number of readers
#'
#'               $Nr.1:       [num] number of readers for the first part of study (before the interim analysis)
#'
#'               $Nr.max:     [num] maximum number of readers
#'
#'               ($Nc.0.max:   [num] maximum number of non-diseased cases (to generate the same data as adaptive method II and make comparison); $Nr.expect: [num] expected number of adaptive readers to reach targeted power)
#' @param randseed    [num] the random seed to indicate the particular doIMRMCfast() intermediate folder/file
#' @param alpha0      [num] the nominal type I error
#' @param beta0       [num] 1-beta0 is the targeted power
#'
#' @return  result.interim [list] with the updated reader size and critical value, and other intermediate results:
#'
#'              $nr.adap: [num] updated total reader sample size
#'
#'              $cv.adap: [num] updated critical value
#'
#'              $...: other intermediate results
#'
#' @export
#' @references Huang Z, Samuelson F, Tcheuko L, Chen W.Statistical Methods in Medical Research. 2020;29(6):1592-1611.
#'
#' @importFrom iMRMC doIMRMC

adaptiveI = function(data.interim, samples, randseed=666, alpha0 = 0.025, beta0 = 0.2){
  ## set the quantile z.alpha and z.beta w.r.t. type I error and targeted power
  z.alpha = qnorm(1-alpha0); z.beta = qnorm(1-beta0)

  ## sample size setting
  nr.init = samples$Nr.init; nr.1 = samples$Nr.1; nr.2 = nr.init-nr.1; nr.max = samples$Nr.max;
  nc.0 = samples$Nc.0; nc.1 = ifelse(!is.null(samples$Nc.1), samples$Nc.1, round(nc.0/samples$R.Nc0toNc1))

  ## calculate {c1,c2,c3,c4}
  vector.c = fun.vector.c(nc.0, nc.1);
  names(vector.c) =  paste("c", 1:length(vector.c), sep = "")

  ## run BDG Ustat model for interim analysis
  #BDG.1 = doIMRMCfast(data.interim, randseed+2)
  wkdir = paste(getwd(), "/", randseed+2, sep = ""); dir.create(wkdir)
  BDG.0 = doIMRMC(data = data.interim, workDir = wkdir)
  BDG.1 = list(Ustat = BDG.0$Ustat,
               moments = rbind(BDG.0$varDecomp$BDG$Ustat$comp$modalityA.modalityB,
                               BDG.0$varDecomp$BDG$Ustat$coeff$modalityA.modalityB[1,]),
               MLEstat = BDG.0$MLEstat)
  row.names(BDG.1$moments) = c("modalityA", "modalityB", "crossAB", "coeff")
  rm(BDG.0); unlink(wkdir, recursive = TRUE)

  ## AUC difference of Modality A and B, the varaince of difference and the z-statistics
  dAUC.1 = BDG.1$Ustat$AUCAminusAUCB[3];
  var.dAUC.1 = BDG.1$Ustat$varAUCAminusAUCB[3]
  z.1 = dAUC.1/sqrt(var.dAUC.1)

  ## moments of the AUC difference
  dAUC.1.moments = BDG.1$moments[1,]+BDG.1$moments[2,] - 2*BDG.1$moments[3,]
  rownames(dAUC.1.moments) = NULL

  ## estimate Vr and Vc from interim analysis # and: var.dAUC.1 = 1/nr.1*Vr.1+Vc.1
  VrVc =  fun.VrVc(dAUC.1.moments, vector.c); Vr.1 = VrVc$Vr; Vc.1 = VrVc$Vc; rm(VrVc)
  Vcr.Ratio.1 = Vc.1/Vr.1

  ## conditional type I error and conditional power given the estimates in interim analysis for right sided test
  CE.inter = 1-pnorm(z.alpha*sqrt(nr.init/(nr.2)*(1+nr.1*Vcr.Ratio.1))-z.1*sqrt(nr.1/nr.2*(1+nr.init*Vcr.Ratio.1)))
  CP.inter = 1-pnorm(z.alpha*sqrt(nr.init/(nr.2)*(1+nr.1*Vcr.Ratio.1))-z.1*nr.init*(1+nr.1*Vcr.Ratio.1)/sqrt(nr.1*(nr.2)*(1+nr.init*Vcr.Ratio.1)))

  ## calculate the limit of adjusted conditional power
  ## when z.1>0 (if z.1<0, CP.adap.low->CP.adap.up and CP.adap.up->CP.adap.low)
  CP.adap.low = pnorm(z.1/sqrt(nr.1*(1+(nr.1+1)*Vcr.Ratio.1))
                      +z.1*sqrt(nr.1/nr.2*(1+nr.init*Vcr.Ratio.1))
                      -z.alpha*sqrt(nr.init/nr.2*(1+nr.1*Vcr.Ratio.1)))
  CP.adap.up = pnorm(z.1/sqrt(nr.1)*sqrt((nr.max-nr.1)/(1+nr.max*Vcr.Ratio.1))
                     +z.1*sqrt(nr.1/nr.2*(1+nr.init*Vcr.Ratio.1))
                     -z.alpha*sqrt(nr.init/nr.2*(1+nr.1*Vcr.Ratio.1)))

  ## adjust the number of readers
  if(z.1<=0){
    nr.adap = nr.init
    s0 = z.alpha*sqrt(nr.init/nr.2*(1+nr.1*Vcr.Ratio.1))-z.1*sqrt(nr.1/nr.2*(1+nr.init*Vcr.Ratio.1))+z.beta
    s1 = sqrt(nr.1)/z.1*s0
  }else{
    if(1-beta0<min(CP.adap.low,CP.adap.up)){
      s0 = s1 = 999
      nr.adap = nr.1+1
    }else if(1-beta0>max(CP.adap.low,CP.adap.up)){
      s0 = s1 = -99
      nr.adap = nr.max
    }else{
      ## resize the number of readers in the interim analysis
      ## and adjust the corresponding critical value
      ## given AUC, Vr and Vc estimated from the interim analysis
      s0 = z.alpha*sqrt(nr.init/nr.2*(1+nr.1*Vcr.Ratio.1))-z.1*sqrt(nr.1/nr.2*(1+nr.init*Vcr.Ratio.1))+z.beta
      s1 = sqrt(nr.1)/z.1*s0
      nr.adap = ceiling((s1^2+nr.1)/(1-Vcr.Ratio.1*s1^2))
    }
  }

  nr.2.adap = nr.adap - nr.1

  ## conditonal power after the number of readers is adjusted
  CP.adap = pnorm(z.1/sqrt(nr.1)*sqrt(nr.2.adap/(1+nr.adap*Vcr.Ratio.1))
                  +z.1*sqrt(nr.1/nr.2*(1+nr.init*Vcr.Ratio.1))
                  -z.alpha*sqrt(nr.init/nr.2*(1+nr.1*Vcr.Ratio.1)))

  ## adjusted critical value
  cv.adap.num = (z.1*(sqrt(nr.1/nr.2.adap*(1+nr.adap*Vcr.Ratio.1))-sqrt(nr.1/nr.2*(1+nr.init*Vcr.Ratio.1)))
                 +z.alpha*sqrt(nr.init/nr.2*(1+nr.1*Vcr.Ratio.1)))
  cv.adap.den = sqrt(nr.adap/nr.2.adap*(1+nr.1*Vcr.Ratio.1))
  cv.adap = cv.adap.num/cv.adap.den; rm(cv.adap.num,cv.adap.den)

  ## test decisions based on normal/t distribution in the interim analysis
  reject.inter.z = as.numeric(z.1 > z.alpha);
  reject.inter.t = as.numeric(BDG.1$Ustat$rejectBDG[3] && z.1 > 0)

  ## return the outputs
  result.interim1 = data.frame(nr.adap, cv.adap)
  result.interim2 = data.frame(Vr.1, Vc.1, Vcr.Ratio.1, dAUC.1, var.dAUC.1, z.1,
                              CE.inter, CP.inter, CP.adap.low, CP.adap.up, CP.adap, s0, s1,
                              reject.inter.z, reject.inter.t)
  names(result.interim2) = c("Vr.inter", "Vc.inter", "Vcr.Ratio.inter", "dAUC.inter", "var.dAUC.inter", "z.inter",
                            "CE.inter", "CP.inter",  "CP.adap.low", "CP.adap.up", "CP.adap", "s0", "s1",
                            "reject.inter.z", "reject.inter.t")
  result.interim = data.frame(result.interim1, dAUC.1.moments, result.interim2)
  return(result.interim)
}

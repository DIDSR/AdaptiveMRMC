#' Adaptive method II where both reader sample and case samples are resized
#'
#' The main function to implement the adaptive method II where both the reader and the cases (non-diseased and diseased) sizes are adapted after the data for the interim analysis has been collected.
#' The critical value and the reader/cases sample sizes are updated simultaneously. The data for interim analysis has been converted to doIMRMC formatted data frame.  Right-sided hypothesis test is considered.
#' Please refer to xxx paper for more detail.
#'
#' @param data.interim [data.frame] with (nc10+nc11)+nr.1*(nc10+nc11)*2 rows and 4 variables including:
#'
#'              $readerID: [Factor] w/ nr.1+1 levels "-1", "reader1", "reader2", ...
#'
#'              $caseID: [Factor] w/ nc10+nc11 levels "negCase1", "negCase2", ..., posCase1", "posCase2", ...
#'
#'              $modalityID: [Factor] w/ 3 levels "truth", "modility A", "modility B"
#'
#'              $score: [num] reader score
#' @param samples  [list] of sample sizes:
#'
#'               $Nc.10:      [num] number of non-diseased cases for the first part of the study (before the interim analysis)
#'
#'               ($Nc.11):    [num] number of diseased cases size for the first part of the study (before the interim analysis)
#'
#'               $R.Nc0toNc1: [num] the ratio of the number of non-diseased cases to that of diseased cases (we assume $R.Nc0toNc1>=1)
#'
#'               $Nc.20:      [num] initial number of non-diseased cases for the second part of the study (after the interim analysis, and we set $Nc.20=$Nc.10)
#'
#'               ($Nc.21):    [num] initial number of diseased cases for the second part of the study (after the interim analysis, , and we set $Nc.21=$Nc.11)
#'
#'               $Nc.0.max:   [num] maximum number of non-diseased cases (to generate the same data as adaptive method I and make comparison)
#'
#'               $Nr.init:    [num] initial total number of readers
#'
#'               $Nr.1:       [num] number of readers for the first part of study (before the interim analysis)
#'
#'               $Nr.max:     [num] maximum number of readers
#'
#'               ($Nr.expect: [num] expected number of adaptive readers to reach targeted power when fixing the cases sizes)
#' @param randseed    [num] the random seed to indicate the particular doIMRMCfast() intermediate folder/file
#' @param alpha0      [num] the nominal type I error
#' @param beta0       [num] 1-beta0 is the targeted power
#'
#' @return  result.interim [list] with the updated reader/cases sizes and critical value, and other intermediate results:
#'
#'              $nr.adap:       [num] updated totoal reader sample size
#'
#'              $nc.20.adap:    [num] updated non-diseased case sample size for the second part of the study (the updated diseased case sample size is nc.21.adap = $nc.20.adap/$R.Nc0toNc1)
#'
#'              $cv.adap:       [num] updated critical value
#'
#'              $...: other intermediate results
#'
#' @export
## @import iMRMC
#' @importFrom iMRMC doIMRMC


adaptiveII = function(data.interim, samples, randseed=666, alpha0 = 0.025, beta0 = 0.2){
  ## set the quantile z.alpha and z.beta w.r.t. type I error and targeted power
  z.alpha = qnorm(1-alpha0); z.beta = qnorm(1-beta0)

  ## sample size setting
  nr.init = samples$Nr.init; nr.1 = samples$Nr.1; nr.2 = nr.init-nr.1; nr.max = samples$Nr.max;
  nc.10 = samples$Nc.10; nc.11 = ifelse(!is.null(samples$Nc.11), samples$Nc.11, round(nc.10/samples$R.Nc0toNc1))
  nc.20 = ifelse(!is.null(samples$Nc.20), samples$Nc.20, samples$Nc.10) ;
  nc.21 = ifelse(!is.null(samples$Nc.21), samples$Nc.21, ifelse(!is.null(samples$Nc.11), samples$Nc.11, round(nc.20/samples$R.Nc0toNc1))) ;
  nc.0.max = samples$Nc.0.max

  ## calculate {c1,c2,c3,c4}
  vector.c1 = fun.vector.c(nc.10, nc.11)
  names(vector.c1) =  paste("c1", 1:length(vector.c1), sep = "")
  vector.c2 = fun.vector.c(nc.20, nc.21)
  names(vector.c2) =  paste("c2", 1:length(vector.c2), sep = "")

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

  ## moments of the AUC difference from the interim analysis
  dAUC.moments.1 = BDG.1$moments[1,]+BDG.1$moments[2,] - 2*BDG.1$moments[3,]
  rownames(dAUC.moments.1) = NULL

  ## estimate Vr1, Vr2, Vc1, Vc2 and Vc3 from interim analysis
  VrVc1 =  fun.VrVc(dAUC.moments.1, vector.c1); VrVc2 =  fun.VrVc(dAUC.moments.1, vector.c2);
  Vr.1 = VrVc1$Vr; Vc.1 = VrVc1$Vc; rm(VrVc1);
  Vr.2 = VrVc2$Vr; Vc.2 = VrVc2$Vc; rm(VrVc2)
  #Vc.3 = Nc.10*Nc.11/(Nc.20*Nc.21)*Vc.1 ## This derivation of Vc.3 is wrong
  Vc.3 = Vc.2

  ## estimate variances from the interim analysis
  sd.1 = sqrt(Vr.1/nr.1+Vc.1) # or sd.1 =  sqrt(BDG.1$Ustat$varAUCAminusAUCB[3])
  sd.cond2 = sqrt(adaptiveII.condvarR2.varR(nr.1, nr.2, Vr.1, Vc.1, Vr.2, Vc.2, Vc.3)$v.condR2)
  sd.total = sqrt(adaptiveII.condvarR2.varR(nr.1, nr.2, Vr.1, Vc.1, Vr.2, Vc.2, Vc.3)$v.R)

  ## AUC difference of Modality A and B, and the z-statistics from the interim analysis
  dAUC.1 = BDG.1$Ustat$AUCAminusAUCB[3]; z.1 = dAUC.1/sd.1

  ## conditional type I error and conditional power given the estimates in interim analysis for right sided test
  CE.inter = 1-pnorm(z.alpha*nr.init/nr.2*sd.total/sd.cond2-z.1*(nr.1*Vr.1+nr.1^2*Vc.1+nr.1*nr.2*Vc.3)/(nr.2*(Vr.1+nr.1*Vc.1))*sd.1/sd.cond2)
  CP.inter = 1-pnorm(z.alpha*nr.init/nr.2*sd.total/sd.cond2-z.1*nr.init/nr.2*sd.1/sd.cond2)

  #### compute the adjusted conditional power over a matrix of grid values of (Nr.1,...Nr.max) X (Nc.10,...Nc.0.max) ####
  #### and figure the maximum and mimimum conditional power based on Formula (xx) of the reference paper
  #### note that the adapted diseased case size is proportional to the adapted non-diseased case size (the ratio equals samples$R.Nc0toNc1)

  ## calculate the constant components in Formula (xx) of the reference paper
  CC1 = z.1*sd.1/(Vr.1+nr.1*Vc.1)
  CC2 = z.1*(nr.1*Vr.1+nr.1^2*Vc.1+nr.1*nr.2*Vc.3)/(nr.2*(Vr.1+nr.1*Vc.1))*sd.1/sd.cond2-z.alpha*nr.init/nr.2*sd.total/sd.cond2

  ## compute the adjusted conditional power matrix and figure out the maximum and minimum
  nr.grid = (nr.1+1):nr.max; nc.grid = nc.10:nc.0.max
  Otherparas =  list(nr.1 = nr.1, nc.10 = nc.10, nc.11 = nc.11,
                     R.Nc0toNc1 = samples$R.Nc0toNc1, Vr.1 = Vr.1, Vc.1 = Vc.1,
                     dAUC.moments.1 = dAUC.moments.1, CC1 = CC1, CC2 = CC2)
  CPkernel.adap.grid = sapply(nc.grid, function(xx){adaptiveII.CondPowerKernel(Otherparas, nr.grid, xx)$CPkernel.adap})
  CP.adap.grid = pnorm(CC1*CPkernel.adap.grid+CC2);
  CP.adap.low = min(CP.adap.grid[!is.na(CP.adap.grid)]); CP.adap.up = max(CP.adap.grid[!is.na(CP.adap.grid)])
  N.sim.na = sum(is.na(CP.adap.grid)) # number of grid items with NA adjusted conditional power
  CP.adap.grid.mod = CP.adap.grid; CP.adap.grid.mod[is.na(CP.adap.grid)] = 999 # fill grid items with NA adjusted conditional power using the value 999

  ## adjust reader and cases sizes
  if(z.1<=0){
    nr.adap = nr.init;
    nc.20.adap = nc.10
  }else{
    if(1-beta0<CP.adap.low){
      nr.adap = nr.1 + 1
      nc.20.adap = nc.10
    }else if(1-beta0>CP.adap.up){
      nr.adap = nr.max
      nc.20.adap = nc.0.max
    }else{
      targetpoint = which(abs(CP.adap.grid.mod-(1-beta0)) == min(abs(CP.adap.grid.mod-(1-beta0))), arr.ind = TRUE)[1,]
      nr.adap = nr.1 + 1 + targetpoint[1] - 1
      nc.20.adap = nc.10 + targetpoint[2] - 1
    }
  }

  ## adjust the critical value, and compute the adapted conditonal power and variance components
  tempt = adaptiveII.CondPowerKernel(Otherparas, nr.adap,nc.20.adap)
  cv.adap = tempt$cv.adap; CP.adap = pnorm(CC1*tempt$CPkernel.adap+CC2);
  Vr.2.adap = tempt$vr2.adap; Vc.2.adap = tempt$vc2.adap;
  sd.cond2.adap = tempt$sdcond2.adap; sd.total.adap = tempt$sdtotal.adap; rm(tempt)

  ## test decisions based on normal/t distribution in the interim analysis
  reject.inter.z = as.numeric(z.1 > z.alpha);
  reject.inter.t = as.numeric(BDG.1$Ustat$rejectBDG[3] && z.1 > 0)

  ## return the outputs
  result.interim1 = data.frame(nr.adap, nc.20.adap, cv.adap)
  result.interim2 = data.frame(Vr.1, Vc.1, Vr.2, Vc.2, Vr.2.adap, Vc.2.adap,
                               dAUC.1, sd.1, sd.cond2, sd.total, z.1,
                               CE.inter, CP.inter, CC1, CC2,
                               CP.adap.low, CP.adap.up, CP.adap, N.sim.na, sd.cond2.adap, sd.total.adap,
                               reject.inter.z, reject.inter.t)
  names(result.interim2) = c( "Vr.1", "Vc.1", "Vr.2", "Vc.2", "Vr.2.adap", "Vc.2.adap",
                              "dAUC.inter", "sd.inter", "sd.cond2", "sd.total", "z.inter",
                              "CE.inter", "CP.inter", "CC1", "CC2",
                              "CP.adap.low", "CP.adap.up", "CP.adap", "N.sim.na", "sd.cond2.adap", "sd.total.adap",
                              "reject.inter.z", "reject.inter.t")
  result.interim = data.frame(result.interim1, dAUC.moments.1, result.interim2)
  return(result.interim)
}

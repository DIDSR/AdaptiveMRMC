#################### Auxiliary Functions for Both Adaptive Methods I & II ####################

#' Calculate the vector (c1, c2, c3, c4)
#'
#' Calculate the vector (c1, c2, c3, c4) based on the numbers of non-diseased and disease cases
#'
#' @param nc0 [num] number of non-diseased cases
#' @param nc1 [num] number of diseasd cases
#' @return [num] vector of (c1, c2, c3, c4)
#' @export
fun.vector.c <- function(nc0, nc1){
  return(c(1/(nc0*nc1),(nc0-1)/(nc0*nc1),(nc1-1)/(nc0*nc1),(nc0-1)*(nc1-1)/(nc0*nc1)))
}

#' Calculate Vr and Vc
#'
#' Calculate Vr and Vc based on the estimated U-statistic moments and the vector (c1, c2, c3, c4)
#'
#' @param moments  [num] vector of the eight estimated U-statistics moments
#' @param vector.c [num] vector of (c1, c2, c3, c4)
#' @return [list] of Vr and Vc
#' @export
fun.VrVc <- function(moments, vector.c){
  if(length(vector.c) != 4 || length(moments) != 8) stop("v.c and moments must have length of 4 and 8 respectively")
  Vr = sum(vector.c*moments[1:4]) - sum(vector.c*moments[5:8])
  Vc = sum(vector.c*moments[5:8]) - moments[8]; Vc = as.numeric(Vc)
  return(list(Vr=Vr,Vc=Vc))
}

#################### Auxiliary Functions for Adaptive Methods II only ####################

#' Conditioanl variance and total variance w.r.t. adaptive method II
#'
#' Calculate the conditional variance of the second part reader-average AUC given the observed first part reader-average AUC
#' and the total variance where the first and second parts cases sizes are different
#' for adaptive method II
#'
#' @param nr1 [num] number of readers before the interim analysis
#' @param nr2 [num] number of readers after the interim analysis
#' @param vr1 [num] Vr based on the numbers of cases before the interim analysis
#' @param vc1 [num] Vc based on the numbers of cases before the interim analysis
#' @param vr2 [num] Vr based on the numbers of cases after the interim analysis
#' @param vc2 [num] Vc based on the numbers of cases after the interim analysis
#' @param vc3 [num] vc3 equals vc2
#' @return [list] of conditiona variance and the total variance
#' @export
adaptiveII.condvarR2.varR <- function(nr1, nr2, vr1, vc1, vr2, vc2, vc3){
  nume.condR2 = vr1*vr2+nr1*vr2*vc1+nr2*vr1*vc2+nr1*nr2*vc1*vc2-nr1*nr2*vc3^2
  deno.condR2 = nr2*(vr1+nr1*vc1)
  v.condR2 = nume.condR2/deno.condR2
  v.R = 1/(nr1+nr2)^2*(nr1*vr1+nr2*vr2+nr1^2*vc1+nr2^2*vc2+2*nr1*nr2*vc3)
  return(list(v.condR2=v.condR2,v.R=v.R))
}

#' Kernel component of calculating the adapted conditional power w.r.t. adaptive method II
#'
#' Calculate the kernel component of the adapted conditional power based on:
#' 1. a possible combination of adapted reader and cases sizes for the second part of the study,
#' 2. reader/cases sizes for the first part of the study,
#' 3. the estimates of the population parameters from the first part of the study.
#' And compute the corresponding adapted critical value.
#' Note that the adapted diseased case size is proportional to the adapted non-diseased case size (the ratio equals samples$R.Nc0toNc1)
#' Please refer to xxx paper for more detail.
#'
#' @param nr.adap        [num]  a possible adapted total reader sizes
#' @param nc20.adap      [num]  a possible adapted non-diseased case sizes
#' @param otherparas     [list] of the following:
#'
#'                  $nr.1:           [num]  number of readers for the first part of study (before the interim analysis)
#'
#'                  $nc.10:          [num]  number of non-diseased cases for the first part of the study (before the interim analysis)
#'
#'                  $nc.11:          [num]  number of diseased cases for the first part of the study (before the interim analysis)
#'
#'                  $ratio.nc0tonc1: [num]  the ratio of the number of non-diseased cases to that of diseased cases
#'
#'                  $Vr.1            [num]  estimated variance component w.r.t. the first part of study
#'
#'                  $Vc.1            [num]  estimated variance component w.r.t. the first part of study
#'
#'                  $dAUC.moments.1  [list] estimated U-statistics moments
#'
#'                  $CC1             [num]  constant component for calculating the kernel power component in Formula (xx) of the reference paper
#'
#'                  $CC2             [num]  constant component for calculating the kernel power component in Formula (xx) of the reference paper
#'
#' @return [list]:
#'
#'           $CPkernel.adap:  [num] the kernel component of the adapted conditional power
#'
#'           $cv.adap:        [num] the corresponding adapted critical value
#'
#'           $...: other potentially necessary outputs
#' @export

adaptiveII.CondPowerKernel = function(otherparas, nr.adap, nc20.adap){
  nr1 = otherparas$nr.1; nc10 = otherparas$nc.10; nc11 = otherparas$nc.11; ratio.nc0tonc1 = otherparas$R.Nc0toNc1
  vr1 = otherparas$Vr.1; vc1 = otherparas$Vc.1; daucmoments1 = otherparas$dAUC.moments.1;
  cc1 = otherparas$CC1; cc2 = otherparas$CC2

  nc21.adap = round(nc20.adap/ratio.nc0tonc1); nr2.adap = nr.adap - nr1

  vectorc2.adap = fun.vector.c(nc20.adap,nc21.adap)
  vr2.adap = fun.VrVc(daucmoments1, vectorc2.adap)$Vr;
  vc2.adap = fun.VrVc(daucmoments1, vectorc2.adap)$Vc;
  #vc3.adap = nc10*nc11/(nc20.adap*nc21.adap)*vc1; ## This derivation of Vc3.adap is wrong
  vc3.adap = vc2.adap
  sdcond2.adap = sqrt(adaptiveII.condvarR2.varR(nr1, nr2.adap, vr1, vc1, vr2.adap, vc2.adap, vc3.adap)$v.condR2)
  sdtotal.adap = sqrt(adaptiveII.condvarR2.varR(nr1, nr2.adap, vr1, vc1, vr2.adap, vc2.adap, vc3.adap)$v.R)

  part1 = nr.adap/nr2.adap*sdtotal.adap/sdcond2.adap
  part2 = (nr1*vr1+nr1^2*vc1+nr1*nr2.adap*vc3.adap)/(nr2.adap*sdcond2.adap)
  cv.adap = 1/part1*(cc1*part2-cc2)

  CPkernel.adap  = (vr1+nr1*vc1-nr1*vc3.adap)/sdcond2.adap

  return(list(CPkernel.adap = CPkernel.adap, cv.adap = cv.adap, cv.adap.part1 = part1, cv.adap.part2 = part2,
              vr2.adap = vr2.adap, vc2.adap = vc2.adap,
              sdcond2.adap = sdcond2.adap,sdtotal.adap = sdtotal.adap))
}

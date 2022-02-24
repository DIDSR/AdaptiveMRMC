########################################################################################################
#' Generate MRMC simulation scores using constrained unequal-variance Roe and Metz model
#'
#' The MRMC simulation scores are generated using the linear mixed effect model (Roe and Metz model) and are fully paired across two modalities (fully-cross design).
#'
#' @param nr [num] number of readers
#' @param nc0 [num] number of signal-absent (non-diseased) cases
#' @param nc1 [num] number of signal-present (diseased) cases
#' @param paras [list] of parameters:
#'
#'               $mu1: [num] grand average
#'
#'               $tauA1: [num] modality effect (for the signal-present class)
#'
#'               $Vr0, $Vtaur0, $Vc0, $Vtauc0, $Vrc0, $Veps0, $b: [num] variance component parameters (see references)
#'
#' @return  dFrame.imrmc [data.frame] with (nc0+nc1)+nr*(nc0+nc1)*2 rows and 4 variables including:
#'
#'              $readerID: [Factor] w/ nr+1 levels "-1", "reader1", "reader2", ...
#'
#'              $caseID: [Factor] w/ nc0+nc1 levels "negCase1", "negCase2", ..., posCase1", "posCase2", ...
#'
#'              $modalityID: [Factor] w/ 3 levels "truth", "modility A", "modility B"
#'
#'              $score: [num] reader score
#' @return S10: [num] [nc0 X nr matrix] score w.r.t. modality A, signal absent
#' @return S11: [num] [nc1 X nr matrix] score w.r.t. modality A, signal present
#' @return S20: [num] [nc0 X nr matrix] score w.r.t. modality B, signal absent
#' @return S21: [num] [nc1 X nr matrix] score w.r.t. modality B, signal present
#'
#' @export
#' @references Roe and Metz, Acad Radiol 1997, 298-303; Hillis, Acad Radiol 2012, 19:1518.
## @example
## @author  Zhipeng Huang, 04/03/2018


mrmcRMscoresFC <- function(nr, nc0, nc1, paras){

  # tau0 = (modality A with normal, modality B with normal)
  # tau1 = (modality A with disease, modality B with disease)
  # re-denotation the parameters
  mu1 = paras$mu1; tau1 = c(paras$tauA1,0); mu0 = 0; tau0 = c(0, 0)
  vr = paras$Vr0; vtr = paras$Vtaur0
  vc = paras$Vc0; vtc = paras$Vtauc0; vrc = paras$Vrc0; ve = paras$Veps0
  b  = paras$b

  # Assign caseIDs and readerIDs
  negCaseIDs <- factor(paste("negCase", 1:nc0, sep = ""))
  posCaseIDs <- factor(paste("posCase", 1:nc1, sep = ""))
  readerIDs <- factor(paste("reader", 1:nr, sep = ""))

  # Create data frame of truth
  dFrame.truth <- data.frame(
    readerID = rep("-1", nc0 + nc1),
    caseID = c(as.character(negCaseIDs), as.character(posCaseIDs)),
    modalityID = rep("truth", nc0 + nc1),
    score = c(rep(0, nc0), rep(1, nc1))
  )

  # effects that are shared by two modalities (but independent of truth states)
  R0 <- rnorm(nr) * sqrt(vr)
  R1 <- rnorm(nr) * sqrt(vr)
  C0 <- rnorm(nc0) * sqrt(vc)
  C1 <- rnorm(nc1) * sqrt(vc)/b
  RC0 <- matrix(rnorm(nr * nc0) * sqrt(vrc), nrow = nc0)
  RC1 <- matrix(rnorm(nr * nc1) * sqrt(vrc)/b, nrow = nc1)

  # modality A, signal absent
  tauR10 <- rnorm(nr) * sqrt(vtr)
  tauC10 <- rnorm(nc0) * sqrt(vtc)
  eps10  <- matrix(rnorm(nc0 * nr) * sqrt(ve), nrow = nc0)
  S10 <- mu0 + tau0[1] + t(replicate(nc0, R0 + tauR10)) + replicate(nr, C0 + tauC10) + RC0 + eps10
  dFrame.modAneg <- data.frame(readerID = rep(readerIDs, rep(nc0, nr)),
                               caseID = rep(negCaseIDs,nr),
                               modalityID = rep("modalityA", nc0 * nr),
                               score = as.vector(S10))
  # modality A, signal present
  tauR11 <- rnorm(nr) * sqrt(vtr)
  tauC11 <- rnorm(nc1) * sqrt(vtc)/b
  eps11  <- matrix(rnorm(nc1 * nr) * sqrt(ve)/b, nrow = nc1)
  S11 <- mu1 + tau1[1] + t(replicate(nc1, R1 + tauR11)) + replicate(nr, C1 + tauC11) + RC1 + eps11
  dFrame.modApos <- data.frame(readerID = rep(readerIDs, rep(nc1, nr)),
                               caseID = rep(posCaseIDs,nr),
                               modalityID = rep("modalityA", nc1 * nr),
                               score = as.vector(S11))

  # modality B, signal absent
  tauR20 <- rnorm(nr) * sqrt(vtr)
  tauC20 <- rnorm(nc0) * sqrt(vtc)
  eps20  <- matrix(rnorm(nc0 * nr) * sqrt(ve), nrow = nc0)
  S20 <- mu0 + tau0[2] + t(replicate(nc0, R0 + tauR20)) + replicate(nr, C0 + tauC20) + RC0 + eps20
  dFrame.modBneg <- data.frame(readerID = rep(readerIDs, rep(nc0, nr)),
                               caseID = rep(negCaseIDs,nr),
                               modalityID = rep("modalityB", nc0 * nr),
                               score = as.vector(S20))

  # modality B, signal present
  tauR21 <- rnorm(nr) * sqrt(vtr)
  tauC21 <- rnorm(nc1) * sqrt(vtc)/b
  eps21  <- matrix(rnorm(nc1 * nr) * sqrt(ve)/b, nrow = nc1)
  S21 <- mu1 + tau1[2] + t(replicate(nc1, R1 + tauR21)) + replicate(nr, C1 + tauC21) + RC1 + eps21
  dFrame.modBpos <- data.frame(readerID = rep(readerIDs, rep(nc1, nr)),
                               caseID = rep(posCaseIDs,nr),
                               modalityID = rep("modalityB", nc1 * nr),
                               score = as.vector(S21))

  dFrame.imrmc <- rbind(
    dFrame.truth,
    dFrame.modAneg,
    dFrame.modBneg,
    dFrame.modApos,
    dFrame.modBpos
  )
  list(S10 = S10, S11 = S11, S20 = S20, S21 = S21, dFrame.imrmc = dFrame.imrmc)
}




#' Analytical variance parameters and power in the AUC space give Roe-Metz simulation parameters
#'
#' Compute the analytical ariance parameters and power in the AUC space give Roe-Metz simulation parameters. Original version was created by Brandon Gallas and Weijie Chen. The original arguments for population parameters are combined as a list argument paras.
#' The function requires the jiRometz_SPsim2.jar file. Make sure the work directory has this jar file.
#'
#' @param Nr [num] or [vector] (vector of) the number of readers
#' @param N0 [num] or [vector] (vector of) the number of signal-absent (non-diseased) cases
#' @param N1 [num] or [vector] (vector of) the number of signal-present (diseased) cases
#' @param paras [list] of parameters:
#'
#'               $mu1: [num] grand average
#'
#'               $tauA1: [num] modality effect (for the signal-present class)
#'
#'               $Vr0, $Vtaur0, $Vc0, $Vtauc0, $Vrc0, $Veps0, $b: [num] variance component parameters (see references)
#'
#' @return  [list] of 3:
#'
#'           $var.da:  [num] analytical variance
#'
#'           $moments: [num] [1 X 8 vector] analytical U-statistics moments
#'
#'           $pow:     [num] analytical power
#' @export
#' @references Roe and Metz, Acad Radiol 1997, 298-303; Hillis, Acad Radiol 2012, 19:1518.

analytical.var.AUC <- function(Nr, N0, N1, paras){
  mu1 = paras$mu1; tau1 = c(paras$tauA1,0); vr = paras$Vr0; vtr = paras$Vtaur0
  vc = paras$Vc0; vtc = paras$Vtauc0; vrc = paras$Vrc0; ve = paras$Veps0; b  = paras$b

  nG <- length(N0)
  if(nG != length(Nr) | nG != length(N1)) stop("The length of Nr, N0, and N1 must be equal.")

  input.file <- paste("_tmp_input_", round(abs(rnorm(1))*1000), ".irm", sep = "")
  strA0 <- c(paste("AR0:", vtr), paste("AC0:", vtc), paste("ARC0:", ve))
  strA1 <- c(paste("AR1:", vtr), paste("AC1:", vtc/b^2), paste("ARC1:", ve/b^2))
  strB0 <- c(paste("BR0:", vtr), paste("BC0:", vtc), paste("BRC0:", ve))
  strB1 <- c(paste("BR1:", vtr), paste("BC1:", vtc/b^2), paste("BRC1:", ve/b^2))
  str0 <- c(paste("R0:", vr), paste("C0:", vc), paste("RC0:", vrc))
  str1 <- c(paste("R1:", vr), paste("C1:", vc/b^2), paste("RC1:", vrc/b^2))
  stru <- c(paste("uA:", mu1+tau1[1]), paste("uB:", mu1+tau1[2]), paste("n0:", sum(N0)),
            paste("n1:", sum(N1)), paste("nr:", sum(Nr)))
  strex <- c("Study Design",
             paste("# of Split-Plot Groups:", nG),
             "Paired Readers: Yes", "Paired Normal: Yes",
             "Paired Disease: Yes", "Number of Experiments: 100",
             "Seed for RNG: 1", "Random Stream: 5000", "MLE analysis: NO")
  writeLines(c(strA0, strA1, strB0, strB1, str0, str1, stru, strex), con = input.file)
  if(Sys.info()['sysname']=="Windows")
    system(paste("java -jar iRometz_SPsim2.jar ", getwd(), "/", input.file, sep=""), show.output.on.console = FALSE)
  if(Sys.info()['sysname']=="Linux")
    system(paste("java -jar iRometz_SPsim2.jar ", getwd(), "/", input.file, sep=""), ignore.stdout = TRUE)

  system(paste("java -jar iRometz_SPsim2.jar ", getwd(), "/", input.file, sep=""), ignore.stdout = T)
  m <- read.csv(paste(substr(input.file, 1, nchar(input.file)-4), "NumericalMoment.csv", sep = ""), header = T)
  #print(m)
  unlink(input.file)
  unlink(paste(substr(input.file, 1, nchar(input.file)-4), "NumericalMoment.csv", sep = ""))
  md <- data.matrix(m[1, 2:9] + m[2, 2:9] - 2 * m[3, 2:9])
  ci <- matrix(0, nrow = nG, ncol = 4)
  Vr <- rep(0, nG); Vc <- rep(0, nG)
  var.da <- 0
  for(g in 1:nG){
    ci[g,] <- c(1.0/(N0[g]*N1[g]), (N0[g] - 1.0)/(N0[g]*N1[g]),
                   (N1[g] - 1.0)/(N0[g]*N1[g]), (N1[g] - 1.0)*(N0[g] - 1.0)/(N0[g]*N1[g]))
    Vr[g] <- ci[g,1] * (md[1]-md[5]) + ci[g,2] * (md[2]-md[6]) + ci[g,3] * (md[3]-md[7]) + ci[g,4] * (md[4]-md[8])
    Vc[g] <- ci[g,1] * md[5] + ci[g,2] * md[6] + ci[g,3] * md[7] - (1-ci[g,4]) * md[8]
    var.da <- var.da + Nr[g]*Vr[g] + Nr[g]^2*Vc[g]
  }
  var.da <- var.da/sum(Nr)^2
  var.total <- 2 * (vr + vtr) + (vc + vtc + vrc + ve) * (1 + 1/b)
  dauc <- pnorm((mu1+tau1[1])/sqrt(var.total)) - pnorm((mu1+tau1[2])/sqrt(var.total))
  pow <- pnorm(dauc/sqrt(var.da) - qnorm(0.975)) + pnorm(-dauc/sqrt(var.da) - qnorm(0.975))
  list(var.da = var.da, moments = md, pow = pow)
}

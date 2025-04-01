###########################################################################################
# Internal function: fit multicomponent FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   nback: number of FMM components to be fitted.
#   timePoints: one single period time points.
#   maxiter: maximum number of iterations for the backfitting algorithm.
#   stopFunction: function to check the criterion convergence for the backfitting algorithm.
#   lengthAlphaGrid, lengthOmegaGrid: precision of the grid of alpha and omega parameters.
#   alphaGrid: grid of alpha  parameters.
#   omegaMin, omegaMax: min and max values for omega, ang grid of omega parameters.
#   (DEPRECATED) numReps: number of times the alpha-omega grid search is repeated.
#   showProgress: TRUE to display a progress indicator on the console.
#   (DEPRECATED) usedApply: paralellized version of apply for grid search
#   gridList: list that contains precalculations to make grid search
#             calculations lighter
# Returns an object of class FMM.
###########################################################################################
fitFMM_back<-function(vData, nback, timePoints = seqTimes(length(vData)),
                      maxiter = nback, stopFunction = alwaysFalse,
                      lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                      alphaGrid = seq(0, 2*pi, length.out = lengthAlphaGrid),
                      omegaMin = 0.0001, omegaMax = 0.999,
                      omegaGrid = exp(seq(log(omegaMin), log(omegaMax),
                                          length.out = lengthOmegaGrid+1))[1:lengthOmegaGrid],
                      numReps = 1, showProgress = TRUE, usedApply = NA,
                      gridList = precalculateBase(alphaGrid = alphaGrid, omegaGrid = omegaGrid,
                                                  timePoints = timePoints)){

  if(!is.na(usedApply)){
    warning("Argument 'usedApply' is deprecated.")
  }

  if(numReps>1){
    warning("Argument 'numReps' is deprecated.")
  }

  nObs <- length(vData)

  if(showProgress){
    totalMarks <- 50
    partialMarkLength <- 2
    progressHeader<-paste(c("|", rep("-", totalMarks), "|\n|"), collapse = "")
    cat(progressHeader)
    completedPercentage <- 0.00001
    previousPercentage <- completedPercentage
  }

  alpha <- rep(NA, nback)
  omega <- rep(NA, nback)
  beta <- rep(NA, nback)
  A <- rep(NA, nback)

  # Object initialization.
  prevFittedFMMvalues <- NULL
  stopCriteria <- "Stopped by reaching maximum iterations ("
  nIter <- 1

  blaschkeProduct <- rep(1,nObs) # Neutral for product
  backFittingData <- vData

  # First Iteration
  for(j in 1:nback){
    # component j fitting using fitFMM_unit function
    fittedFMM <- fitFMM_unit(backFittingData, timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                             lengthOmegaGrid = lengthOmegaGrid, alphaGrid = alphaGrid, omegaMin = omegaMin,
                             omegaMax = omegaMax, omegaGrid = omegaGrid, gridList = gridList)
    alpha[j] <- fittedFMM@alpha; omega[j] <- fittedFMM@omega
    aj <- (1-fittedFMM@omega)/(1+fittedFMM@omega) * exp(1i*(fittedFMM@alpha+pi))

    backFittingData <- Re(gsignal::hilbert(backFittingData - getFittedValues(fittedFMM)) / mobius(aj, timePoints))
    blaschkeProduct <- blaschkeProduct * mobius(aj, timePoints)

    # showProgress
    if(showProgress){
      completedPercentage <- completedPercentage + 100/(nback*maxiter)
      if(ceiling(previousPercentage) < floor(completedPercentage)){
        progressDone <- paste(rep("=",sum((seq(ceiling(previousPercentage), floor(completedPercentage), by = 1)
                                           %% partialMarkLength == 0))), collapse = "")
        cat(progressDone)
        previousPercentage <- completedPercentage
      }
    }
  }
  DM <- calculateDesignMatrix(alpha = alpha, omega = omega, timePoints = timePoints)
  prevFittedFMMvalues <- as.numeric(DM %*% qr.solve(DM, vData))

  # Backfitting algorithm: iteration (only if maxiter > 1)
  if(maxiter>1){
    for(i in 2:maxiter){
      # Backfitting algorithm: component
      for(j in 1:nback){
        # Standard residuals.
        DM <- calculateDesignMatrix(alpha = alpha[-j], omega = omega[-j], timePoints = timePoints)
        fittedFMMvalues <- as.numeric(DM %*% qr.solve(DM, vData))
        stdResiduals <- vData - fittedFMMvalues

        # Transformed standard residuals (Blaschke product).
        aj <- (1-omega[j])/(1+omega[j]) * exp(1i*(alpha[j]+pi))
        blaschkeProduct <- blaschkeProduct / mobius(aj, timePoints)
        backFittingData <- Re(gsignal::hilbert(stdResiduals)/blaschkeProduct)

        # component j fitting using fitFMM_unit function
        fittedFMM <- fitFMM_unit(backFittingData, timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                                                  lengthOmegaGrid = lengthOmegaGrid, alphaGrid = alphaGrid, omegaMin = omegaMin,
                                                  omegaMax = omegaMax, omegaGrid = omegaGrid, gridList = gridList)
        # Update Blaschke product
        alpha[j] <- fittedFMM@alpha; omega[j] <- fittedFMM@omega
        aj <- (1-omega[j])/(1+omega[j]) * exp(1i*(alpha[j]+pi))
        blaschkeProduct <- blaschkeProduct * mobius(aj, timePoints)

        # showProgress
        if(showProgress){
          completedPercentage <- completedPercentage + 100/(nback*maxiter)
          if(ceiling(previousPercentage) < floor(completedPercentage)){
            progressDone <- paste(rep("=",sum((seq(ceiling(previousPercentage), floor(completedPercentage), by = 1)
                                               %% partialMarkLength == 0))), collapse = "")
            cat(progressDone)
            previousPercentage <- completedPercentage
          }
        }
      }

      # Check stop criterion
      # Fitted values as sum of all components
      DM <- calculateDesignMatrix(alpha = alpha, omega = omega, timePoints = timePoints)
      fittedFMMvalues <- as.numeric(DM %*% qr.solve(DM, vData))

      if(!is.null(prevFittedFMMvalues)){
        if(PV(vData, prevFittedFMMvalues) > PV(vData, fittedFMMvalues)){
          fittedFMMvalues <- prevFittedFMMvalues
          stopCriteria <- "Stopped by reaching maximum R2 ("
          break
        }
        if(stopFunction(vData, fittedFMMvalues, prevFittedFMMvalues)){
          stopCriteria <- "Stopped by the stopFunction ("
          break
        }
      }
      prevFittedFMMvalues <- fittedFMMvalues
    }
    nIter <- i

    # showProgress
    if(showProgress){
      if(completedPercentage < 100){
        completedPercentage <- 100
        nMarks <- ifelse(ceiling(previousPercentage) < floor(completedPercentage),
                         sum((seq(ceiling(previousPercentage),
                                  floor(completedPercentage), by = 1) %% partialMarkLength == 0)), 0)
        if (nMarks > 0) {
          cat(paste(rep("=",nMarks), collapse = ""))
          previousPercentage <- completedPercentage
        }
      }
      cat("|\n", paste(stopCriteria, nIter, sep = ""),"iteration(s))","\n")
    }
  }

  # linear parameters are recalculated by linear regression; DM is the design matrix
  DM <- calculateDesignMatrix(alpha = alpha, omega = omega, timePoints = timePoints)
  coefs <- qr.solve(DM, vData)
  fittedFMMvalues <- as.numeric(DM %*% coefs)

  # A, betas and M estimates
  M <- coefs[1]
  for (i in 1:nback) {
    A[i] = sqrt(coefs[2*i]^2+coefs[2*i+1]^2)
    beta[i] = atan2(-coefs[2*i+1], coefs[2*i]) %% (2*pi)
  }

  # Residual sum of squares
  SSE <- sum((vData - fittedFMMvalues)^2)

  # Reorder waves by explained variability
  iPV <- sapply(1:nback, function(x){PV(vData, predict(lm(vData ~ DM[,c(2*x, 2*x+1)])))})

  iPVF <- iPV*0
  waveOrder<-which.max(iPV)
  while(length(waveOrder)<nback){
    iPV <- sapply((1:nback)[-waveOrder], function(x){PV(vData, predict(lm(vData ~ DM[,c(2*waveOrder, 2*waveOrder+1, 2*x, 2*x+1)])))})
    waveOrder<-c(waveOrder,(1:nback)[-waveOrder][which.max(iPV)])
  }

  A <- A[waveOrder]; alpha <- alpha[waveOrder]; beta <- beta[waveOrder]; omega <- omega[waveOrder]

  # Returns an object of class FMM.
  return(FMM(
    M = M,
    A = A,
    alpha = alpha,
    beta = beta,
    omega = omega,
    timePoints = timePoints,
    summarizedData = vData,
    fittedValues= fittedFMMvalues,
    SSE = SSE,
    R2 = PVj(vData, timePoints, alpha, beta, omega),
    nIter = nIter
  ))
}

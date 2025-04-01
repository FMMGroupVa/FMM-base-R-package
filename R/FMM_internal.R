################################################################################
# Auxiliary internal functions
# Functions:
#   step1FMM:        M, A and beta initial parameter estimations.
#   step2FMM:        second step of FMM fitting process.
#   refineFMM:       fitFMM from a previous objectFMM.
#   PV:              percentage of variability explained.
#   PVj:             percentage of variability explained by each component of
#                    FMM model.
#   seqTimes:        to build a sequence of equally time points spaced in range
#                    [0,2*pi].
#   calculateCosPhi: to calculate components' cos(phi(t)).
################################################################################

################################################################################
# Internal function: Returns the RSS for an estimation given (alpha,omega).
# Arguments:
#   - optBase: A list containing precalculated values to avoid redundant operations:
#       1. "base": inv(X'X)X', where X = [1, cos(tStar), sin(tStar)]
#       2. "alpha"
#       3. "omega"
#       4. "cost": cos(tStar)
#       5. "sint": sin(tStar)
#
#   - vData: Numeric vector representing the data to be fitted with the FMM model.
#
# Returns:
#   A numeric vector of length 3 with the following elements: alpha, omega, RSS
################################################################################

step1FMM <- function(optBase, vData) {
  # Linear parameters estimation and RSS
  pars <- optBase[["base"]] %*% vData
  residualSS <- sum((vData - pars[1] - pars[2]*optBase[["cost"]] - pars[3]*optBase[["sint"]])^2)
  return(c(optBase[["alpha"]], optBase[["omega"]], residualSS))
}

################################################################################
# Internal function: Profile likelihood (alpha, omega)
# Arguments:
#   parameters: alpha, omega initial parameter estimations
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   omegaMin: min value for omega (must be >0).
#   omegaMax: max value for omega (must be <1).
################################################################################

step2FMM <- function(parameters, vData, timePoints, omegaMin, omegaMax){
  if(parameters[2] >= omegaMin  &  parameters[2] <= omegaMax){
    # Linear parameters estimation
    nonlinearMob = 2*atan(parameters[2]*tan((timePoints-parameters[1])/2))
    DM <- cbind(rep(1, length(timePoints)), cos(nonlinearMob), sin(nonlinearMob))
    pars <- stats::.lm.fit(DM, vData)$coefficients
    # Return RSS
    return(sum((vData - pars[1] - pars[2]*cos(nonlinearMob) - pars[3]*sin(nonlinearMob))^2))
  }else{
    return(Inf)
  }
}

################################################################################
# Internal function: to calculate the percentage of variability explained by
#   the FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   pred: fitted values.
################################################################################
PV <- function(vData, pred){
  return(1 - sum((vData - pred)^2)/sum((vData - mean(vData))^2))
}

################################################################################
# Internal function: to calculate the percentage of variability explained by
#   each component of FMM model
# Arguments:
#   vData: data to be fitted an FMM model.
#   timePoints: one single period time points.
#   alpha, beta, omega: vectors of corresponding parameter estimates.
################################################################################
PVj <- function(vData, timePoints, alpha, beta, omega){
  # Fitted values of each wave
  DM <- calculateDesignMatrix(alpha = alpha, omega = omega, timePoints = timePoints)

  # The percentage of variability explained up to wave i is determined
  cumulativePV <- sapply(1:length(alpha), function(x){
    PV(vData, predict(lm(vData ~ DM[,1:(2*x+1)]-1)))
  })

  # individual percentage of variability is the part that adds to the whole
  return(c(cumulativePV[1], diff(cumulativePV)))
}

################################################################################
# Internal function: to build a sequence of equally time points spaced
#                    in range [0,2*pi).
# Arguments:
#   nObs: secuence length.
################################################################################
seqTimes <- function(nObs){
  return(seq(0, 2*pi, length.out = nObs+1)[1:nObs])
}

################################################################################
# Internal function: to calculate a single components' cos(phi(t)).
# Arguments:
#   alpha, beta, omega: parameters.
#   timePoints: time points in which the FMM model is computed.
# Returns a matrix of each component's cos(phi(t)) as columns.
################################################################################
calculateSingleCosPhi <- function(alpha, beta, omega, timePoints){
  return(cos(beta + 2*atan(omega*tan((timePoints - alpha)/2))))
}

################################################################################
# Internal function: to calculate components' cos(phi(t)).
# Arguments:
#   alpha, beta, omega: parameters.
#   timePoints: time points in which the FMM model is computed.
# Returns a matrix of each component's cos(phi(t)) as columns.
################################################################################
calculateCosPhi <- function(alpha, beta, omega, timePoints){
  calculateSingleCosPhi <- function(alpha, beta, omega, timePoints){
    return(cos(beta + 2*atan(omega*tan((timePoints - alpha)/2))))
  }
  return(sapply(1:length(alpha), FUN = function(x){
    calculateSingleCosPhi(alpha = alpha[x], beta = beta[x], omega = omega[x], timePoints = timePoints)
  }))
  # return(mapply(FUN = calculateSingleCosPhi, alpha = alpha, beta = beta, omega = omega, timePoints=timePoints))
}

################################################################################
# Internal function: to calculate design matrix with fixed alphas and omegas
# Arguments:
#   alpha, omega: parameters.
#   timePoints: time points in which the FMM model is computed.
# Returns the design matrix (fixed alphas and omegas)
################################################################################
calculateDesignMatrix <- function(alpha, omega, timePoints){
  design_matrix <- matrix(1, nrow = length(timePoints), ncol = 1)
  for (i in seq_along(alpha)) {
    cos_col <- calculateCosPhi(alpha[i], 0, omega[i], timePoints)
    sin_col <- calculateCosPhi(alpha[i], -pi/2, omega[i], timePoints) # sin(f) = cos(f-pi/2)
    design_matrix <- cbind(design_matrix, cos_col, sin_col)
  }
  return(design_matrix)
}

################################################################################
# Internal function: to precalculate inv(M'M)M' for M=[1, cos(t*), sin(t*)].
# Arguments:
#   alphagrid, omegaGrid: search grid.
#   timePoints: time points in which the FMM model is computed.
# Returns a list where each element is a list with elements:
#   base: inv(M'M)M',
#   alpha, omega,
#   cost: cos(tStar), sint: sin(tStar)
################################################################################

precalculateBase <- function(alphaGrid, omegaGrid, timePoints){
  # Expanded grid: each row contains a pair (alpha, omega)
  grid <- expand.grid(alphaGrid, omegaGrid)
  optBase <- apply(grid, 1, FUN = function(x){
    x <- as.numeric(x)
    DM <- calculateDesignMatrix(x[1], x[2], timePoints)
    mat <- t(DM) %*% DM

    # Check if matrix is invertible, too low omegas and extreme timepoints
    # configurations lead to computationally singular systems
    if (det(mat) < 10e-12){
      return(NULL)
    }else{
      return(list(base = solve(mat)%*%t(DM),
                  alpha = x[1], omega = x[2],
                  cost = cos(2*atan(x[2]*tan((timePoints-x[1])/2))),
                  sint = sin(2*atan(x[2]*tan((timePoints-x[1])/2)))))
    }
  }, simplify = FALSE)
  optBase <- Filter(Negate(is.null), optBase)
  return(optBase)
}

################################################################################
# Internal function: to calculate the Mobius transform in e^it with argument a
# Arguments:
#   a: parameter.
#   timePoints: time points in which the Mobius transform is computed.
# Returns the transform
################################################################################
mobius <- function(a, t){
  return((exp(1i*t)-a) / (1-Conj(a)*exp(1i*t)))
}


################################################################################
# Internal function: to check if the arguments passed to fitFMM are correct
# Arguments: See fitFMM function arguments
################################################################################

checkArguments <- function(vData, nPeriods, timePoints, nback, maxiter,
                           betaOmegaRestrictions, omegaMin, omegaMax,
                           lengthAlphaGrid, lengthOmegaGrid,
                           omegaGrid, numReps, parallelize){
  # Check grid values
  if(omegaMin<=0) stop("Incorrect arguments: omegaMin must be greater than 0")
  if(omegaMax>=1) stop("Incorrect arguments: omegaMax must be lower than 1")
  if(omegaMin>=omegaMax) stop("Incorrect arguments: omegaMin must be lower than omegaMax")
  if(lengthAlphaGrid <= 0) stop("Incorrect arguments: lengthAlphaGrid must be greater than 0")
  if(lengthOmegaGrid <= 0) stop("Incorrect arguments: lengthOmegaGrid must be greater than 0")
  if(!is.null(omegaGrid) & length(omegaGrid) == 0) stop("Incorrect arguments: omegaGrid is empty")
  if(!is.null(omegaGrid) & length(omegaGrid[omegaGrid>omegaMin & omegaGrid<omegaMax]) == 0)
    stop("Incorrect arguments: omegaGrid has no values in range [omegaMin, omegaMax]")

  # Check data input
  if(length(vData) < 5) stop("The minimum number of observations should be 5")
  if(nPeriods > 1 & length(vData) %% nPeriods != 0) stop("Data length is not a multiple of nPeriods")
  if(sd(vData, na.rm=TRUE) == 0) stop("Data are constant")

  # Check timePoints input
  if(!is.null(timePoints)){
    timePoints <- sort(timePoints)
    if(any(timePoints < 0) | any(timePoints > 2*pi)) stop("timePoints must be between 0 and 2*pi")
    if(nPeriods == 1 & length(timePoints) != length(vData))
      stop("timePoints must have the same length as one-period data")
    if(max(diff(c(utils::tail(timePoints,1), timePoints)) %% (2*pi)) > pi/4 & length(vData) > 10)
      warning("Detected large gaps between time points, which may cause unstable or nonsensical solutions.")
  }


  # Warn deprecated arguments
  if(length(unique(betaOmegaRestrictions)) == nback && numReps>1)
    warning("Argument 'numReps' is deprecated for unrestricted fittings. For a finer grid search, increase 'lengthAlphaGrid' and 'lengthOmegaGrid'.")
  if(length(unique(betaOmegaRestrictions)) == nback && parallelize == TRUE)
    warning("Argument 'parallelize' is deprecated for unrestricted fittings.")
}

################################################################################
# Internal function: to check if the solution is valid
# Arguments: object of class FMM, omegaMin and omegaMax
################################################################################
checkSolution <- function(fittedFMM, omegaMin, omegaMax){

  # Omegas
  if(any(getOmega(fittedFMM) <= omegaMin)) warning("Detected unusually small omegas. Check your model parameters and input data.")

  # Amplitudes
  if(any(getA(fittedFMM) < 0)) stop("Invalid solution: check function input parameters.")

  slack <- 0.1*diff(range(fittedFMM@data))
  if(any(getOmega(fittedFMM) <= 0.1 &
         2*getA(fittedFMM) > diff(range(fittedFMM@data)) + 4*sd(fittedFMM@data-fittedFMM@fittedValues) + 2*slack
         )){
    warning("Unusually large amplitudes detected. This may indicate issues with model parameters or input data. If omega values are too low or the dataset is sparse, consider increasing omegaMin or interpolating your data for better stability.")
  }
  # if(any(getOmega(fittedFMM) <= 0.1)){
  #
  #   residuals <- fittedFMM@data
  #   amplitudeFlag <- FALSE
  #   for (i in 1:length(fittedFMM@alpha)) {
  #     nonlinearMob = 2*atan(fittedFMM@omega[i]*tan((fittedFMM@timePoints-fittedFMM@alpha[i])/2))
  #     DM <- cbind(fittedFMM@timePoints*0+1, cos(nonlinearMob), sin(nonlinearMob))
  #     pars <- stats::.lm.fit(DM, residuals)$coefficients
  #     sigmaHat <- sd(residuals - pars[1] - pars[2]*cos(nonlinearMob) - pars[3]*sin(nonlinearMob))
  #     Ai <- fittedFMM@A[i]
  #     if((pars[1]+Ai > max(residuals)+2*sigmaHat+slack |
  #        pars[1]-Ai < min(residuals)-2*sigmaHat-slack) &
  #        fittedFMM@omega[i] <= 0.1) amplitudeFlag = TRUE
  #     residuals <- residuals - pars[1] - pars[2]*cos(nonlinearMob) - pars[3]*sin(nonlinearMob)
  #   }
  #   if(amplitudeFlag) warning("Unusually large amplitudes detected. This may indicate issues with model parameters or input data. If omega values are too low and the dataset is sparse, consider increasing omegaMin or interpolating your data for better stability.")
  # }
}

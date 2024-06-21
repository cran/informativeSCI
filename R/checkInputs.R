# Implementation of validation functions for the input values


###############################################################################
# Validation function for part of the input values
# of the informSCI-function, the inExactSCI -function, 
# the notInExactSCI-function and the explore_q-function.

checkPartOfInput1 <- function(gMCP, g, weights, mu_0, alpha){
  
  if(is.null(gMCP)  &&  is.null(g)){
    stop("Either gMCP or g and weights must be specified!")  
  }
  if(is.null(gMCP) && is.null(weights)){
    stop("Either gMCP or g and weights must be specified!")
  }
  if(!is.null(gMCP) && !is.null(g)){
    message("Only one parameter of either gMCP or g needs to be specified!
            The transition matrix of gMCP is used for the algorithmn.
            The matrix g is ignored!") 
  }
  if(!is.null(gMCP) && !is.null(weights)){
    message("Only one parameter of either gMCP or weights needs to be 
            specified! The weights of gMCP are used for the algorithmn. 
            The entered weights vector is ignored!") 
  }
  if(is.null(gMCP) && !(is.null(g))){ 
    dimOf_g <- dim(g)
    if(length(unique(dimOf_g)) > 1){
      stop("The transition matrix g must be a square matrix!")
    }
    if(dimOf_g[1] != length(weights)){
      stop("The dimensions of g and weights do not match!")
    }
    if(any(g < 0)){
      stop("The transition weights must be non-negative!")
    }
    for(i in 1:dim(g)[1]){
      if(sum(g[i, ]) > 1){
        stop("The sum of the transition weights of each hypothesis must not
           exceed 1!")
      }
    }
    if(any(weights < 0)){
      stop("The weights of the hypotheses must be non-negative!")
    }
    if(sum(weights) > 1){
      stop("The sum of the weights of the hypotheses must not exceed 1!")
    }
    m <- length(weights)
    if(m == 0){
      stop("The graph must contain at least one hypothesis") 
    }
  }else{
    m <- length(getWeights(gMCP))
  }
  if(!(length(mu_0) == 1 || length(mu_0) == m)){
    stop("The dimension of mu_0 must be 1 or equal to the number of
         hypotheses!")
  }
  if(any(mu_0 == -Inf)|| any(mu_0 == Inf)){
    stop("The entries of mu_0 must be real-valued!")
  }
  if(alpha <= 0 || alpha >= 1){
    stop("The significance level alpha must be strictly between 0 and 1!")
  }
}
###############################################################################
# Validation function for part of the input values
# of the informSCI-function, the inExactSCI-function and the 
# notInExactSCI-function 

checkPartOfInput2 <- function(m, q, estimates, Z, pValues, SE, I, randomShifts,
                              shifts, tolTrueSCI){
  if(!(length(q) == 1 || length(q) == m)){
    stop("The dimension of q must be 1 or equal to the number of hypotheses!")
  }
  if(any(q > 1) || any(q < 0)){
    stop("Each information weight must be between 0 and 1 (inclusive)!")
  }

  if(is.null(estimates) && is.null(Z) && is.null(pValues)){
    stop("Either estimates, Z or pValues must be specified.") 
  }
  if(!is.null(estimates) && !is.null(Z)){
    message("Only one parameter of either estimates, Z or pValues needs to be
            specified! Z is ignored, the algorithm uses estimates!") 
  }
  if(!is.null(estimates) && !is.null(pValues)){
    message("Only one parameter of either estimates, Z or pValues needs to be
            specified! pValues is ignored, the algorithm uses estimates!")
  }
  if(is.null(estimates) && !is.null(Z) && !is.null(pValues)){
    message("Only one parameter of either estimates, Z and pValues needs to be
            specified! pValues is ignored, the algorithm uses Z!")
  }
  if(!is.null(estimates) && m != length(estimates)){
    stop("The dimensions of the input values do not match.")
  }
  if(any(estimates == -Inf)|| any(estimates == Inf)){
    stop("The entries of estimates must be real-valued!")
  }
  if(!is.null(Z) && m != length(Z)){
    stop("The dimensions of the input values do not match.")
  }
  if(any(Z == -Inf)|| any(Z == Inf)){
    stop("The entries of Z must be real-valued!")
  }
  if(!is.null(pValues) && m != length(pValues)){
    stop("The dimensions of the input values do not match.")
  }
  if(any(pValues >= 1)|| any(pValues <= 0)){
    stop("The entries of pValues must be strictly between 0 and 1!")
  }
  if(is.null(SE) && is.null(I)){
    stop("Either SE or I must be specified!")
  }
  if(!is.null(SE) && !is.null(I)){
    message("Only one parameter of either SE or I needs to be specified!
            I is ignored, the algorithm uses SE!")
  }
  if(!is.null(SE) && length(SE) != 1 && length(SE) != m){
    stop("The dimension of SE must be 1 or equal to the number of hypotheses!")
  }
  if(any(SE <= 0)){
    stop("The entries of SE must be positive!")
  }
  if(any(SE == Inf)){
    stop("The entries of SE must be less than infinity!")
  }
  if(!is.null(I) && length(I) != 1 && length(I) != m){
    stop("The dimension of I must be 1 or equal to the number of hypotheses!")
  }
  if(any(I <= 0)){
    stop("The entries of I must be positive!")
  }
  if(any(I == Inf)){
    stop("The entries of I must be less than infinity!")
  }
  if(randomShifts%%1 != 0){
    stop("randomShifts must be an integer!")
  }
  if(any(shifts < 0)){
    stop("Each entry of shifts must be non-negative!") 
  }
  if(any(shifts == Inf)){
    stop("The entries of shift must be less than infinity!") 
  }
  if(!is.null(shifts)){ 
    if(is.null(dim(shifts))){
      if(length(shifts) != m){
        stop("The shift vectors must have length m!")  
      }
      if(sum(shifts) == 0){
        stop("Each shift vector must contain at least one positive entry!")
      }
    }else{
      if(dim(shifts)[2] != m){
        stop("The shifts matrix must have m columns!")
      }
      if(any(rowSums(shifts)==0)){
        stop("Each shift vector must contain at least one positive entry!")
      }
    }
  }
  if(tolTrueSCI <= 0){
    stop("The parameter tolTrueSCI must be positive!") 
  }
}
###############################################################################
# Validation function for part of the input values
# of the informSCI-function and the explore_q-function. 

checkPartOfInput3 <- function(timesSmallerEps, maxIter, maxIterBisec){
  if(timesSmallerEps <= 0){
    stop("The parameter timesSmallerEps must be greater than zero!")
  }
  if(maxIter < 1){
    stop("The maximum number of iterations for the entire algorithm must be
         at least 1!") 
  }
  if(maxIterBisec < 1){
    stop("The maximum number of iterations for the bisection must be
         at least 1!") 
  }
}
###############################################################################
# Validation function for the input values of the informSCI-function

checkSCIValues <- function(gMCP, g, weights, q, mu_0, estimates, Z, pValues, SE,
                           I, alpha, timesSmallerEps, randomShifts, shifts,
                           tolTrueSCI, maxIter, maxIterBisec) {
  checkPartOfInput1(gMCP, g, weights, mu_0, alpha)
  
  if(is.null(gMCP)){
    m <- length(weights)
  }else{
    m <- length(getWeights(gMCP))
  }
  checkPartOfInput2(m, q, estimates, Z, pValues, SE, I, randomShifts, shifts,
                    tolTrueSCI)
  checkPartOfInput3(timesSmallerEps, maxIter, maxIterBisec)
}
###############################################################################
# Validation function for the input values of the inExactSCI-function and the
# notInExactSCI-function
checkPrecisionTestInput <- function(L, randomShifts, shifts, tolTrueSCI, gMCP,
                                    g, weights, q, estimates, Z, pValues, SE, I,
                                    mu_0, alpha){
  checkPartOfInput1(gMCP, g, weights, mu_0, alpha)
  
  if(is.null(gMCP)){
    m <- length(weights)
  }else{
    m <- length(getWeights(gMCP))
  }
  
  checkPartOfInput2(m, q, estimates, Z, pValues, SE, I, randomShifts, shifts,
                    tolTrueSCI)
}
###############################################################################
# Validation function for part of the input values of the explore_q- function
# and the covMatrixManyToOne-function
checkPartOfInput4 <- function(sampleSizes, sampleSizeControl, varObs){
  if(any(c(sampleSizeControl, sampleSizes) <= 0)){
    stop("The sample size of each group needs to be greater than zero!")
  }
  if(any(c(sampleSizeControl, sampleSizes)%%1 != 0)){
    stop("The sample sizes must be natural numbers!")
  }
  if(!is.null(varObs) && varObs <= 0){
    stop("varObs must be positive!")
  }
}
###############################################################################
# Validation function for the input values of the explore_q-function:
checkSimulationValues <- function(gMCP, g, weights, trueParam, sigma, qFixed,
                                  mu_0, alpha, addHyp, allRej, atLeastOneRej,
                                  qInterval, qStepSize, qGrid, numSim,
                                  sampleSizes, sampleSizeControl, varObs,
                                  timesSmallerEps, maxIterSCI, maxIterBisec){
  
  checkPartOfInput1(gMCP, g, weights, mu_0, alpha)
  checkPartOfInput3(timesSmallerEps, maxIterSCI, maxIterBisec)
  
  if(is.null(gMCP)){
    m <- length(weights)
  }else{
    m <- length(getWeights(gMCP))
  }
  
  dimensions <- c(m, length(trueParam), nrow(sigma), ncol(sigma))
  if(length(unique(dimensions)) > 1){
    stop("The length of the vectors weights, trueParam and the number of 
         rows and columns of sigma do not match!")
  }
  if(any(abs(trueParam)==Inf)){
    stop("The entries of trueParam must be real-valued!")
  }
  if(!is.null(sigma)){
    if(any(t(sigma)!= sigma)){
      stop("The covariance matrix must be symmetric!")  
    }
    sigma.eigenValues <- eigen(sigma)$values
    if(any(sigma.eigenValues < 0)){
      stop("The covariance matrix must be positive semi-definite!")  
    }
  }

  if(! all(qFixed[ ,1] %in% 1:m)){ 
    stop("The first column of the qFixed-parameter may contain only indices 
         which are corresponding to the considered hypotheses!")
  }
  if(length(unique(qFixed[ ,1])) != length(qFixed[ ,1])){
    stop("The values of the fixed information weights may not be specified more
         than once!")
  }
  if(any(qFixed[ ,2] < 0) || any(qFixed[ ,2] > 1)){
    stop("Each fixed information weight must be between 0 and 1 (inclusive)!")
  } 
  if(any(addHyp[ ,1] %in% 1:m)){
    stop("Each hypothesis of the addHyp-matrix must have an identifier greater
         than m!")
  }
  if(length(unique(addHyp[ ,1])) < nrow(addHyp)){
    stop("Each hypothesis of addHyp needs an unique identifier!")
  }
  if(! all(addHyp[ ,2] %in% 1:m)){
    stop("The second column of the addHyp-parameter may contain only indices 
         which are corresponding to the parameters of interest!")
  }
  for(set in allRej){
    if(any(set < 1)|| any(set%%1 != 0)){
      stop("Each set of the allRej-list may contain only indices which
           are corresponding to the considered hypotheses!")
    }
  }
  for(set in atLeastOneRej){
    if(any(set < 1)|| any(set%%1 != 0)){
      stop("Each set of the atLeastOneRej-list may contain only indices
           which are corresponding to the considered hypotheses!")
    }
  }
  if(dim(qFixed)[1] != m){
    if(is.null(qGrid)){
      if(qInterval[1] < 0 || qInterval[2] > 1){
        stop("The simulation interval for q must be a subset of the closed 
          interval going from 0 to 1!") 
      }
      if(qInterval[1] >  qInterval[2]){
        stop("The entered lower bound for the information weight must be less than 
         or equal to the entered upper bound!") 
      }
      if(qStepSize <= 0){
        stop("The step size must be positive!") 
      }
    }else{
      message("The parameter qGrid is used. qInterval and qStepSize are
            ignored!")
      
      if(any(qGrid < 0) || any(qGrid> 1)){
        stop("The values for the varying information weights given by qGrid must 
            be between 0 and 1 (inclusive)!")
      }
    }
  }
  if(numSim < 1 || numSim%%1 !=0){
    stop("The number of simulations must be a natural number!")
  }
  if(is.null(sigma) && any(c(is.null(sampleSizes), is.null(sampleSizeControl),
                             is.null(varObs)))){
    stop("Either sigma or all three parameters sampleSizes, sampleSizeControl 
         and varObs need to be specified!")  
  }
  if(!is.null(sigma) && !all(c(is.null(sampleSizes), is.null(sampleSizeControl),
                               is.null(varObs)))){
    message("The parameter sigma is used. The parameters sampleSizes,
            sampleSizeControl and varObs are ignored!")
  }
  checkPartOfInput4(sampleSizes, sampleSizeControl, varObs)
}
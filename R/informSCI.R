# Implementation of the informative SCI algorithm and helper functions


#' Calculation of Lower Informative SCI-bounds
#'
#' @description
#' The function calculates informative lower SCI-bounds for a given graph of m
#' hypotheses and given information weights. m is a natural number.
#'
#' @details
#' It is assumed that there are m parameters of interest
#' \eqn{\vartheta_1,\dots,\vartheta_m}. For each parameter there is a null
#' hypothesis defined as \eqn{H_i^{{\mu_0}_i}:\vartheta_i\leq{\mu_0}_i}.
#' The bounds \eqn{{\mu_0}} correspond to \code{mu_0}. The parameter \code{gMCP} 
#' or the parameters \code{g} and \code{weights} define a graphical multiple 
#' test procedure for the hypotheses.
#' 
#' The algorithm further assumes that for each hypothesis there is an unbiased 
#' point estimator which is (asymptotically) normal. 
#' The \code{informSCI}-algorithm is based on the p-values from the
#' corresponding (asymptotic) z-tests.
#' 
#' The algorithm terminates when the Chebyshev distance of two successive 
#' calculated approximations is less than \code{eps}
#' \code{timesSmallerEps}-times in succession or if the maximum number of
#' iterations is reached. 
#' 
#' The function also tries to find information about the precision of the final
#' calculated approximation of the true lower informative SCI-bounds by
#' calling the \code{\link{inExactSCI}}- and the
#' \code{\link{notInExactSCI}}-functions. 
#' 
#' For further details see the given references.
#'
#' @param gMCP An object of class \code{\link[gMCP]{graphMCP}} indicating the
#' underlying graphical test.
#' @param g  Numeric square matrix of transition weights for the graphical test
#' with m rows and m columns. The i-th row of the entered matrix defines the
#' arrows starting from the i-th hypothesis. Each entry has to be between
#' 0 and 1 and each row must sum to a number less than or equal to 1. The 
#' diagonal elements must be zero. Entering \code{g} and \code{weights} can be
#' used as the input as an alternative to specifying \code{gMCP}. 
#' @param weights Numeric vector of weights of dimension m. It defines the 
#' initial proportion of significance level which is assigned to each null 
#' hypothesis. Entering \code{g} and \code{weights} can be used as the input as
#' an alternative to specifying \code{gMCP}.
#' @param q  A numeric vector of dimension 1 or m defining the information 
#' weights for each hypothesis. The entries have to be between 0 and 1
#' (inclusive). If \code{q} has dimension 1, the same information weight is 
#' used for each hypothesis.
#' @param mu_0 A numeric vector of dimension 1 or m defining the bounds of the
#' null hypotheses of the underlying graphical test. If \code{mu_0} has 
#' dimension 1, the same value is used for each null hypothesis.
#' @param estimates An m-dimensional numeric vector of unbiased point estimates
#' for the parameters of interest. Each estimator is assumed to be
#' (asymptotically) normal.
#' @param Z An m-dimensional numeric vector of z-scores for testing the null 
#' hypotheses. \code{Z} can be entered as an alternative to \code{estimates}.
#' @param pValues An m-dimensional numeric vector of p-values from (asymptotic)
#' z-tests for testing the null hypotheses. \code{pValues} can be entered as an 
#' alternative to \code{estimates} and \code{Z}.
#' @param SE A positive numeric vector of dimension 1 or m indicating the 
#' standard errors of the point estimators. If \code{SE} has dimension 1, the
#' same value is used for each estimator.
#' @param I A positive numeric vector indicating the information of the
#' estimators. It can be entered as an alternative to the vector \code{SE}.
#' The relationship \code{SE}\eqn{=1/}\code{I}\eqn{^{1/2}} is assumed. 
#' @param alpha A numeric defining the overall significance level for the
#' graphical test (i.e. SCIs will have coverage probability of at least 
#' \eqn{1-}\code{alpha}. The parameter must be strictly between 0 and 1.
#' @param eps A numeric indicating the desired strict upper bound on the 
#' Chebyshev distance between two successive calculated approximations (the
#' Chebyshev distance is induced by the maximum norm).
#' @param timesSmallerEps A positive integer indicating how many times the 
#' Chebyshev distance of two successive calculated approximations should be 
#' less than \code{eps} in succession. Here we use the convention 
#' \eqn{|-\infty-(-\infty)|:=0}. 
#' @param randomShifts A positive integer indicating how many random directions
#' of length \code{tolTrueSCI} should be generated. At the end of the algorithm 
#' the parameter is passed to the \code{\link{inExactSCI}}- and the 
#' \code{\link{notInExactSCI}}-functions to determine whether the approximation 
#' \code{L} of the true lower SCI-bounds is precise or imprecise.
#' @param shifts A matrix with m columns and any number of rows. Each entry must 
#' be non-negative. Each row is a direction in the m-dimensional real space.
#' Each row must have at least one positive entry. At the end of the algorithm 
#' the parameter is passed to the \code{\link{inExactSCI}}- and the 
#' \code{\link{notInExactSCI}}-functions to determine whether the approximation 
#' \code{L} of the true lower SCI-bounds is precise or imprecise. 
#' \code{randomShifts} must be a positive integer or \code{shifts} must contain
#' at least one row. It is recommended to choose \code{randomShifts}\eqn{>0}
#' or enter a \code{shifts}-matrix with at least one row. Entering both
#' parameters is also permissible.
#' @param tolTrueSCI The randomly generated shift-vectors and the row vectors 
#' in \code{shifts} are rescaled to have length \code{tolTrueSCI}. It is
#' recommended to choose \code{tolTrueSCI} greater than
#' \eqn{\sqrt{m}\cdot}\code{eps} and greater than \code{tolBisec}.
#' @param maxIter Maximum number of iterations for determining the lower 
#' informative SCI bounds.
#' @param maxIterBisec Maximum number of iterations of the bisection method
#' which is used during the algorithm for finding roots.
#' @param tolBisec A non-negative numeric indicating the error tolerance of the
#' bisection method which is used for finding roots.
#' @param calculateCSCI A boolean indicating whether compatible bounds should 
#' also be calculated.
#' @param checkInput A boolean specifying whether the entered values should be 
#' checked.
#'
#' @import gMCP
#' @importFrom stats qnorm
#' @return The function returns a list containing the calculated lower 
#' informative SCI-bounds as well as compatible lower SCI-bounds 
#' (if \code{calculateCSCI==TRUE}) to allow a comparison. Additionally, the 
#' returned list contains elements which can give some information about the
#' precision of the calculated lower informative SCI-bounds compared to the
#' true informative SCI-bounds.
#' \itemize{
#' \item \code{L}: A numeric vector of dimension m of the lower informative 
#' SCI-bounds
#' \item \code{rejecHyp}: A boolean vector of dimension m indicating the rejected
#' hypotheses of the multiple test induced by the informative SCI-bounds
#' \item \code{diffApprox}: A non-negative numeric indicating the Chebyshev distance
#' between the calculated last two approximations of the lower informative
#' SCI-bounds.
#' \item \code{timesApprSmallerEps}: A natural number between 0 and 
#' \code{timesSmallerEps} indicating how many times the Chebyshev distance of 
#' two successive calculated approximations in succession  was less than
#' \code{eps} when the algorithm terminated.
#' \item \code{numIter}: A natural number indicating the number of iterations
#' required by the algorithm. 
#' \item \code{accuracyL}: A string containing information about the collected 
#' information about the precision of the calculated lower informative 
#' SCI-bounds compared with the true lower SCI-bounds.
#' }
#' 
#' If \code{calculateCSCI=TRUE}:
#' \itemize{
#' \item \code{informSCIcompatible}: A boolean vector of dimension m indicating
#' whether each informative bound is compatible to the test decision
#' about its corresponding null hypothesis by the underlying graphical test.
#' \item \code{cSCI}: A numeric vector of dimension m of compatible lower
#' SCI-bounds from \code{\link[gMCP]{simConfint}}.
#' \item \code{rejecHypGraph}: A boolean vector of dimension m indicating the 
#' rejected hypotheses of the underlying graphical test.
#' }
#' @export
#' @references F. Bretz, W. Maurer, W. Brannath, M. Posch: A graphical approach
#' to sequentially rejective multiple test procedures. Statistics in Medicine
#' 28.4 (2009), pp. 586-604.
#' 
#' K. Strassburger, F. Bretz: Compatible simultaneous lower confidence bounds 
#' for the Holm procedure and other Bonferroni based closed tests. Statistics 
#' in Medicine 27.4 (2008), pp. 4914-4927.
#' 
#' S. Schmidt, W. Brannath: Informative Simultaneous Confidence Intervals
#' in Hierarchical Testing. Methods of Information in Medicine 53.4 (2014), 
#' pp. 278–283.

#' @seealso \code{\link[gMCP]{gMCP}} \code{\link[gMCP]{simConfint}} 
#' \code{\link{explore_q}}
#' @examples 
#' informSCI(gMCP=BonferroniHolm(3), q=0.3, mu_0=-0.5 ,estimates=c(0,2,-1),
#' SE=0.1467)
#' @examples
#' Z <- (c(0,2,-1)-(-0.5))/0.1467
#' informSCI(gMCP=BonferroniHolm(3), q=0.3, mu_0=-0.5, Z=Z, I=1/(0.1467^2),
#' randomShifts=100)
#' @examples
#' informSCI(g=matrix(c(0,0,1,0),2,2), weights=c(1,0), q=c(0.0068,1),
#' mu_0=c(-1,0), pValues=c(0.0002,0.01), SE=c(0.31,1.11), alpha=0.025, 
#' shifts=rbind(c(1,0),c(0,1),c(1,1)))
#' @examples
#' informSCI(g=matrix(c(0,0,1,0),2,2), weights=c(1,0), q=c(0.0068,1),
#' mu_0=c(-1,0), pValues=c(0.0002,0.01), I=1/c(0.31,1.11)^2, alpha=0.025, 
#' shifts=rbind(c(1,0),c(0,1),c(1,1)), calculateCSCI = FALSE)
#' 
informSCI <- function(gMCP = NULL, g = NULL, weights = NULL, q, mu_0 = 0,
                      estimates = NULL, Z = NULL, pValues = NULL, SE = NULL,
                      I = NULL, alpha = 0.05, eps = 1/10^5, timesSmallerEps = 3,
                      randomShifts = 0, shifts = NULL, 
                      tolTrueSCI = sqrt(ifelse(!is.null(gMCP),
                                               length(getWeights(gMCP)),
                                               length(weights)))*eps,
                      maxIter = 1000, maxIterBisec = 1000, tolBisec = 1/10^5,
                      calculateCSCI = TRUE, checkInput = TRUE)
  {
  ############# check and process parameters: ###########
  
  if(checkInput){
    checkSCIValues(gMCP, g, weights, q, mu_0, estimates, Z, pValues, SE, I, 
                   alpha, timesSmallerEps, randomShifts, shifts, tolTrueSCI,
                   maxIter, maxIterBisec)  
  }
  
  if(!(is.null(gMCP))){
    g <- getMatrix(gMCP)
    weights <- getWeights(gMCP) 
  }
  if(!is.null(I) && is.null(SE)){
    SE <- 1/sqrt(I)
  }
  if(!is.null(Z) && is.null(estimates)){
    estimates <- Z*SE + mu_0
  }
  if(!is.null(pValues) && is.null(estimates) && is.null(Z)){
    estimates <- qnorm(1-pValues) * SE + mu_0 
  }
  
  m <- length(weights) # number of hypotheses 
  
  # Check whether the same q, mu_0 and/ or SE should be used for all hypotheses:
  if(length(q) == 1){
    q <- rep(q,m)
  }
  if(length(mu_0) == 1){
    mu_0 <- rep(mu_0,m)
  }
  if(length(SE) == 1){
    SE <- rep(SE,m)
  }
  ############# calculation of SCI  ###########
  
  output <- list()
  if(all(q==0)){ # perform graphical test:
    output <- graphTestAndCompSCI(g, weights, estimates, SE, mu_0, alpha)  
  }else{# calculation of lower informative SCI-bounds
    out_SCIwBisec <- SCIwithBisection(g, weights, q, estimates, SE, mu_0, alpha,
                                      eps, timesSmallerEps, maxIter,
                                      maxIterBisec, tolBisec) 
    L <- out_SCIwBisec$L 
    
    exactApprox <- inExactSCI(L = L, randomShifts = randomShifts,
                              shifts = shifts, tolTrueSCI = tolTrueSCI, g = g,
                              weights = weights,  q = q, estimates = estimates,
                              SE = SE,  mu_0 = mu_0, alpha = alpha,
                              checkInput = FALSE) 
    
    noExactApprox <- notInExactSCI(L = L, randomShifts = randomShifts,
                                   shifts = shifts, tolTrueSCI = tolTrueSCI,
                                   g = g, weights = weights,  q = q, 
                                   estimates = estimates, SE = SE, mu_0 = mu_0,
                                   alpha = alpha, checkInput = FALSE) 
    
    if(calculateCSCI){
      original_gMCP <- graphTestAndCompSCI(g, weights, estimates, SE, mu_0,
                                           alpha)#The underlying graphical 
      # test is performed and compatible lower SCI bounds are calculated.
      original_gMCPcalculated <- TRUE
    }else{
      original_gMCP <- "notCalculated"
      original_gMCPcalculated <- FALSE  
    }
    
    output <- makeReturnList(out_SCIwBisec, mu_0, exactApprox, noExactApprox,
                             original_gMCPcalculated, original_gMCP)  
  }
  return(output)
}
###############################################################################
########################### helper functions ##################################
###############################################################################

###############################################################################
# This function calculates p-values for right-tailed Z-Tests.
#' @importFrom stats pnorm

pvalueRTzTest <- function(estimates, mu, SE){
  return(1 - pnorm((estimates - mu)/SE))
}

###############################################################################
# function for calculating lower confidence bounds with coverage probability
# 1-gamma[i] for \vartheta_i, 1 <= i <= m.

#' @importFrom stats qnorm
inversePValueZTest <- function(estimates, SE, gamma){
  return(estimates - SE * qnorm(1 - gamma))
}
###############################################################################
# function for adapting the information weights if the sums of the rows of the
# transition weights matrix are less than 1. 

informWeight_adapted <- function(mu, mu_0 , q, trWeights.rowsums){
  greaterZero <- as.integer(mu > mu_0) * as.integer(q > 0) 
  return((1 - (1-q^{mu-mu_0}) * trWeights.rowsums) * greaterZero)
}

###############################################################################
# This function defines one side of the key equation which needs to be solved
# during the informSCI-algorithm 

partOfKeyEquation <- function(z, mu_0, eta.mu, q, g.rowsums){
  return(eta.mu * ifelse(z-mu_0 > 0,
                         informWeight_adapted(z, mu_0, q, g.rowsums), 1))
}
###############################################################################
#' Bisection function
#' 
#' @description
#' Bisection function to find solutions of the key equation of the
#' \code{informSCI}-algorithm. 
#' 
#' @details
#' The function tries to find a solution of the key equation of the
#' \code{informSCI}-algorithm which is equivalent to  determining the
#' intersection point of \code{f_1} and \code{f_2}.
#' The function uses the bisection method and tries to determine the root
#' of the function \code{f_1-f_2}. Note that by definition of the key equation
#' and the assumptions of the \code{informSCI}-algorithm \code{f_1-f_2} is 
#' a continuous strictly increasing function. Because of the assumptions on
#' \code{a} and \code{b} \code{f_1-f_2} has a non-positive function value in
#' point \code{a} and non-negative function value in point \code{b}. Thus,
#' \code{f_1-f_2} has exactly one root in the closed interval \eqn{[a,b]}.
#' 
#' The bisection method repeatedly halves the interval between \code{a} and
#' \code{b}. The function stops when the root is found or when the maximum 
#' number of iterations is reached or when the interval is less than \code{tol}.
#' 
#' @param f_1 Left side of the key equation as a function in one variable.
#' @param f_2 Right side of the key equation as a function in one variable.
#' @param a A real value indicating the left bound of the search region. 
#' \eqn{f_1(a)\leq f_2(a)} must hold true.
#' @param b A real value indicating the right bound of the search region. 
#' \eqn{f_1(b)\geq f_2(b)} must hold true.
#' @param maxIter A positive integer defining the maximum number of iterations.
#' @param tol A non-negative numeric indicating the error tolerance.
#'
#' @keywords internal
#' @return Returns intersection point. In the case that no intersection point
#' is found, the left side of the final interval is returned, rather than the
#' midpoint. The returned point is a lower approximation of the solution of the
#' key equation.
#'
funcBisec <- function(f_1, f_2, a, b, maxIter = 1000, tol = 1/10^3){
  # finding intersection between f_1 and f_2 is equivalent to finding root of f:
  f <- function(x){f_1(x)-f_2(x)}
  f.a <- f(a)
  f.b <- f(b)
  if(f.a > 0 || f.b < 0 ){ 
    stop("Bisection cannot be executed! The function value f_1(a) must be less
         than or equal to f_2(a) and f_1(b) must be greater than or equal
         to f_2(b)!") 
  }
  if(f.a == 0){
    return(a)
  }
  
  check <- f_2(b) # for numerical stability 
  if(f.b == 0 && check != 0){
    return(b)
  }
  
  for(i in 1:maxIter){
    if(abs(b - a) < tol){ 
      break
    }
    xmid <- (a + b)/2
    ymid <- f(xmid)
    if (ymid < 0) { # ymid < 0 means f_1(x_mid) < f_2(xmid)
      a <- xmid
    }
    else if(ymid == 0){ # i.e. f_1(x_mid) == f_2(x)
      check <-  f_2(xmid) 
      if(check != 0){ # for numerical stability
        a <- xmid
        break  
      }else{
        b <- xmid
      }
    }else{ # i.e. f_1(x_mid) > f_2(x)
      b <- xmid
    }
  }
  return(a)
}
###############################################################################
#' Function for determining the (monotone part of the) local significance levels
#'
#' @description
#' Function for determining the monotone part (\code{eta.mu}) of the local 
#' significance levels for the key equation of the informative SCI algorithm.
#' The function creates dual graphs and rejects some of its hypotheses to
#' obtain the local significance levels.
#' 
#' @details
#' m = number of hypotheses.
#' 
#' The function is not suitable if for all \eqn{1\leq i\leq m} it holds
#' \code{q[i]==0} and \code{mu[i]>mu_0[i]}.
#'  
#' @param mu A real-valued vector (-\code{Inf} is also allowed) of dimension m 
#' indicating which dual graph should be created and which null hypotheses 
#' should be rejected. \code{mu[i]>mu_0[i]} iff the corresponding hypothesis is 
#' rejected, \eqn{1\leq i\leq m}.
#' @param g A numeric square matrix of transition weights for the graphical
#' test procedure.
#' @param weights A numeric vector of dimension m of initial weights for the
#' graphical test procedure.
#' @param alpha Overall level of the graphical test procedure.
#' @param q A numeric vector of dimension m of information weights.
#' @param mu_0 A numeric vector of dimension m of bounds of the null hypotheses.
#'
#' @return Returns a numeric vector of dimension m (\code{eta.mu}) used for
#' solving the key equation of the \code{informSCI} algorithm. It contains the
#' local levels in \code{mu} divided by \code{q^{max(mu-mu_0,0)}} or divided by
#' adapted information weights (only if \code{q[i]>0}).
#' 
#' @keywords internal
weightsGTP <- function(mu, g, weights, alpha, q, mu_0){ 
  
  m <- length(weights)
  
  q.zero = which(q == 0)
  q.eqZeroIndicator <- rep(1,m)
  q.eqZeroIndicator[q.zero] <- 0
  
  indices.pos <- which(mu-mu_0 > 0) 
  numOfPosMu <- length(indices.pos) # number of hypotheses that should be
  # rejected
  indices.nonpos <- setdiff(1:m, indices.pos)
  
  if(numOfPosMu == 0){
    return(weights * alpha) 
    
  }else if(numOfPosMu == 1){
    i.p <- indices.pos[1] # single index for which mu[i]>mu_0[i] is true
    eta.mu <- numeric(m)
    alpha_posInd <- weights[i.p] * alpha 
    eta.mu[i.p] <- alpha_posInd * q.eqZeroIndicator[i.p] 
    eta.mu[indices.nonpos] <- weights[indices.nonpos] * alpha +
      (1-q[i.p]^{mu[i.p]-mu_0[i.p]}) * g[i.p,indices.nonpos] * alpha_posInd
    return(eta.mu)
    
  }else{
    indices.sorted <- c(indices.pos,indices.nonpos)
    
    # modify graph to get graph G^mu:
    
    g.red <- g[indices.pos, ] # only relevant transition weights are considered
    g.red <- (1-q[indices.pos]^{mu[indices.pos]-mu_0[indices.pos]}) * g.red 
    
    g.red <- g.red[ ,indices.sorted] # g.red is rearranged to simplify
    # following calculations
    d <- diag(q.eqZeroIndicator[indices.pos])
    g.mu <- cbind(g.red, d) # shifted hypotheses are added
    
    # g.mu is a block matrix comprised of three matrices in a row. The first
    # matrix contains all transition weights going from any H_i^{mu_0_i}
    # fulfilling mu_i > mu_0_i to any H_j^{mu_0_j} with mu_j > mu_0_j
    # (i.e. i,j in indices.pos). The second matrix contains all transition
    # weights going from any H_i^{mu_0_i} with mu_i > mu_0_i to any H_k^{mu_k}
    # with mu_k <= mu_0_k (i.e. i in indices.pos, k in indices.nonpos).
    # The third matrix consists of the transition weights going from
    # any H_i^{mu_0_i} with mu_0_i > mu_i to H_i^{mu_i}. Here we use
    # q.eqZeroIndicator[i] as transition weights instead of
    # q[i]^{mu_i-mu_0_i}*q.eqZeroIndicator[i] to stabilize the calculations.
    eta.mu <- c(alpha * weights[indices.sorted], numeric(numOfPosMu))
    
    # update significance levels:
    for(l in 1:numOfPosMu){ 
      
      eta.mu[(l+1):(m+l)] <- eta.mu[(l+1):(m+l)]+eta.mu[l] * g.mu[l,(l+1):(m+l)] 
      
      # update transition weights:
      if(l < numOfPosMu){ # when the last hypothesis has been rejected, no
        # further update of the transition weights is necessary
        
        for(j in (l+1):numOfPosMu){ # only the transition weights between 
          # non-rejected hypotheses need to be updated 
          for(k in (l+1):(m+j)){ # for m+j+1 <= k <= m+numOfPosMu the
            # transition weights remain 0. Thus no update is necessary
            
            if(g.mu[j,l] * g.mu[l,j] == 1 || j == k){
              g.mu[j,k] <- 0
            }else{
              g.mu[j,k] <- (g.mu[j,k] + g.mu[j,l] * g.mu[l,k])/
                (1 - g.mu[j,l] * g.mu[l,j])
            } 
          }
        }
      }
    }
    eta.mu <-eta.mu[(numOfPosMu+1):(numOfPosMu+m)] 
    
    eta.mu <-eta.mu[order(c(indices.nonpos,indices.pos))] # rearrange the 
    # entries of the vector into the original arrangement
    
    #eta.mu gives alpha^{mu}/q^{max(0,mu-mu_0)} or
    # alpha^{mu}/informWeight_adapted(mu, mu_0 , q, rowSums(g)) where alpha^{mu}
    # is defined as in the underlying paper. If q[i] == 0 and mu[i]< mu_0[i] it
    # holds eta.mu[i] = alpha_i^{mu}. If q[i] == 0 and mu[i]>=mu_0[i] it holds 
    # eta.mu[i] = 0.
    return(eta.mu)  
  }
}
###############################################################################
# The function calculates the lower informative SCI-bounds. The function is not
# suited when all entries of q are equal to zero. 

SCIwithBisection <- function(g, weights, q, estimates, SE, mu_0, alpha = 0.05,
                             eps = 1/10^5, timesSmallerEps = 3, maxIter = 1000,
                             maxIterBisec = 1000, tolBisec = 1/10^5){
  
  m <- length(weights)
  g.rowsums <- rowSums(g)
  localLevels <- alpha * weights
  
  q.zero = which(q == 0)
  
  # determine lower estimate for SCI:
  
  stabilConstant <- 1/10^3 # for numerical stability, to prevent certain
  # rounding errors from sabotaging the result
  
  mu <- pmin(inversePValueZTest(estimates, SE, localLevels), mu_0) -
    stabilConstant 
  mu_new <- mu

  eta.mu <- localLevels
  
  # determine upper estimate for SCI by calculating the un-adjusted bounds
  biggestUpEstSCI <- inversePValueZTest(estimates, SE, alpha*rep(1,m))
  
  counterConvergence <- 0
  numIter <- 0
  lowEstSCI <- numeric(m)
  
  for(i in 1:maxIter){
    indexset <- which(eta.mu > 0) # If eta.mu[j] == 0 and q[j] > 0 or
    # if eta.mu[j] == 0, q[j] == 0 and mu[j] != mu_0[j], no improvement of
    # mu[j]=-Inf is possible. Hence, it remains mu_new[j] = mu[j] = -Inf.
    # Furthermore, if q[j] == 0 and mu[j] == mu_0[j], no improvement is possible
    # either. In this case it holds that eta.mu[j]==0 due to the update of 
    # eta.mu at the end of this function. In conclusion, if eta.mu[j] == 0,
    # no update step is done. 
    
    for(j in indexset){
      
      # if q[j] == 0 or q[j] == 1, no bisection is needed: 
      if(q[j] == 0){ 
        mu_new[j] <- min(inversePValueZTest(estimates[j], SE[j], eta.mu[j]),
                         mu_0[j])
      }else if(q[j] == 1){
        mu_new[j] <- inversePValueZTest(estimates[j], SE[j], eta.mu[j])
      }else{
        # for using bisection we need a lower bound and an upper bound
        # of the root:
        if(mu[j] == -Inf){
          lowEstSCI[j] <- min(inversePValueZTest(estimates[j], SE[j], eta.mu[j]),
                           mu_0[j])-stabilConstant
        }else{
          lowEstSCI[j] <- mu[j]
        }
        f_1 <- function(z){pvalueRTzTest(estimates[j], z, SE[j])}
        f_2 <- function(z){partOfKeyEquation(z, mu_0[j], eta.mu[j], q[j],
                                             g.rowsums[j])}
        mu_new[j] <- funcBisec(f_1, f_2, lowEstSCI[j], biggestUpEstSCI[j],
                               maxIterBisec, tolBisec) 
      }
    }
    
    numIter <- numIter + 1
    
    diffAppr <- max(abs(mu_new - mu),na.rm=TRUE)
    if(diffAppr < eps){ 
      counterConvergence <- counterConvergence + 1
      if(counterConvergence >= timesSmallerEps){ # checks for convergence
        break
      }
    }else{
      counterConvergence <- 0
    }
    mu <- mu_new
    
    # In the case that q[j]==0 and mu[j]==mu_0[j], i.e. H_j^mu_0_j has been
    # rejected, the whole level of H_j^mu_0_j should be given to the other
    # hypotheses. To achieve this, we pass a value greater than mu_0[j] to the
    # weightsGTP-function due to the construction of the weightsGTP-function.
    # It does not matter how much larger the value is. We input mu[j]+1.
    # The other components of mu remain unchanged:
    indicesToGetCorrectWeights <- intersect(which(mu==mu_0),q.zero)
    muToGetCorrectWeights <- mu
    muToGetCorrectWeights[indicesToGetCorrectWeights] <-
      muToGetCorrectWeights[indicesToGetCorrectWeights] + 1
    
    eta.mu <- weightsGTP(muToGetCorrectWeights, g, weights, alpha, q, mu_0)
  }
  
  # generate return list
  returnList <- list("L" = mu_new,
                     "diffAppr" = diffAppr,
                     "counterConvergence" = counterConvergence, 
                     "numIter" = numIter)
  return(returnList)
  
}
###############################################################################
# The function performs a graph based multiple test procedure for a given graph
# and given estimates. Furthermore it calculates corresponding compatible SCI
# bounds. The function uses the gMCP-package
#'
#' @import gMCP
#'
graphTestAndCompSCI <- function(g, weights, estimates, SE, mu_0, alpha){

  
  pValues <- pvalueRTzTest(estimates = estimates, mu = mu_0, SE)
  graph <- matrix2graph(g, weights)
  
  result <- gMCP(graph = graph, pvalues = pValues, alpha = alpha) # perform 
  # graphical test
  rejHyp <- getRejected(result)
  indices.rejHyp <- which(rejHyp == TRUE)
  finalLevel <- getWeights(result) * alpha 
  
  if(all(rejHyp == TRUE)){
    L <- pmax(inversePValueZTest(estimates,SE,gamma = alpha*weights), mu_0) 
  }else{
    L <- inversePValueZTest(estimates, SE, gamma = finalLevel)
    L[indices.rejHyp] <- mu_0[indices.rejHyp] # bounds for rejected hypotheses
  }
  names(L) <- sprintf("L%d", 1:length(weights))
  returnList <- list("rejecHyp" = rejHyp, "L" = L)
  return(returnList)
}
###############################################################################
#' Checking Precision of Approximations
#' 
#' @description
#' The functions checks whether information about the precision of an 
#' approximation for the informative lower SCI-bounds can be collected.
#' 
#' @details
#' The function checks if it can be determined whether \code{L} can be shifted
#' by a randomly generated rescaled direction or by a rescaled direction in the
#' shift matrix such that it lies in the true SCI. If this is possible,
#' the approximation is precise.
#' (The random directions are generated in such a way that all entries are 
#' positive.)
#' 
#' Let m be the dimension of \code{L}. m also describes the number of 
#' hypotheses of interest.
#' 
#' @param L An m-dimensional non-negative vector whose entries are the lower
#'  bounds of an approximation of the informative SCI.
#' @param randomShifts A positive integer indicating how many random directions
#' of length \code{tolTrueSCI} should be generated. 
#' @param shifts A matrix with m columns and any number of rows. Each entry must 
#' be non-negative. Each row is a direction in the m-dimensional real space.
#' Each row must have at least one positive entry. \code{randomShifts} should 
#' be a positive integer or \code{shifts} should contain at least one row.
#' @param tolTrueSCI The randomly generated shift-vectors and the row vectors
#' in \code{shifts} are rescaled to have length \code{tolTrueSCI}.
#' @param gMCP An object of class \code{\link[gMCP]{graphMCP}} indicating the
#' underlying graphical test.
#' @param g  Numeric square matrix of transition weights for the graphical test
#' with m rows and m columns. The i-th row of the entered matrix defines the
#' arrows starting from the i-th hypothesis. Each entry has to be between
#' 0 and 1 and each row must sum to a number less than or equal to 1. The 
#' diagonal elements must be zero. Entering \code{g} and \code{weights} can be
#' used as the input as an alternative to specifying \code{gMCP}. 
#' @param weights Numeric vector of weights of dimension m. It defines the 
#' initial proportion of significance level which is assigned to each null 
#' hypothesis. Entering \code{g} and \code{weights} can be used as the input as
#' an alternative to specifying \code{gMCP}.
#' @param q  A numeric vector of dimension 1 or m defining the information 
#' weights for each hypothesis. The entries have to be between 0 and 1
#' (inclusive). If \code{q} has dimension 1, the same information weight is 
#' used for each hypothesis.
#' @param estimates An m-dimensional numeric vector of unbiased point estimates
#' for the parameters of interest. Each estimator is assumed to be
#' (asymptotically) normal.
#' @param Z An m-dimensional numeric vector of z-scores for testing the null 
#' hypotheses. \code{Z} can be entered as an alternative to \code{estimates}.
#' @param pValues An m-dimensional numeric vector of p-values from (asymptotic)
#' z-tests for testing the null hypotheses. \code{pValues} can be entered as an 
#' alternative to \code{estimates} and \code{Z}.
#' @param SE A positive numeric vector of dimension 1 or m indicating the 
#' standard errors of the point estimators. If \code{SE} has dimension 1, the
#' same value is used for each estimator.
#' @param I A positive numeric vector indicating the information of the
#' estimators. It can be entered as an alternative to the vector \code{SE}.
#' The relationship \code{SE}\eqn{=1/}\code{I}\eqn{^{1/2}} is assumed.  
#' @param mu_0 A numeric vector of dimension 1 or m defining the bounds of the
#' null hypotheses of the underlying graphical test. If \code{mu_0} has 
#' dimension 1, the same value is used for each null hypothesis.
#' @param alpha A numeric defining the overall significance level for the
#' graphical test (i.e. SCIs will have coverage probability of at least 
#' \eqn{1-}\code{alpha}. The parameter must be strictly between 0 and 1.
#' @param checkInput A boolean specifying whether the entered values should be 
#' checked.
#'
#' @import mvtnorm
#' @importFrom stats qnorm
#' @return Returns \code{TRUE} if we can determine that the approximation is
#' indeed precise. Returns \code{FALSE} if we cannot determine that the
#' approximation is precise. (The approximation may still be precise.)
#' @export
#'
#' @seealso \code{\link{informSCI}} \code{\link{explore_q}}
#' @examples
#' g <- matrix(c(0,0,1,0),2,2)
#' weights <- c(1,0)
#' q <- c(0.0068,1)
#' mu_0 <- c(-1,0)
#' pValues <- c(0.0002,0.01)
#' SE <- c(0.31,1.11)
#' alpha <- 0.025
#' L <- informSCI(g=g, weights=weights, q=q, mu_0=mu_0, pValues=pValues, SE=SE,
#' alpha=alpha, eps=1/10^5, tolBisec=1/10^5)$L
#' # When the randomShifts- or shift-parameter in the informSCI-function is
#' # specified, the inExactSCI-function is called by the informSCI-function.
#' # It is also possible to analyse the accuracy of a calculated L (or an 
#' # approximation of the lower informative SCI-bounds) by directly using 
#' # the inExactSCI-function:
#' inExactSCI(L=L, randomShifts=100, tolTrueSCI=1/10^5, g=g, weights=weights,
#' q=q, pValues=pValues, SE=SE, mu_0=mu_0, alpha=alpha)

inExactSCI <- function(L, randomShifts = 0, shifts = NULL, tolTrueSCI,
                       gMCP = NULL, g = NULL, weights = NULL, q,
                       estimates = NULL, Z = NULL, pValues = NULL, SE = NULL,
                       I = NULL, mu_0, alpha, checkInput = TRUE){
  
  if((randomShifts <= 0) && is.null(shifts)){
    return(FALSE)
  }

  checkPrecisionTestInput(L, randomShifts, shifts, tolTrueSCI, gMCP, g, weights,
                          q, estimates, Z, pValues, SE, I, mu_0, alpha)
  
  # check and process parameters:
  
  if(!(is.null(gMCP))){
    g <- getMatrix(gMCP)
    weights <- getWeights(gMCP) 
  }
  m <- length(weights) # number of hypotheses 
  
  if(!is.null(I) && is.null(SE)){
    SE <- 1/sqrt(I)
  }
  if(!is.null(Z) && is.null(estimates)){
    estimates <- Z*SE + mu_0
  }
  if(!is.null(pValues) && is.null(estimates) && is.null(Z)){
    estimates <- qnorm(1-pValues) * SE + mu_0 
  }
  
  # Check whether the same q, mu_0 and/ or SE should be used for all hypotheses:
  if(length(q) == 1){
    q <- rep(q,m)
  }
  if(length(mu_0) == 1){
    mu_0 <- rep(mu_0,m)
  }
  if(length(SE) == 1){
    SE <- rep(SE,m)
  }
  
  if(randomShifts > 0){
    randShifts <- abs(rmvnorm(n = randomShifts, mean = numeric(m)))
    while(any(rowSums(randShifts)==0)){
      randShifts <- abs(rmvnorm(n = randomShifts, mean = numeric(m)))  
    }
    
    shifts <- rbind(shifts,randShifts)
  } 
  
  if(is.null(dim(shifts))){ # needed when shifts is only a vector (only one
    # direction is considered
    shifts <- matrix(shifts, nrow =  1, ncol = m) 
  }
  
  g.rowsums <- rowSums(g)
  numShifts <- dim(shifts)[1]
  
  # rescaling of the rows of shifts:
  for(i in 1:numShifts){
    shifts.norm <- norm(shifts[i, ], type = '2')
    shifts[i, ] <- shifts[i, ] * tolTrueSCI/shifts.norm
  }
  counter <- 1
  precApprox <- FALSE
  
  while(counter <= numShifts && precApprox == FALSE){
    L_shifted <- L + shifts[counter, ]
    eta.L_shifted <- weightsGTP(L_shifted, g, weights, alpha, q, mu_0)
    
    equation_L_shifted <- pvalueRTzTest(estimates, L_shifted, SE)-
      partOfKeyEquation(L_shifted, mu_0, eta.L_shifted, q, g.rowsums)
    
    # The indices for which the local level in L_shifted is zero,
    # must be excluded:
    indicesLocLevelEqZeroPartOne <- which(eta.L_shifted == 0)
    indicesLocLevelEqZeroPartTwo <- intersect(which(L_shifted > mu_0),
                                              which(q==0))
    indicesLocLevelEqZero <- union(indicesLocLevelEqZeroPartOne,
                                   indicesLocLevelEqZeroPartTwo)
    indicesPosLocLevel <- setdiff(1:m, indicesLocLevelEqZero) # not all indices 
    # from indicesLocLevelEqZero are mathematically necessary to exclude,
    # but the exclusion of all indices from this vector contributes to
    # numerical stability
    
    if(all(equation_L_shifted[indicesPosLocLevel] > 0)){
      precApprox <- TRUE
      break
    }
    counter <- counter + 1
  }
  return(precApprox)
  
}
###############################################################################
#' Checking Precision of Approximations
#'
#' @description The function checks whether information about the precision of 
#' an approximation for the informative lower SCI-bounds can be collected.
#'
#' @details
#' The function checks if it can be determined whether \code{L} can be shifted
#' by a rescaled randomly generated direction or by a rescaled direction in the
#' shift matrix such that it describes valid lower informative SCI bounds.
#' If this is possible, the approximation \code{L} is imprecise.
#' (The random directions are generated in such a way that all entries are 
#' positive.)
#' 
#' @param L An m-dimensional non-negative vector whose entries are the lower
#'  bounds of an approximation of the informative SCI.
#' @param randomShifts A positive integer indicating how many random directions
#' of length \code{tolTrueSCI} should be generated. 
#' @param shifts A matrix with m columns and any number of rows. Each entry must 
#' be non-negative. Each row is a direction in the m-dimensional real space.
#' Each row must have at least one positive entry.  \code{randomShifts} should 
#' be a positive integer or \code{shifts} should contain at least one row.
#' @param tolTrueSCI The randomly generated shift-vectors and the row vectors
#' in \code{shifts} are rescaled to have length \code{tolTrueSCI}.
#' @param gMCP An object of class \code{\link[gMCP]{graphMCP}} indicating the
#' underlying graphical test.
#' @param g  Numeric square matrix of transition weights for the graphical test
#' with m rows and m columns. The i-th row of the entered matrix defines the
#' arrows starting from the i-th hypothesis. Each entry has to be between
#' 0 and 1 and each row must sum to a number less than or equal to 1. The 
#' diagonal elements must be zero. Entering \code{g} and \code{weights} can be
#' used as the input as an alternative to specifying \code{gMCP}. 
#' @param weights Numeric vector of weights of dimension m. It defines the 
#' initial proportion of significance level which is assigned to each null 
#' hypothesis. Entering \code{g} and \code{weights} can be used as the input as
#' an alternative to specifying \code{gMCP}.
#' @param q  A numeric vector of dimension 1 or m defining the information 
#' weights for each hypothesis. The entries have to be between 0 and 1
#' (inclusive). If \code{q} has dimension 1, the same information weight is 
#' used for each hypothesis.
#' @param estimates An m-dimensional numeric vector of unbiased point estimates
#' for the parameters of interest. Each estimator is assumed to be
#' (asymptotically) normal.
#' @param Z An m-dimensional numeric vector of z-scores for testing the null 
#' hypotheses. \code{Z} can be entered as an alternative to \code{estimates}.
#' @param pValues An m-dimensional numeric vector of p-values from (asymptotic)
#' z-tests for testing the null hypotheses. \code{pValues} can be entered as an 
#' alternative to \code{estimates} and \code{Z}.
#' @param SE A positive numeric vector of dimension 1 or m indicating the 
#' standard errors of the point estimators. If \code{SE} has dimension 1, the
#' same value is used for each estimator.
#' @param I A positive numeric vector indicating the information of the
#' estimators. It can be entered as an alternative to the vector \code{SE}.
#' The relationship \code{SE}\eqn{=1/}\code{I}\eqn{^{1/2}} is assumed.  
#' @param mu_0 A numeric vector of dimension 1 or m defining the bounds of the
#' null hypotheses of the underlying graphical test. If \code{mu_0} has 
#' dimension 1, the same value is used for each null hypothesis.
#' @param alpha A numeric defining the overall significance level for the
#' graphical test (i.e. SCIs will have coverage probability of at least 
#' \eqn{1-}\code{alpha}. The parameter must be strictly between 0 and 1.
#' @param checkInput A boolean specifying whether the entered values should be 
#' checked.
#'
#' @import mvtnorm
#' @importFrom stats qnorm
#' @return Returns \code{TRUE} if we can determine that the approximation is
#' imprecise. Returns \code{FALSE} if we cannot determine that the 
#' approximation is imprecise. (The approximation may still be imprecise.) Note
#' that \code{inExactSCI} and \code{notInExactSCI} could both return 
#' \code{FALSE}.
#' @export
#' @seealso \code{\link{informSCI}} \code{\link{explore_q}}
#' @examples 
#' g <- matrix(c(0,0,1,0),2,2)
#' weights <- c(1,0)
#' q <- c(0.0068,1)
#' mu_0 <- c(-1,0)
#' pValues <- c(0.0002,0.01)
#' SE <- c(0.31,1.11)
#' alpha <- 0.025
#' L <- informSCI(g=g, weights=weights, q=q, mu_0=mu_0, pValues=pValues, SE=SE,
#' alpha=alpha, eps=1/10, tolBisec=1/10)$L
#' # When the randomShifts- or shift-parameter in the informSCI-function is
#' # specified, the notInExactSCI-function is called by the informSCI-function.
#' # It is also possible to analyse the accuracy of a calculated L (or an 
#' # approximation of the lower informative SCI-bounds) by directly using 
#' # the notInExactSCI-function:
#' notInExactSCI(L=L, randomShifts=100, tolTrueSCI=1/10^5, g=g, weights=weights, 
#' q=q, pValues=pValues, SE=SE, mu_0=mu_0, alpha=alpha)
notInExactSCI <- function(L, randomShifts = 0, shifts = NULL,
                          tolTrueSCI, gMCP = NULL, g = NULL, weights = NULL, q,
                          estimates = NULL, Z = NULL, pValues = NULL, SE = NULL,
                          I = NULL, mu_0, alpha, checkInput = TRUE){
  
  if((randomShifts <= 0) && is.null(shifts)){
    return(FALSE)
  }
  # check and process parameters: 
  checkPrecisionTestInput(L, randomShifts, shifts, tolTrueSCI, gMCP, g, weights,
                          q, estimates, Z, pValues, SE, I, mu_0, alpha)
  
  if(!(is.null(gMCP))){
    g <- getMatrix(gMCP)
    weights <- getWeights(gMCP) 
  }
  
  if(!is.null(I) && is.null(SE)){
    SE <- 1/sqrt(I)
  }
  if(!is.null(Z) && is.null(estimates)){
    estimates <- Z*SE + mu_0
  }
  if(!is.null(pValues) && is.null(estimates) && is.null(Z)){
    estimates <- qnorm(1-pValues) * SE + mu_0 
  }
  
  m <- length(weights) # number of hypotheses 
  
  # Check whether the same q, mu_0 and/ or SE should be used for all hypotheses:
  if(length(q) == 1){
    q <- rep(q,m)
  }
  if(length(mu_0) == 1){
    mu_0 <- rep(mu_0,m)
  }
  if(length(SE) == 1){
    SE <- rep(SE,m)
  }
  
  if(randomShifts > 0){
    randShifts <- abs(rmvnorm(n = randomShifts, mean = numeric(m)))
    while(any(rowSums(randShifts)==0)){
      randShifts <- abs(rmvnorm(n = randomShifts, mean = numeric(m)))  
    }
    
    shifts <- rbind(shifts,randShifts)  
  } 
  
  if(is.null(dim(shifts))){ # needed when shifts is only a vector (only one
    # direction is considered
    shifts <- matrix(shifts, nrow =  1, ncol = m) 
  }
  
  g.rowsums <- rowSums(g)
  numShifts <- dim(shifts)[1]
  
  # rescaling of shifts:
  for(i in 1:numShifts){
    shifts.norm <- norm(shifts[i, ], type = '2')
    shifts[i, ] <- shifts[i, ] * tolTrueSCI/shifts.norm
  }
  
  counter <- 1
  nonPrecApprox <- FALSE
  
  while(counter <= numShifts && nonPrecApprox == FALSE){
    
    L_shifted <- L + shifts[counter, ]
    eta.L_shifted <- weightsGTP(L_shifted, g, weights, alpha, q, mu_0)
    
    partEq <- partOfKeyEquation(L_shifted, mu_0, eta.L_shifted, q, g.rowsums)
    
    equation_L_shifted <- pvalueRTzTest(estimates, L_shifted, SE)-partEq
    
    if(all(equation_L_shifted <= 0) && all(partEq > 0)){
      nonPrecApprox <- TRUE
      break
    }
    counter <- counter + 1
  }
  return(nonPrecApprox)
}
###############################################################################
###############################################################################
# function for creating the return-list of the informSCI-function
makeReturnList <- function(out_SCIwBisec, mu_0, exactApprox, noExactApprox,
                           original_gMCPcalculated, original_gMCP){
  
  
  L <- out_SCIwBisec$L
  names(L) <- sprintf("L%d", 1:length(L)) 
  rejecHyp <- (L >= mu_0) 
  names(rejecHyp) <- sprintf("H%d", 1:length(L)) 
  returnList <- list("L" = L,
                     "rejecHyp" = rejecHyp,
                     "diffApprox" = out_SCIwBisec$diffAppr, 
                     "timesApprSmallerEps" =
                       out_SCIwBisec$counterConvergence,
                     "numIter" = out_SCIwBisec$numIter)
  if(exactApprox){
    returnList$accuracyL = "Precise approximation!"  
  }else if(noExactApprox){ 
    returnList$accuracyL = "No precise approximation!" 
  }else{
    returnList$accuracyL = "No further information could be determined."  
  }
  if(original_gMCPcalculated){
    returnList$informSCIcompatible <- 
      all(rejecHyp==original_gMCP$rejecHyp)
    returnList$cSCI = original_gMCP$L
    returnList$rejecHypGraph = original_gMCP$rejecHyp
  }
  return(returnList)
}
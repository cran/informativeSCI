# Implementation of simulation algorithm to determine optimal
# information weights

#' Exploration of the Information Weights
#' 
#' @description
#' The function calculates various statistical quantities giving some 
#' information about the behavior of informative lower SCI-bounds 
#' (\link{informSCI}) and its induced test for a given graphical test procedure
#' with m hypotheses. The simulation is done for different information weights
#' of the hypotheses. These statistical quantities are intended to be used for
#' determining information weights that represent the best possible trade-off
#' between the number of rejections and the expected size of the informative 
#' lower informative SCI-bounds. The statistical quantities can also be
#' calculated for the graphical test and the related compatible lower 
#' SCI-bounds, which allows a comparison between the two strategies. 
#' 
#'
#' @details
#' It is assumed that there are m parameters of interest
#' \eqn{\vartheta_1,\dots,\vartheta_m}. For each parameter there is a null
#' hypothesis defined as \eqn{H_i^{{\mu_0}_i}:\vartheta_i\leq{\mu_0}_i}.
#' The bounds \eqn{{\mu_0}} correspond to \code{mu_0}. The underlying graphical
#' test (specified by \code{gMCP} or \code{g} and \code{weights}) is based on 
#' these hypotheses.
#' 
#' The function simulates estimations of point estimators for the parameter of
#' interest \eqn{\vartheta_1,\dots, \vartheta_m}. The estimators follow a 
#' multivariate normal distribution with mean \code{trueParam} and covariance 
#' matrix \code{sigma}. The function repeatedly calls the 
#' \code{\link{informSCI}}-function.
#'
#' The algorithm only optimizes for a single parameter, which is used for all
#' non-fixed information weights. 
#' The parameter is chosen from a grid specified by \code{qInterval} and
#' \code{qStepsize}. The constructed grid contains all values which are between
#' \code{qInterval[1]} and \code{qInterval[2]} and can be written as 
#' \code{qInterval[1]}\eqn{+k\cdot}\code{qStepsize} where k is a natural number.
#' Alternatively, the parameter is chosen directly from \code{qGrid}.
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
#' @param trueParam A numeric vector of dimension m defining the assumed true
#' parameters \eqn{\vartheta_i, 1\leq i\leq m}.
#' @param sigma A covariance matrix of dimension \eqn{m\times m}. \code{sigma} 
#' indicates the covariance matrix of the point estimators for the parameter
#' of interest. Can be missing in the case of a many-to-one comparison.
#' Then, \code{sampleSizes}, \code{sampleSizeControl} and \code{varObs}
#' must be specified. 
#' @param qFixed A numeric matrix with l rows and 2 columns, where l is an
#' integer between 0 and m. The matrix describes the fixed information weights
#' of the simulation. The first column indicates the indices of the hypothesis
#' for which the information weight should be fixed during the simulation
#' (i.e. the entries of the first column must be natural numbers between
#' 1 and m). The second column contains the fixed values of their respective
#' fixed information weights (i.e. the entries of the second column must be
#' between 0 and 1 (inclusive)). It is permissible for all information weights 
#' to be fixed  (i.e. \code{qFixed} has m rows) or none to be fixed
#' (i.e. \code{qFixed} has 0 rows).
#' @param mu_0 A numeric vector of dimension 1 or m defining the bounds of the
#' null hypotheses of the underlying graphical test. If \code{mu_0} has 
#' dimension 1, the same value is used for each null hypothesis.
#' @param alpha A numeric defining the overall significance level for the
#' graphical test (i.e. SCIs will have coverage probability of at least 
#' \eqn{1-}\code{alpha}. The parameter must be strictly between 0 and 1.
#' @param addHyp A numeric matrix with k rows and 3 columns (k can be 0)
#' The matrix indicates for which (further) shifted hypotheses the rejection
#' probability is to be calculated. Every row describes one hypothesis.
#' The first entry is a natural number greater than m identifying the 
#' hypothesis. The second entry of each row is the index of the corresponding
#' parameter of interest. The third entry is the right border of the hypothesis.
#' @param allRej A list of vectors. Each vector in the list contains the indices
#' of subfamilies of the family of all hypotheses, including the \code{addHyp}.
#' The indices of the null hypotheses of the underlying graph range from 1 to m.
#' The indices for \code{addHyp} are given by the first column of \code{addHyp}.
#' For each such family, the probability of rejecting all hypotheses at the same
#' time is calculated.
#' @param atLeastOneRej A list of vectors. Each vector in the list contains the
#' indices of subfamilies of the family of all hypotheses, including 
#' the \code{addHyp}. The indices of the null hypotheses of the underlying graph
#' range from 1 to m. The indices for \code{addHyp} are given by the first 
#' column of \code{addHyp}. For each such family, the probability of rejecting
#' at least one hypothesis is calculated.
#' @param qGrid A numeric vector indicating the values of the non-fixed
#' information weights for the simulation. The entries must be between 0 and 1
#' (inclusive). 
#' @param qInterval A numeric vector of dimension 2 specifying the minimum
#' and maximum values allowed for the varying information weights.
#' \code{qInterval} and \code{qStepsize} can be used as the input as an
#' alternative to specifying \code{qGrid}. If all are entered, \code{qGrid}
#' is used and \code{qInterval} and \code{qStepSize} are ignored.  
#' @param qStepSize  A positive numeric defining the step size for the varying
#' information weights. \code{qInterval} and \code{qStepsize} can be used as 
#' the input as an alternative to specifying \code{qGrid}.
#' @param numSim A natural number indicating how many simulations are to be
#' performed.
#' @param sampleSizes A numeric vector indicating the sample size of each 
#' non-control group, in the many-to-one case. Not required if \code{sigma}
#' is entered. 
#' @param sampleSizeControl A numeric indicating the sample size of the control
#' group, in the many-to-one case. Not required if \code{sigma} is entered.
#' @param varObs A positive numeric indicating the variance of the individual
#' observations, in the many-to-one case. Not required if \code{sigma} is 
#' entered.
#' @param exploreGraph A boolean indicating whether the simulation should be
#' also done for the underlying graphical test and the corresponding compatible
#' lower SCI-bounds.
#' @param eps A numeric for the \code{informSCI}-algorithm indicating the 
#' desired strict upper bound on the Chebyshev distance between two successive 
#' calculated approximations (the Chebyshev distance is induced by 
#' the maximum norm).
#' @param timesSmallerEps A positive integer for the \code{informSCI}-algorithm 
#' indicating how many times the Chebyshev distance of two successive
#' calculated approximations should be less than \code{eps} in succession. Here
#' we use the convention \eqn{-\infty- (-\infty):=0}.
#' @param maxIterSCI Maximum number of iterations for determining the lower
#' informative SCI-bounds.
#' @param maxIterBisec Maximum number of iterations of the bisection method
#' which is used during the \code{informSCI}-algorithm for finding roots.
#' @param tolBisec A non-negative numeric indicating the error tolerance of
#' the bisection method which is used for finding roots in the 
#' \code{informSCI}-algorithm.
#'
#' @import mvtnorm
#' @return The function returns a list containing several statistical quantities
#' to use for the informative lower SCI-bounds to find the best possible
#' trade-off between the number of rejections and the expected size of the 
#' informative lower SCI-bounds. In the case that  \code{exploreGraph=TRUE},
#' the returned list also contains the same quantities for the (original)
#' graphical test and related compatible bounds. This allows a comparison.
#' \itemize{
#' \item \code{rejecHyp}: A matrix containing for several hypotheses the 
#' empirical rejection probability by the informative confidence bounds. 
#' The first m rows correspond to the hypotheses of the graph. The other rows
#' correspond to the hypotheses specified by \code{addHyp}. Each row indicates
#' the rejection probability for different values of the information weights.
#' \item \code{meanISCI}: A matrix containing in its columns the empirical mean
#' of the lower informative confidence bounds for different information weights.
#' Only the lower bounds which are greater than \code{-Inf} are used for the
#' empirical mean. 
#' \item \code{impISCI}: A matrix containing in its columns the empirical 
#' average distance between the lower informative confidence bounds and 
#' \code{mu_0} for different information weights. Only the lower bounds which
#' are greater than \code{-Inf} are used for the empirical average distance.
#' \item \code{biasISCI}: A matrix containing in its columns the empirical
#' average distance between the lower informative confidence bounds and the
#' true parameters \code{trueParam} for different information weights. Only the
#' lower bounds which are greater than \code{-Inf} are used for the empirical
#' average distance.
#' \item \code{numISCIfinite}: A matrix containing in its columns how many times
#' the lower informative confidence bounds were each greater than \code{-Inf} 
#' for different information weights.
#' \item \code{rejecAllHyp}: A matrix containing in its columns for each family
#' from \code{allRej} the empirical probability of rejecting all of the 
#' hypotheses from the family with the induced test at the same time for 
#' different information weights. 
#' \item \code{rejecAtLeastHyp}: A matrix containing in its columns for each 
#' family from \code{atLeastOneRej} the empirical probability of rejecting
#' at least one of the hypotheses from the family with the induced test for 
#' different information weights.
#' }
#' 
#' If \code{exploreGraph=TRUE}:
#' \itemize{
#' \item \code{rejecHypGraph}: A vector containing for each of the null
#' hypotheses of the graph and of the additional hypotheses (specified by
#' \code{addHyp}) its empirical rejection probability by the original graph.
#' \item \code{meanCSCI}: A vector containing, for each parameter 
#' \eqn{\vartheta_i, 1\leq i\leq m} the empirical mean of the lower compatible
#' confidence bounds. Only the lower bounds which are greater than \code{-Inf} 
#' are used for the empirical mean. 
#' \item \code{impCSCI}: A vector containing, for each parameter, the empirical
#' average distance between the lower compatible confidence bounds and
#' \code{mu_0}. Only the lower bounds which are greater than \code{-Inf} are
#' used.
#' \item \code{biasCSCI}: A vector containing, for each parameter,
#' the empirical average distance between the lower compatible confidence bounds
#' and the true parameters \code{trueParam}. Only the lower bounds which are 
#' greater than \code{-Inf} are used.
#' \item \code{numCSCIfinite}: A vector containing, for each parameter, how 
#' many times the compatible lower confidence bounds were each greater 
#' than \code{-Inf}.
#' \item \code{rejecAllHypCSCI}: A vector containing, for each family from 
#' \code{allRej}, the empirical probability of rejecting all of the hypotheses
#' from the family with the (original) graphical test.
#' \item \code{rejecAtLeastHypCSCI}: A vector containing, for each family from
#' \code{atLeastOneRej}, the empirical probability of rejecting at least one
#' of the hypotheses from the family with the (original) graphical test.
#' }
#' @references 
#' S. Schmidt, W. Brannath: Informative simultaneous confidence intervals
#' for the fallback procedure. Biometrical Journal 57.4 (2015), pp. 712–719.
#' @export
#'
#' @seealso \code{\link{informSCI}} \code{\link[gMCP]{gMCP}}
#'  \code{\link[gMCP]{simConfint}} 
#' 
#' @examples 
#' explore_q(gMCP=BonferroniHolm(3), trueParam=c(1.5,1,0.2),
#' sigma=diag(3)*0.2, qFixed=matrix(c(2,3,0.3,0.3),2,2), mu_0=c(-0.5,0,0),
#' addHyp=matrix(c(4,1,0),1,3),allRej =list(c(1,2), c(4,2)), 
#' atLeastOneRej=list(c(2,3)),numSim=100)
#' @examples
#' explore_q(g=matrix(c(0,0,1,0),2,2), weights=c(1,0), trueParam=c(0.5,2), 
#' mu_0=c(-1,0), alpha=0.025, qGrid=c(1/10*c(1:10),c(0.97,0.98,0.99)), 
#' numSim=100, sampleSizes=c(89,95), sampleSizeControl=77, varObs=10)
#' 

explore_q <- function(gMCP = NULL, g = NULL, weights = NULL, trueParam,
                      sigma = NULL, qFixed = matrix(0,0,2), mu_0 = 0,
                      alpha = 0.05, addHyp =  matrix(0,0,3), allRej = NULL,
                      atLeastOneRej = NULL, qGrid = NULL, qInterval = c(0,1), 
                      qStepSize = 1/10,  numSim = 1000, sampleSizes = NULL,
                      sampleSizeControl = NULL, varObs = NULL,
                      exploreGraph = TRUE, eps = 1/10^5, timesSmallerEps = 3,
                      maxIterSCI = 1000, maxIterBisec = 1000,
                      tolBisec = 1/10^3){ 
  
  ########## check and process parameters: ###########
  
  checkSimulationValues(gMCP, g, weights, trueParam, sigma, qFixed, mu_0, alpha,
                        addHyp, allRej, atLeastOneRej, qInterval, qStepSize,
                        qGrid, numSim,  sampleSizes, sampleSizeControl, varObs,
                        timesSmallerEps, maxIterSCI, maxIterBisec)
  
  
  if(!(is.null(gMCP))){
    g <- getMatrix(gMCP)
    weights <- getWeights(gMCP) 
  }
  
  m <-length(weights) # number of hypotheses 
  
  # Check whether the same mu_0 should be used for all hypotheses:
  if(length(mu_0) == 1){
    mu_0 <- rep(mu_0,m)
  } 
  
  # Check whether the simulation should be done for a many-to-one-comparison.
  # In this case, sigma is calculated:
  if(is.null(sigma)){
    sigma <- sigmaManyToOne(sampleSizes, sampleSizeControl, varObs,
                            checkInput = FALSE)
  }
  
  if(is.null(allRej)){
    allRej <- list(c(1:m))  
  }
  if(is.null(atLeastOneRej)){
    atLeastOneRej <- list(c(1:m))  
  }
  
  
  ######## set the values of q for which a simulation will be performed: #######
  
  if(dim(qFixed)[1] == m){
    numGridPoints <- 1  # If all entries of q are fixed, a simulation will
    # only be performed for exactly this q.
    qGrid <- NaN # The information weights will not be varied because all
    # entries are fixed
  }else{
    if(is.null(qGrid)){
      a <- qInterval[1]
      b <- qInterval[2]  
      numGridPoints <- floor((b-a)/qStepSize) + 1 
      qGrid <- (a + c(0:(numGridPoints-1)) * qStepSize)
      
    }else{
      numGridPoints <- length(qGrid)
    }  
  }
  
  
  ########### variables are defined for the output: ##########
  
  
  # variables for informative bounds and the induced multiple test procedure:
  
  # matrix for rejection probability of the hypotheses of the original graph
  # and of the additional hypotheses (characterized by addHyp):
  numAddHyp <- nrow(addHyp)
  rejecPropHyp <- matrix(0, m + numAddHyp, numGridPoints) 
  
  # matrix for mean of lower informative SCI bounds:
  informB.mean <- matrix(0, m, numGridPoints)
  
  # matrix for number of simulations for which the informative bounds are
  # greater than minus infinity:
  informB.finite <- matrix(0, m, numGridPoints)
  
  # for each familiy of hypotheses from the allRej-list, the following matrix 
  # should contain the probability of rejecting all of these hypotheses with the
  # induced test at the same time:
  numAllRej <- length(allRej) 
  rejecProbAllHyp <- matrix(0, numAllRej, numGridPoints)
  
  # for each family of hypotheses from the atLeastOneRej-list, the following
  # matrix should contain the probability of rejecting at least one of these
  # hypotheses with the induced test:
  numAtLeOneRej <- length(atLeastOneRej) 
  rejecProbAtLeastHyp <- matrix(0, numAtLeOneRej, numGridPoints) 
  
  #############################################################################
  # for comparison with original graphical test procedure and compatible lower
  # SCI bounds:
  
  # array for rejection probability of the hypotheses of the original graph
  # and of the additional hypotheses (characterized by addHyp):
  # procedure:
  rejecPropHypCSCI <- numeric(m + numAddHyp)
  
  # array for mean of lower compatible SCI bounds:
  compB.mean <- numeric(m)
  
  # array for number of simulations for which the compatible bounds are greater
  # than minus infinity:
  compB.finite <- numeric(m)
  
  # for each family of hypotheses from the allRej-list, the following vector
  # should contain the probability of rejecting all of these hypotheses by the
  # original graph at the same time:
  rejecProbAllHypCSCI <- numeric(numAllRej)
  
  # for each family of hypotheses from the atLeastOneRej-list, the following
  # vector should contain the probability of rejecting at least of these
  # hypotheses by the original graph:
  rejecProbAtLeastHypCSCI <- numeric(numAtLeOneRej) 
  
  ############ simulation: ############
  
  means_simulations <- rmvnorm(n = numSim, mean = trueParam, sigma = sigma)
  # perform simulations of the estimators
  
  SE <- sqrt(diag(sigma))
  
  for(l in 1:numSim){
    for(j in 1:numGridPoints){
      q <- qGrid[j] * rep(1,m)
      q[qFixed[ ,1]] <- qFixed[ ,2] # the simulation is only performed for
      # the non-fixed information weights.
      
      sci <- informSCI(g = g, weights = weights, q = q, mu_0 = mu_0,
                       estimates = means_simulations[l, ], SE = SE, 
                       alpha = alpha, eps = eps,
                       timesSmallerEps = timesSmallerEps,
                       maxIter = maxIterSCI, maxIterBisec = maxIterBisec,
                       tolBisec = tolBisec, calculateCSCI = FALSE,
                       checkInput = FALSE) 
      
      informB <- sci$L
      
      # create desired output values related to the informative SCI and the
      # associated induced test:
      
      
      rejecPropHyp[1:m,j] <- rejecPropHyp[1:m,j] + as.numeric(informB >= mu_0)
      
      if(numAddHyp >= 1){
        rejecPropHyp[(m+1):(m+numAddHyp),j] <- 
          rejecPropHyp[(m+1):(m+numAddHyp),j] +
          as.numeric(informB[addHyp[ ,2]] >= addHyp[ ,3])  
      }
      
      
      indicesRealB <- which(informB > -Inf)
      informB.mean[indicesRealB,j] <- informB.mean[indicesRealB,j] +
        informB[indicesRealB]
      
      informB.finite[indicesRealB,j]  <-
        informB.finite[indicesRealB,j] + 1
      
      rejectedHyp <- c(informB >= mu_0, informB[addHyp[ ,2]] >= addHyp[ ,3])
      
      if(numAllRej >= 1){
        # check whether for each family in allRej all hypotheses of this family
        # are rejected by the induced test:
        rejAll <- sapply(allRej,
                         function(set){as.numeric(all(rejectedHyp[set]))}) 
        rejecProbAllHyp[ ,j] <- rejecProbAllHyp[ ,j] + array(unlist(rejAll))
      }
      if(numAtLeOneRej >= 1){
        # check whether for each family in atLeastOneRej at least one
        # hypothesis of this family is rejected by the induced test:
        rejAtL <- sapply(atLeastOneRej,
                         function(set){as.numeric(any(rejectedHyp[set]))})  
        rejecProbAtLeastHyp[ ,j] <- rejecProbAtLeastHyp[ ,j] +
          array(unlist(rejAtL)) 
      }
    }
  }
  
  
  if(exploreGraph){
    for(l in 1:numSim){
      # create desired output values related to the graphical test procedure and 
      # compatible SCI:
      origiTest <- graphTestAndCompSCI(g, weights, 
                                       estimates = means_simulations[l, ], SE,
                                       mu_0, alpha)
      
      compB <- origiTest$L 
      
      rejecPropHypCSCI[1:m] <- rejecPropHypCSCI[1:m] +
        as.numeric(origiTest$rejecHyp == TRUE)
      
      if(numAddHyp >= 1){
        rejecPropHypCSCI[(m+1):(m+numAddHyp)] <-
          rejecPropHypCSCI[(m+1):(m+numAddHyp)] + 
          as.numeric(compB[addHyp[ ,2]] >= addHyp[ ,3])
      }
      
      indicesCompBReal <- which(compB > -Inf)
      compB.mean[indicesCompBReal] <- compB.mean[indicesCompBReal] +
        compB[indicesCompBReal] 
      
      compB.finite[indicesCompBReal] <- compB.finite[indicesCompBReal] + 1
      
      origiGraphRejHyp <- c(origiTest$rejecHyp,compB[addHyp[ ,2]] >= addHyp[ ,3])
      
      if(numAllRej >= 1){
        # check whether for each family in allRej all hypotheses of this set
        # are rejected by the original graph:
        rejAll <- sapply(allRej,
                         function(set){as.numeric(all(origiGraphRejHyp[set]))}) 
        
        rejecProbAllHypCSCI <- rejecProbAllHypCSCI +
          array(unlist(rejAll)) 
      }
      if(numAtLeOneRej >= 1){
        # check whether for each family in atLeastOneRej at least one
        # hypothesis of this family is rejected by the original graph:
        rejAtL <- sapply(atLeastOneRej,
                         function(set){as.numeric(any(origiGraphRejHyp[set]))})  
        rejecProbAtLeastHypCSCI <-rejecProbAtLeastHypCSCI + 
          array(unlist(rejAtL)) 
      }
      
    }
  }
  
  ######### final calculations and creation of return list: ##############
  
  # variables related to the informative bounds:
  rejecPropHyp <-rejecPropHyp/numSim 
  informB.mean <- informB.mean/informB.finite 
  rejecProbAllHyp <- rejecProbAllHyp/numSim 
  rejecProbAtLeastHyp <- rejecProbAtLeastHyp/numSim 
  
  # variables related to compatible bounds:
  rejecPropHypCSCI <- rejecPropHypCSCI/numSim 
  compB.mean <- compB.mean/compB.finite
  rejecProbAllHypCSCI <- rejecProbAllHypCSCI/numSim 
  rejecProbAtLeastHypCSCI <-rejecProbAtLeastHypCSCI/numSim 
  
  # create return list:
  listISCI <- makeListISCI(m, qGrid, rejecPropHyp, addHyp, informB.mean, mu_0,
                           trueParam, informB.finite, allRej,
                           rejecProbAllHyp, atLeastOneRej, rejecProbAtLeastHyp)
  
  if(exploreGraph){
    listCSCI <- makeListCSCI(m, rejecPropHypCSCI, addHyp, compB.mean, mu_0,
                             trueParam, compB.finite, allRej,
                             rejecProbAllHypCSCI, atLeastOneRej,
                             rejecProbAtLeastHypCSCI)
    return(c(listISCI, listCSCI))
  }else{
    return(listISCI)
  }
}
###############################################################################
########################### helper functions ##################################
###############################################################################

#' Calculation of the Covariance Matrix for a Many-to-one-Comparison
#'
#' @description
#' The function calculates the covariance matrix for many-to-one-comparisons.
#' The covariance matrix is calculated for the point estimators, 
#' each defined by the difference between the empirical mean of one of
#' the experimental groups and the empirical mean of the control group.
#' 
#' @param sampleSizes A numeric vector indicating the sample size of each
#' non-control group.
#' @param sampleSizeControl A numeric indicating the sample size of the control
#' group.
#' @param varObs A positive numeric indicating the variance of the individual
#' observations.
#' @param checkInput A boolean specifying whether the entered values should be 
#' checked.
#'
#' @return Returns covariance matrix.
#' @export
#'
#' @examples sigmaManyToOne(sampleSizes=c(89,95), sampleSizeControl=77,
#'  varObs=10)
sigmaManyToOne <- function(sampleSizes, sampleSizeControl, varObs,
                           checkInput = TRUE){
  # check input:
  if(checkInput){
    checkPartOfInput4(sampleSizes, sampleSizeControl, varObs)  
  }
  # calculate covariance matrix:
  m <- length(sampleSizes)
  covMatrix <- matrix(varObs/sampleSizeControl, m, m)
  diag(covMatrix) <- (1/sampleSizes * varObs) + diag(covMatrix)
  return(covMatrix)
}

###############################################################################

# function for creating the elements of the return-list of the 
# explore_q-function which are related to the informative SCI-bounds:

makeListISCI <- function(m, qGrid, rejecPropHyp, addHyp, informB.mean, mu_0,
                         trueParam, informB.finite, allRej,
                         rejecProbAllHyp, atLeastOneRej, rejecProbAtLeastHyp){
  
  rownames(rejecPropHyp) <- c(sprintf("P(L_%d >= mu_0[%d])", 1:m, 1:m),
                              sprintf("P(L_%d >= %d)",addHyp[ ,2], addHyp[ ,3]))
  colnames(rejecPropHyp) <- as.character(qGrid)
  
  rownames(informB.mean) <- sprintf("E(L_%d)", 1:m)
  colnames(informB.mean) <- as.character(qGrid)
  
  averageDistanceToMu_0 <- informB.mean - mu_0
  rownames(averageDistanceToMu_0) <- sprintf("E(L_%d)-mu_0[%d]", 1:m, 1:m)
  colnames(averageDistanceToMu_0) <- as.character(qGrid)
  
  averageDistanceToTrueParam <- informB.mean - trueParam 
  rownames(averageDistanceToTrueParam) <-
    sprintf("E(L_%d)-trueParam[%d]", 1:m, 1:m)
  colnames(averageDistanceToTrueParam) <- as.character(qGrid)
  
  
  rownames(informB.finite) <- sprintf("L%d:", 1:m)
  colnames(informB.finite) <- as.character(qGrid)
  
  if(length(allRej) >= 1){
    namesSets <- names(allRej)
    if(is.null(namesSets)){
      emptyNames <- c(1:length(allRej))  
    }else{
      emptyNames <- which(namesSets=="")  
    }
    for(i in emptyNames){
      namesSets[i] <- paste(as.character(allRej[[i]]),collapse=",")
    }
    rownames(rejecProbAllHyp) <- namesSets
    colnames(rejecProbAllHyp) <- as.character(qGrid)  
  }
  if(length(atLeastOneRej) >= 1){
    namesSets <- names(atLeastOneRej)
    if(is.null(namesSets)){
      emptyNames <- c(1:length(atLeastOneRej))
    }else{
      emptyNames <- which(namesSets=="")  
    }
    for(i in emptyNames){
      namesSets[i] <- paste(as.character(atLeastOneRej[[i]]),collapse=",")
    }
    rownames(rejecProbAtLeastHyp) <- namesSets
    colnames(rejecProbAtLeastHyp) <- as.character(qGrid)  
  }
  returnList <- list("rejecHyp" = rejecPropHyp, "meanISCI" = informB.mean,
                     "impISCI" = averageDistanceToMu_0,
                     "biasISCI"= averageDistanceToTrueParam,
                     "numISCIfinite" = informB.finite,
                     "rejecAllHyp"= rejecProbAllHyp,
                     "rejecAtLeastHyp" = rejecProbAtLeastHyp)
  return(returnList)
}
###############################################################################
# function for creating the elements of the return-list of the 
# explore_q-function which are related to the compatible SCI bounds:

makeListCSCI <- function(m, rejecPropHypCSCI, addHyp, compB.mean, mu_0,
                         trueParam, compB.finite, allRej,
                         rejecProbAllHypCSCI, atLeastOneRej,
                         rejecProbAtLeastHypCSCI){
  
  names(rejecPropHypCSCI) <- c(sprintf("P(L_%d >= mu_0[%d])", 1:m, 1:m),
                               sprintf("P(L_%d >= %d)", addHyp[ ,2], 
                                       addHyp[ ,3]))
  
  names(compB.mean) <- sprintf("E(L_%d)", 1:m)
  
  compB.averageDistanceToMu_0 <- compB.mean - mu_0
  names(compB.averageDistanceToMu_0) <-
    sprintf("E(L_%d)-mu_0[%d]", 1:m, 1:m)
  
  compB.averageDistanceToTrueParam <- compB.mean - trueParam
  names(compB.averageDistanceToTrueParam) <-
    sprintf("E(L_%d)-trueParam[%d]", 1:m, 1:m)
  
  names(compB.finite) <- sprintf("L%d:", 1:m)
  
  if(length(allRej) >= 1){
    namesSets <- names(allRej)
    if(is.null(namesSets)){
      emptyNames <- c(1:length(allRej))  
    }else{
      emptyNames <- which(namesSets=="")  
    }
    for(i in emptyNames){
      namesSets[i] <- paste(as.character(allRej[[i]]),collapse=",")
    }
    names(rejecProbAllHypCSCI) <- namesSets
  }
  if(length(atLeastOneRej) >= 1){
    namesSets <- names(atLeastOneRej)
    if(is.null(namesSets)){
      emptyNames <- c(1:length(atLeastOneRej))  
    }else{
      emptyNames <- which(namesSets=="")  
    }
    for(i in emptyNames){
      namesSets[i] <- paste(as.character(atLeastOneRej[[i]]),collapse=",")
    }
    names(rejecProbAtLeastHypCSCI) <- namesSets
  }
  
  returnList <- list("rejecHypGraph" = rejecPropHypCSCI,
                     "meanCSCI" = compB.mean,
                     "impCSCI" = compB.averageDistanceToMu_0,
                     "biasCSCI"= compB.averageDistanceToTrueParam,
                     "numCSCIfinite" = compB.finite,
                     "rejecAllHypCSCI" = rejecProbAllHypCSCI,
                     "rejecAtLeastHypCSCI" = rejecProbAtLeastHypCSCI)
  return(returnList)
}
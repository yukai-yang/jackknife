#################################################################################
## package name: jackknife
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
#################################################################################


#' Implementing the jackknife method of estimating standard errors
#'
#' Implementing the jackknife method of estimating standard errors of some linear statistics
#'
#' The function takes several arguments including functions and returns a list containing the results.
#'
#' @param vx a vector of observations
#' @param tfunc the function which takes the rho argument and returns the statistic
#' @param dfunc the delta function which takes a tuple or window of the observations and returns a vector of deltas
#' @param im the tuple or window size
#' @param vw the weights for deleting or downweighting the blocks
#' @param ... other arguments for both \code{tfunc} and \code{dfunc} if any
#'
#' @return a list containing the results
#'
#' The object is a list containing the following components:
#' \item{T}{the value of the statistic}
#' \item{jack}{a vector or matrix containing the jackknifed statistics}
#' \item{sig2}{the jackknifed estimate of the variance of the statistic}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#'
#' @keywords jackknife
#'
#' @examples
#' vw1 = 1
#' vw2 = c(0.25,0.75,1,0.75,0.25)
#'
#' # the AR(1) model
#' iN = 100
#' vx <- rnorm(iN)
#'
#' # for first order autocorrelation
#' # dfunc vec 2, y_t and y_t-1
#' dfunc <- function(yy, my, vy) return(prod(yy - my)/vy)
#'
#' tfunc <- function(rho, my, vy) return(rho)
#'
#' # results
#' ret = jackknife(as.vector(vx), tfunc, dfunc, im=2, vw=vw1, my=mean(vx), vy=var(vx))
#' ret
#'
#' ret = jackknife(as.vector(vx), tfunc, dfunc, im=2, vw=vw2, my=mean(vx), vy=var(vx))
#' ret
#'
#' @export
jackknife <- function(vx, tfunc, dfunc, im, vw, ...){
  iN = length(vx)
  nn = iN - im + 1
  il = length(vw)

  ftmp <- function(nter) return(dfunc(vx[nter:(nter+im-1)], ...))

  deltas = sapply(1:nn, ftmp)

  if(is.matrix(deltas)){
    deltas = t(deltas)
    rho = colMeans(deltas)
  }else rho = mean(deltas)

  # compute the statistic T = T(rho)
  TT = tfunc(rho, ...)

  # the two norms
  nw1 = sum(vw)
  nw22 = c(vw%*%vw)

  jack <- function(nter){
    weight = rep(1, nn); weight[nter:(nter+il-1)] = 1-vw
    # return T_j
    return(tfunc((weight %*% deltas)/(nn - nw1), ...))
  }

  TTT = sapply(1:(nn-il+1), jack)
  if(is.matrix(TTT)){
    JSig2 = rowMeans((TTT - rowMeans(TTT))^2) * (nn-nw1)^2 / nw22 / nn
  }else JSig2 = mean((TTT - mean(TTT))^2) * (nn-nw1)^2 / nw22 / nn

  return(list(T=TT, jack=TTT, sig2=JSig2))
}

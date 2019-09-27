# tfunc, the statistic function T which is a function of rho
# dfunc, delta function, input a m tuple and return a scalar
# im the window size or tuple size
# vw the weights for deleting or downweighting blocks with length l
# ... other arguments for tfunc and dfunc
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

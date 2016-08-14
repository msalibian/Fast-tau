
FastTau <- function(x, y, N=500, kk=2, tt=5, rr=2, approximate=0, seed=123)
{
  # fast-tau algorithm for linear regression
  #
  # tau-estimate is tuned to have 95% efficiency, and 50% bdp,
  # using Optimal rho-function
  #
  # INPUT:
  # 	y : response vector (n x 1)
  # 	x : covariates matrix (n x p), possibly including intercept column
  # 	N : number of elemental starts, e.g. 500 
  #	  kk : number of initial IRWLS steps on each candidate
  #	  tt : number of best solutions to RWLS-iterate until convergence
  #   rr : number of iterations in scale approximation in case approximate=1
  #   approximate : if 0, fully compute S-scale in objective function evaluation, 
  #               otherwise approximate 
  # OUTPUT:
  # 	res$beta : tau-estimate of regression coefficients
  #	res$scale : tau-estimate of residual scale
  
  if (tt<1) stop("parameter tt should be at least 1")
  
  
  # save existing random seed
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  
  if(!missing(seed)) set.seed(seed)
  
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  c1 <- .4046
  b1 <- .5
  c2 <- 1.09
  b2 <- .1278
  
  RWLStol <- 1e-11
  
  bestbetas <- matrix(0, p, tt)
  bestscales <- 1e20 * rep(1, tt)
  besttauscales <- 1e20 * rep(1, tt)
  worsti <- 1
  rworst <- y
  
  for (i in 1:N) {
    # find a p-subset in general position.
    singular <- 1; itertest <- 1
    while (singular==1 && itertest<100) {
      ranset <- randomset(n,p)
      xj <- x[ranset,]
      yj <- y[ranset] 
      bj <- as.matrix(qr.coef(qr(xj),yj))
      singular <- any(!is.finite(bj))
      itertest <- itertest + 1
    }
    if (itertest==100) stop("Too many degenerate subsamples")
    
    # perform kk steps of IRLS on elemental start
    if (kk > 0) {
      tmp <- IWLSiteration(x, y, bj, 0, kk, RWLStol, b1, c1, c2)
      betarw <- tmp$betarw
      resrw <- y - x %*% betarw
      scalerw <- tmp$scalerw
    }
    else {
      betarw <- bj
      resrw <- y - x %*% betarw
      scalerw <- median(abs(resrw))/.6745
    }
    
    # long-term memory vector, for finding a special extra candidate at the end :
    if (i > 1) LTMvec = LTMvec + abs(resrw)
    else LTMvec = abs(resrw)
    
    
    # check whether new subsample yields one of best t tau-objective values
    
    if (!approximate) { # compute actual scale, but use tau-conditions!
      scaletest1 <- mean(rhoOpt(resrw / bestscales[worsti],c1)) < b1
      scaletest2 <- sum(rhoOpt(resrw / bestscales[worsti],c2)) < sum(rhoOpt(rworst/bestscales[worsti],c2))
      if (scaletest1 || scaletest2) {
        # if conditions fulfulled, compute objective value
        snew <- Mscale(resrw, b1, c1, scalerw)
        taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
        if (taunew < besttauscales[worsti]) {
          # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
          besttauscales[worsti] <- taunew
          bestscales[worsti] <- snew
          bestbetas[,worsti] <- betarw
          worsti <- which.max(besttauscales) 
          rworst <- y - x %*% bestbetas[,worsti]
        }
      }
    }
    else { # or just compute approximations (and don't bother with the conditions)
      snew = scalerw;
      if (rr>0) {
        for (kstep in 1:rr) { 
          snew <- sqrt( snew^2 * mean( rhoOpt(resrw/snew,c1) ) / b1 )
        }
      }
      taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
      if (taunew < besttauscales[worsti]) {
        # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
        besttauscales[worsti] <- taunew
        bestscales[worsti] <- snew
        bestbetas[,worsti] <- betarw
        worsti <- which.max(besttauscales) 
        rworst <- y - x %*% bestbetas[,worsti]
      }
    }
  }
  
  # consider an extra subsample, made up of badly fit observations
  
  IXLTM <- order(LTMvec, decreasing=T)
  singular <- 1 
  extrasize <- p
  while (singular==1) {
    xs <- x[IXLTM[1:extrasize],]
    ys <- y[IXLTM[1:extrasize]]
    bbeta <- as.matrix(qr.coef(qr(xs),ys))
    singular <- any(!is.finite(bbeta))
    extrasize <- extrasize + 1
  }
  
  # perform kk steps of IRLS on elemental start
  if (kk > 0) {
    tmp <- IWLSiteration(x, y, bbeta, 0, kk, RWLStol, b1, c1, c2)
    betarw <- tmp$betarw
    resrw <- y - x %*% betarw
    scalerw <- tmp$scalerw
  }
  else {
    betarw <- bbeta
    resrw <- y - x %*% betarw
    scalerw <- median(abs(resrw))/.6745
  }
  
  # check whether this candidate yields one of best t tau-objective values
  
  if (!approximate) { # compute actual scale, but use tau-conditions!
    scaletest1 <- mean(rhoOpt(resrw / bestscales[worsti],c1)) < b1
    scaletest2 <- sum(rhoOpt(resrw / bestscales[worsti],c2)) < sum(rhoOpt(rworst/bestscales[worsti],c2))
    if (scaletest1 || scaletest2) {
      # if conditions fulfulled, compute objective value
      snew <- Mscale(resrw, b1, c1, scalerw)
      taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
      if (taunew < besttauscales[worsti]) {
        # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
        besttauscales[worsti] <- taunew
        bestscales[worsti] <- snew
        bestbetas[,worsti] <- betarw
        worsti <- which.max(besttauscales) 
        rworst <- y - x %*% bestbetas[,worsti]
      }
    }
  }
  else { # or just compute approximations (and don't bother with the conditions)
    snew = scalerw;
    if (rr>0) {
      for (kstep in 1:rr) { 
        snew <- sqrt( snew^2 * mean( rhoOpt(resrw/snew,c1) ) / b1 )
      }
    }  
    taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
    if (taunew < besttauscales[worsti]) {
      # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
      besttauscales[worsti] <- taunew
      bestscales[worsti] <- snew
      bestbetas[,worsti] <- betarw
      worsti <- which.max(besttauscales) 
      rworst <- y - x %*% bestbetas[,worsti]
    }
  }
  
  superbesttauscale <- 1e20
  
  # RWLS-iterate each of the best tt candidates until convergence, and retain the best result
  for (i in 1:tt) {
    tmp <- IWLSiteration(x, y, bestbetas[,i], bestscales[i], 500, RWLStol, b1, c1, c2)
    resrw <- y - x %*% tmp$betarw
    tauscalerw <- tmp$scalerw * sqrt(mean(rhoOpt(resrw/tmp$scalerw,c2)))
    if (tauscalerw < superbesttauscale) {
      superbesttauscale <- tauscalerw
      superbestbeta <- tmp$betarw
      superbestscale <- tmp$scalerw
    }
  }
  
  superbestscale <- Mscale(y - x%*%superbestbeta, b1, c1, superbestscale)
  superbesttauscale <- superbestscale * sqrt(mean(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2)))
  
  # Olive's two extra candidates:
  
  # add LS candidate
  betaLS <- as.matrix(qr.coef(qr(x),y))
  resLS <- y - x %*% betaLS
  scaleLS <- median(abs(resLS))/.6745  
  scaletest1 <- mean(rhoOpt(resLS / superbestscale,c1)) < b1
  scaletest2 <- sum(rhoOpt(resLS / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
  if (scaletest1 || scaletest2) {
    snew <- Mscale(resLS, b1, c1, scaleLS)
    taunew <- snew * sqrt(mean(rhoOpt(resLS/snew,c2)))
    if (taunew < superbesttauscale) {
      superbestscale <- snew
      superbestbeta <- betaLS
      superbesttauscale <- taunew
    }   
  }
  
  # add HB candidate
  IXmed <- order(abs(y - median(y)))
  xhalf <- x[IXmed[1:floor(n/2)],]
  yhalf <- y[IXmed[1:floor(n/2)]]
  bbeta <- as.matrix(qr.coef(qr(xhalf),yhalf))
  # + 10 C-steps
  tmp <- IWLSiteration(x, y, bbeta, 0, 10, RWLStol, b1, c1, c2)
  betaHB <- tmp$betarw
  resHB <- y - x %*% betaHB
  scaleHB <- tmp$scalerw
  scaletest1 <- mean(rhoOpt(resHB / superbestscale,c1)) < b1
  scaletest2 <- sum(rhoOpt(resHB / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
  if (scaletest1 || scaletest2) {
    snew <- Mscale(resHB, b1, c1, scaleHB)
    taunew <- snew * sqrt(mean(rhoOpt(resHB/snew,c2)))
    if (taunew < superbesttauscale) {
      superbestbeta <- betaHB
      superbesttauscale <- taunew
    }   
  }
  
  # restore seed existing before call
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  
  return(list( beta = superbestbeta, scale = superbesttauscale/sqrt(b2) ))
  
}


# -------------------------------------------------------------------

IWLSiteration <- function(x, y, inib, iniscale, maxiter, tol, b, c1, c2)
{  
  # approximate IRWLS iteration; pass maxiter=500, say, if convergence is desired
  # e.g. tol = 1e-11
  
  n <- nrow(x)
  p <- ncol(x)
  
  res <- y - x %*% inib
  if (iniscale == 0)
    scale <- median(abs(res))/.6745
  else
    scale <- iniscale
  
  oldbeta <- inib
  
  betadiff <- 2*tol
  iter <- 0
  while ((betadiff > tol) && (iter < maxiter)) {
    scale <- sqrt( scale^2 * mean( rhoOpt(res/scale,c1) ) / b )
    scaledres <- res/scale
    Wn.teller <- sum(WtellerOpt(scaledres,c2))
    Wn.noemer <- sum(psixOpt(scaledres,c1))
    Wn <- Wn.teller / Wn.noemer
    weights <- (Wn * fwOpt(scaledres,c1) + fwOpt(scaledres,c2))
    sqweights <- sqrt(weights)
    xw <- x * as.vector(sqweights)
    yw <- y * sqweights
    newbeta <- qr.coef(qr(xw),yw)
    if (any(!is.finite(newbeta))) {
      newbeta <- inib
      scale <- iniscale
      break
    }
    betadiff <- sqrt(sum((oldbeta - newbeta)^2))
    res <- y - x %*% newbeta
    oldbeta <- newbeta
    iter <- iter + 1
  }
  
  return( list( betarw = newbeta, scalerw = scale ) )
  
}

#--------------------------------------------------------------------------  

Mscale <- function(u, b, c, initialsc) 
{
  # from Kristel's fastSreg
  if (initialsc==0)
    initialsc = median(abs(u))/.6745
  maxit <- 100
  sc <- initialsc
  i <- 0 
  eps <- 1e-10
  err <- 1
  while  (( i < maxit ) & (err > eps)) {
    sc2 <- sqrt( sc^2 * mean(rhoOpt(u/sc,c)) / b)
    err <- abs(sc2/sc - 1)
    sc <- sc2
    i <- i+1
  }
  
  return(sc)
  
}

#---------------------------------------------------------------------------------------

randomset <- function(tot,nel) {
  ranset <- rep(0,nel)
  for (j in 1:nel) {
    num <- ceiling(runif(1)*tot)
    if (j > 1) {
      while (any(ranset==num)) 
        num <- ceiling(runif(1)*tot)
    }
    ranset[j] <- num
  }
  return(ranset)
}

# --------------------------------------------------------------------

rhoOpt <- function(x, cc)
{
  tmp <- x^2 / 2 / (3.25*cc^2)
  tmp2 <- (1.792 - 0.972 * x^2 / cc^2 + 0.432 * x^4 / cc^4 - 0.052 * x^6 / cc^6 + 0.002 * x^8 / cc^8) / 3.25
  tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
  tmp[abs(x) > 3*cc] <- 1
  tmp
  
}

# --------------------------------------------------------------------

fwOpt <- function(x, cc)
{
  tmp <- (-1.944 / cc^2 + 1.728 * x^2 / cc^4 - 0.312 * x^4 / cc^6 + 0.016 * x^6 / cc^8) / 3.25
  tmp[abs(x) < 2*cc] <- 1 / (3.25*cc^2)
  tmp[abs(x) > 3*cc] <- 0
  tmp
  
}

# --------------------------------------------------------------------

psiOpt <- function(x, cc)
{
  tmp <- x / (3.25*cc^2)
  tmp2 <- (-1.944 * x / cc^2 + 1.728 * x^3 / cc^4 - 0.312 * x^5 / cc^6 + 0.016 * x^7 / cc^8) / 3.25
  tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
  tmp[abs(x) > 3*cc] <- 0
  tmp
  
}

# --------------------------------------------------------------------

psixOpt <- function(x, cc)
{
  tmp <- x^2 / (3.25*cc^2)
  tmp2 <- (-1.944 * x^2 / cc^2 + 1.728 * x^4 / cc^4 - 0.312 * x^6 / cc^6 + 0.016 * x^8 / cc^8) / 3.25
  tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
  tmp[abs(x) > 3*cc] <- 0
  tmp
  
}
# --------------------------------------------------------------------

WtellerOpt <- function(x, cc)
{
  tmp <- (3.584 - 0.864 * x^4 / cc^4 + 0.208 * x^6 / cc^6 - 0.012 * x^8 / cc^8) / 3.25
  tmp[abs(x) < 2*cc] <- 0
  tmp[abs(x) > 3*cc] <- 2
  tmp
  
}

# --------------------------------------------------------------------



#' Variance component model prediction
#' 
#' @param y a vector of response variable
#' @param K a list of covariance kernels
#' @param W a vector of parameters (beta, sigma^2)
#' @param X matrix of covariate, intercept included
#' @param rt return a table? if false, a named vector is returned
.vpd <- function(y, K=NULL, W=NULL, X=NULL, rt=1)
{
    y <- unname(y)
    N <- NROW(y)

    ## fixed effect?
    if(!is.null(X))
    {
        xb <- X %*% W[1:NCOL(X)]
        W <- W[-(1:NCOL(X))]
    }
    else
        xb <- 0
    e <- y - xb

    ## use null model?
    if(is.null(W))
        W <- c(EPS=sum(e^2) / N)

    ## prepand noisy kernel
    C <- c(list(EPS=diag(N)), cv=K)
    
    ## make predictions
    V <- unname(Reduce(`+`, mapply(`*`, C, W, SIMPLIFY=FALSE)))
    alpha <- solve(V, e)                # V^{-1}e
    f <- V - diag(W[1], length(e))      # W[1] is ESP

    ## prediction: conditional mean and covariance
    h <- f %*% alpha + xb
    mht <- mean(h)
    
    mse <- mean((y - h)^2)
    cyh <- if(sd(h) == 0) 0 else cyh <- cor(y, h)
    rsq <- cyh^2
    
    ## negative likelihood
    ldt <- with(determinant(V), modulus * sign) / N
    eae <- sum(alpha * e) / N           # e^T V^{-1} e
    nlk <- eae + ldt
    
    ## return
    rpt <- c(rsq=rsq, mse=mse, nlk=nlk, cyh=cyh, mht=mht, ssz=N)
    if(rt == 1)
        rpt <- data.frame(key=names(rpt), val=rpt, row.names=names(rpt))
    if(rt == 2)
        rpt <- data.frame(t(rpt))
    rpt
}

#' Test GCTA
#'
#' Generate a training and testing sample of given size N and (genomic) feature size
#' P. Develop a LMM model of two kernels (linear + quadratic) and 3 fixed effects on
#' the training data, evaluate this model on the testing data.
#' 
#' @param N the sample size
#' @param P the number of genomic features (i.e., SNP)
#' @param e the size of noise
#' 
#' @return a list, containing the estimated fixed effect and variance components,
#' the performance on both training and testing data.
#' @export
gcta.test <- function(N=500, P=2000, e=5.0)
{
    ## ------------------------------ generation ------------------------------ ##
    ## design matrix for fix effect
    X.dvp <- cbind(X00=1, X01=rbinom(N, 1, .3), X02=rbinom(N, 1, .5), X03=rnorm(N))
    X.evl <- cbind(X00=1, X01=rbinom(N, 1, .3), X02=rbinom(N, 1, .5), X03=rnorm(N))
    bts <- c(X00=-0.6, X01=0.5, X02=1.0, X03=-1.0) # beta
    M.dvp <- X.dvp %*% bts                         # mean
    M.evl <- X.evl %*% bts                         # mean
    
    ## design matrix for random effect
    Z.dvp <- matrix(rpois(N * P, 2), N, P)
    Z.evl <- matrix(rpois(N * P, 2), N, P)
    K.dvp <- list(LN1=tcrossprod(scale(Z.dvp)) / P) # 1st order
    K.dvp$LN2 <- K.dvp$LN1^2                        # 2nd order

    K.evl <- list(LN1=tcrossprod(scale(Z.evl)) / P) # 1st order
    K.evl$LN2 <- K.evl$LN1^2                        # 2nd order

    eps <- c(EPS=e)                     # true noise
    vcs <- c(LN1=1.3, LN2=1.0)          # true effect

    ## matrix of covariance
    S.dvp <- vcs[1] * K.dvp$LN1 + vcs[2] * K.dvp$LN2 # Sigma
    S.evl <- vcs[1] * K.evl$LN1 + vcs[2] * K.evl$LN2 # Sigma
    
    ## response
    y.dvp <- MASS::mvrnorm(1, M.dvp, S.dvp) + rnorm(N, 0, sqrt(eps))
    y.evl <- MASS::mvrnorm(1, M.evl, S.evl) + rnorm(N, 0, sqrt(eps))

    ## ------------------------------ working fit ------------------------------ ##
    ## working design matrix and kernels
    X.use <- X.dvp # [, -1]
    K.use <- list(EPS=diag(N), LN1=K.dvp$LN1, LN2=K.dvp$LN1^2)

    md4 <- gcta.reml(y.dvp, K.use[-1], X.use[, -1]) # call GCTA
    par <- rbind(gct=md4$par, ref=c(bts, eps, vcs))

    ## ------------------------------ testing ------------------------------ ##
    X.use <- X.evl # [, -1]
    K.use <- list(EPS=diag(N), LN1=K.evl$LN1, LN2=K.evl$LN1^2)
    pd4 <- .vpd(y.evl, K=K.use[-1], md4$par, X.use)
    
    ret <- list(par=par, dvp=md4$rpt, evl=pd4)
    ret
}

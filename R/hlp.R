## helper function to make kernels
mks <- function(zmx, vcs)
{
    ncv <- length(vcs)
    cv1 <- tcrossprod(scale(zmx, TRUE, TRUE)) / NCOL(zmx)
    ret <- lapply(seq.int(ncv), function(p)
    {
        kp <- cv1^p
        kp <-  kp / mean(diag(kp))
        kp <- with(eigen(kp, TRUE), tcrossprod(vectors %*% diag(pmax(values, 0)), vectors))
        kp
    })
    names(ret) <- names(vcs)
    ret
}
## helper ofunction to combine kernels
cmb <- function(kns, vcs) Reduce(`+`, mapply(`*`, kns, vcs, SIMPLIFY=FALSE))

## helper function to simulate MVN
mvn <- function(mcv)
{
    N <- nrow(mcv)
    e <- eigen(mcv, symmetric = TRUE)
    y <- rnorm(N)
    e$vectors %*% diag(sqrt(pmax(e$values, 0)), N) %*% y
}

#' Test GCTA
#'
#' Generate a training and testing sample of given size N and (genomic) feature size
#' P. Develop a LMM model of two kernels (linear + quadratic) and 3 fixed effects on
#' the training data, evaluate this model on the testing data.
#' 
#' @param N the sample size
#' @param P the number of genomic features (i.e., SNP)
#' @param eps the size of noise
#' @param vcs true variances components to generate data
#' @param use number of components used for model development.
#' 
#' @return a list, containing the estimated fixed effect and variance components,
#' the performance on both training and testing data.
#' @export
gcta.test <- function(N=500, P=2000, eps=2.0, vcs=c(0.5, 1.0, 2.0), use=3)
{
    ## ------------------------------ generation ------------------------------ ##
    ## design matrix for fix effect
    X.dvp <- cbind(X00=1, X01=rbinom(N, 1, .3), X02=rbinom(N, 1, .5), X03=rnorm(N))
    X.evl <- cbind(X00=1, X01=rbinom(N, 1, .3), X02=rbinom(N, 1, .5), X03=rnorm(N))
    bts <- c(X00=-0.6, X01=0.5, X02=1.0, X03=-1.0) # beta
    M.dvp <- X.dvp %*% bts                         # mean
    M.evl <- X.evl %*% bts                         # mean
    names(eps) <- "EPS"                            # true noise
    names(vcs) <- sprintf("K%02d", seq_along(vcs))
    
    ## design matrix for random effect
    Z.dvp <- matrix(rpois(N * P, 2), N, P)
    K.dvp <- mks(Z.dvp, vcs)
    Z.evl <- matrix(rpois(N * P, 2), N, P)
    K.evl <- mks(Z.evl, vcs)
    
    ## matrix of covariance
    S.dvp <- cmb(K.dvp, vcs)            # Sigma
    S.evl <- cmb(K.evl, vcs)            # Sigma

    ## response
    y.dvp <- M.dvp + mvn(S.dvp) + rnorm(N, 0, sqrt(eps))
    y.evl <- M.evl + mvn(S.evl) + rnorm(N, 0, sqrt(eps))

    ## ------------------------------ working fit ------------------------------ ##
    ## working design matrix and kernels
    X.use <- X.dvp[, -1]
    K.use <- K.dvp[seq.int(use)]

    md1 <- gcta.reml(y.dvp, K.use, X.use, zbd=0) # GCTA
    md2 <- gcta.reml(y.dvp, K.use, X.use, zbd=1) # GCTA
    ref <- c(bts, eps, vcs)
    par <- rbind(md1$par, md2$par)
    
    ret <- list(ref=ref, par=par)
    ret
}

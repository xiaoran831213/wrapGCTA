## GCTA wrappers
.find.gcta <- function(download=TRUE)
{
    os <- Sys.info()['sysname']
    if(grepl("Windows", os, TRUE))
        os <- 'win'
    else if(grepl("Darwin", os, TRUE))
        os <- 'mac'
    else
        os <- 'lnx'

    ex <- file.path('gcta', os, 'gcta64')
    if(os == 'win')
        ex <- paste0(ex, '.exe')

    ## must be the one and only executables for a specific OS
    if(!file.exists(ex) && download)
    {
        ## download GCTA64 1.9 beta
        cat("Download GCTA64 for ", os, "\n", sep="")

        ## URL
        url <- "https://cnsgenomics.com/software/gcta/gcta_1.91.7beta"
        url <- paste0(url, switch(os, win="_win", mac="_mac", ""), ".zip")
        cat("URL: ", url, '\n', sep="")

        ## download
        zf <- tempfile(fileext=".zip")
        download.file(url, zf)

        ## unzip
        xf <- grep("gcta[0-9][0-9]", unzip(zf, list=TRUE)$Name, value=TRUE)
        unzip(zf, xf, junkpaths=TRUE, overwrite=TRUE, exdir=dirname(ex))
        unlink(zf)

        if(os == 'lnx' || os == 'mac')
            Sys.chmod(ex, mode = "0777", use_umask = TRUE)

        if(!file.exists(ex))
            stop("could not fetch GCTA64.")
    }
    ex
}

#' Clean GCTA Temporaries
#'
#' @export
gcta.clean <- function() unlink('gcta', TRUE, TRUE)


#' GCTA REML
#' 
#' @param y vector of response
#' @param K list of N x N kernel matrices
#' @param X design matrix of fixed effect, numerically coded
#' @param maxit maximum number of iterations
#' @param rm.temp remove temporary files (relatedness matrices in GCTA format.)
#'
#' @return list of estimates, report, predicted effect, and log.
#' @export
gcta.reml <- function(y, K, X=NULL, maxit=100, rm.temp=TRUE)
{
    ## a function to write tab-delimited table
    WT <- function(x, f) write.table(x, f, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

    ## printf
    PF <- function(...) cat(sprintf(...))

    ## temporary working directory for GCTA
    wd <- paste0("gcta", tempfile("", ''))
    if(!dir.create(wd, FALSE, TRUE))
        stop('failed to create directory', GRM.dir)
    N <- NROW(y)

    ## save GRM
    GRM <- K
    if(is.null(names(GRM)))
        names(GRM) <- sprintf('G%02d', seq(length(GRM)))
    K <- c(list(EPS=diag(N)), K)
    
    GRM.dir <- file.path(wd, 'grm')
    if(!dir.create(GRM.dir))
        stop('failed to create directory', GRM.dir)
    GRM.path <- file.path(GRM.dir, names(GRM))
    for(i in seq_along(GRM))
    {
        .saveGRM(GRM.path[i], GRM[[i]])
        PF("save GRM: %d %16s %8d %8d\n", i, GRM.path[i], nrow(GRM[i]), ncol(GRM[i]))
    }
    mgrm.path <- file.path(wd, 'grm.lst')
    write(GRM.path, mgrm.path)

    ## write phenotype
    id <- .makeID(GRM[[1]])
    phe.path <- file.path(wd, 'phe.txt')
    WT(cbind(id, y), phe.path) # phenotype

    ## covariate, assume X is numerically coded 
    if(!is.null(X))
    {
        qcv.path <- file.path(wd, 'qcv.txt')
        WT(cbind(id, X), qcv.path)
    }

    ## compose GCTA command
    exe <- .find.gcta()                 # binary
    fun <- '--reml'                     # REML
    if(length(GRM) == 1)                # GRM(s)
        mgm <- paste('--grm', GRM.path[[1]])
    else
        mgm <- paste('--mgrm', mgrm.path)
    phe <- paste('--pheno', phe.path)

    if(!is.null(X))
        qcv <- paste('--qcovar', qcv.path)
    else
        qcv <- NULL
    
    out <- paste('--out', file.path(wd, 'out'))

    opt <- paste('--reml-pred-rand', '--reml-no-lrt', '--thread-num 1', '--reml-est-fix')
    if(!is.null(maxit))
        opt <- paste(opt, '--reml-maxit', maxit)

    cmd <- paste(exe, fun, mgm, phe, qcv, out, opt)

    ## execute command
    t0 <- Sys.time()
    ext <- system(cmd)
    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    if(ext != 0)
    {
        if(rm.temp)
            unlink(wd, TRUE, TRUE)
        stop(cmd, ' GCTA existed with non-zero.')
    }

    ## prepend intercept (GCTA does this by itself)
    X <- cbind(X00=rep(1, N), X)

    ## parse the output
    out <- within(.gcta.parse(wd),
    {
        names(par) <- c(colnames(X), names(K))
        names(se1) <- names(par)
    })

    ## prediction
    y <- y - X %*% out$par[seq(ncol(X))]
    vcs <- out$par[-seq(ncol(X))]
    
    ## summary
    h <- rowSums(out$blp[-1:-3])
    v <- Reduce(`+`, mapply(`*`, K, vcs, SIMPLIFY=FALSE))
    u <- chol(v)
    a <- chol2inv(u)

    ## mse <- mean((y - h)^2)           # mean squre error
    mse <- mean(out$blp[, 3]^2)         # estimate residual
    cyh <- cor(y, h)                    # correlation
    rsq <- cyh^2
    ## nlk <- .5 * sum(crossprod(y, a) * y) + sum(log(diag(u))) + .5 * N * log(2 * pi)
    nlk <- crossprod(y, a) %*% y + 2 * sum(log(diag(u)))
    nlk <-  nlk / N                     # NLK

    ## h <- y - a %*% y / diag(a)
    ## loo <- mean((y - h)^2)

    rpt <- data.frame(
        key=c('mse', 'nlk', 'cyh', 'rsq', 'rtm', 'ssz'),
        val=c(mse, nlk, cyh, rsq, td, N))
    rownames(rpt) <- rpt$key

    ## remove temporary files
    if(rm.temp)
        unlink(wd, TRUE, TRUE)

    ## pack up and return
    c(list(cmd=cmd, ext=ext, rpt=rpt), out)
}

#' Parse GCTA output
#'
#' @param wd the working directory for GCTA
.gcta.parse <- function(wd)
{
    ## --- parse *.hsq for variance ---
    hsq.path <- file.path(wd, 'out.hsq')
    hsq <- readLines(hsq.path)

    ## variance components, and their names
    reg <- "(^Source|^V[(])"
    vcs <- grep(reg, hsq, value=TRUE)
    hsq <- grep(reg, hsq, value=TRUE, invert=TRUE)
    vcs <- read.delim(text=vcs)
    names(vcs) <- c('par', 'val', 'se1')
    vcs <- within(vcs, par <- sub("V[(](.*)[)]", "\\1", par))
    vcs <- subset(vcs, !grepl('Vp$', par))
    vcs <- vcs[c(nrow(vcs), 1:(nrow(vcs) - 1)), ]
    vcs$par[1] <- 'EPS'

    ## fixed effects
    reg <- "(^Fix_eff|^[-|0-9.])"
    fix <- grep(reg, hsq, value=TRUE)
    hsq <- grep(reg, hsq, value=TRUE, invert=FALSE)
    fix <- read.delim(text=fix)
    fix <- cbind(par=sprintf('f%02d', seq(0, nrow(fix)-1)), fix)
    names(fix) <- c('par', 'val', 'se1')
    
    ## --- parse *.indi.blp for predicted genetic effect ---
    blp.path <- file.path(wd, 'out.indi.blp')
    blp <- read.table(blp.path)
    blp <- blp[, c(1, 2, seq(4, ncol(blp), 2))] # effects
    blp <- blp[, c(1:2, ncol(blp), seq(3, ncol(blp) - 1))]
    names(blp)[-(1:2)] <- vcs$par

    ## --- parse *.log for time elapsed
    log.path <- file.path(wd, 'out.log')
    log <- readLines(log.path)
    
    ## return
    est <- rbind(fix, vcs)
    par <- est$val
    se1 <- est$se1
    list(par=par, se1=se1, hsq=hsq, blp=blp, log=log)
}

## GCTA wrappers
.find.gcta <- function(download=TRUE)
{
    ## find OS type
    os <- Sys.info()['sysname']
    if(grepl("Windows", os, TRUE))
        os <- 'win'
    else if(grepl("Darwin", os, TRUE))
        os <- 'mac'
    else
        os <- 'lnx'

    ## try gcta on OS search path
    ex <- Sys.which('gcta64')
    if(file.exists(ex))
    {
        ## cat("Find gcta64 at ", ex, "\n", sep="")
        return('gcta64')
    }
    if(file.exists(Sys.which('gcta')))
    {
        ## cat("Find gcta at ", ex, "\n", sep="")
        return('gcta')
    }
    ## try project path
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
#' @param zbd TRUE to lower bound the variance component to 0 (def=TRUE).
#' @param itr maximum number of iterations
#' @param alg integer to specify the alorithm to be used, 0=Average information,
#' 1=Fisher Scoring, 2=EM, corresponding to option --reml-alg.
#' @param quiet TRUE to execute GCTA quietly.
#' @param rm.temp remove temporary files (relatedness matrices in GCTA format.)
#' after computation.
#'
#' @return list of estimates, report, predicted effect, and log.
#' @export
gcta.reml <- function(y, K, X=NULL, zbd=TRUE, itr=100, alg=1, quiet=TRUE, rm.temp=TRUE)
{
    WT <- function(x, f) write.table(x, f, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
    PF <- function(...) cat(sprintf(...))

    ## temporary working directory for GCTA
    wd <- paste0("gcta", tempfile("", ''))
    if(!dir.create(wd, FALSE, TRUE))
        stop('failed to create directory', wd)
    N <- NROW(y)

    ## write GRM files
    if(is.null(names(K)))
        names(K) <- sprintf('G%02d', seq_along(K))
    
    dir.grm <- file.path(wd, 'grm')
    if(!dir.create(dir.grm)) stop('failed to create directory', dir.grm)
    fns.grm <- file.path(dir.grm, names(K))
    for(i in seq_along(K))
    {
        .saveGRM(fns.grm[i], K[[i]])
        PF("Write GRM: %d %16s %8d %8d\n", i, fns.grm[i], N, N)
    }

    ## write phenotype
    id <- .makeID(K[[1]])
    fns.phe <- file.path(wd, 'phe.txt')
    WT(cbind(id, y), fns.phe) # phenotype

    ## covariate, assume X is numerically coded
    ## number of fix effect
    nfx <- if(length(X) > 0) NCOL(X) else 0
    if(nfx > 0L)
    {
        fns.qcv <- file.path(wd, 'qcv.txt')
        WT(cbind(id, X), fns.qcv)
    }
    
    ## RUN GCTA, recursively
    tm <- 0
    vcs <- rep(1, length(K))            # init VCs as positive
    eps <- 1
    fns.mgr <- file.path(wd, 'mgr')     # list of GRMs

    ## compose GCTA command
    cmd <- paste(.find.gcta(), "--reml")

    ## the list of kernels, drop those corresponding to non-positive
    ## variance components estimated by previouse GCTA RUN
    write(fns.grm, fns.mgr)
    cmd <- paste(cmd, '--mgrm',  fns.mgr)
    cmd <- paste(cmd, '--pheno', fns.phe)
    if(length(X) > 0L)
        cmd <- paste(cmd, '--qcovar', fns.qcv)
    cmd <- paste(cmd, '--out', file.path(wd, 'out'))

    ## iteration limit
    if(!is.null(itr))
        cmd <- paste(cmd, '--reml-maxit', itr)
    ## algorithms (0=AI, 1=FS, 2=EM)
    if(!is.null(alg))
        cmd <- paste(cmd, '--reml-alg', alg)
    if(!is.null(zbd) && !zbd)
        cmd <- paste(cmd, '--reml-no-constrain')

    ## other options
    opt <- paste('--reml-pred-rand', '--reml-no-lrt', '--thread-num 1',
                 '--reml-est-fix')
    cmd <- paste(cmd, opt)

    ## execute command
    t0 <- Sys.time()
    ext <- system(cmd, intern=FALSE, ignore.stdout=quiet)
    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    tm <- tm + td

    if(ext != 0)                    # error ?
    {
        err <- grep('^error:', readLines(file.path(wd, 'out.log')), TRUE, value=TRUE)
        if(rm.temp)
            unlink(wd, TRUE, TRUE)
        stop(cmd, ' GCTA existed with non-zero:', err)
    }

    ## success
    out <- gcta.parse(wd)
    par <- out$par
    names(par)[seq(1 + nfx)] <- if(nfx > 0) c("X00", colnames(X)) else "X00"
    names(par)[2 + nfx] <- "EPS"
    names(par)[seq(3 + nfx, length(par))] <- names(K)
    cat(paste(sprintf("%-9s", round(par, 5), collapse=" ")), "\n", sep="")

    ## remove temporary files
    if(rm.temp)
        unlink(wd, TRUE, TRUE)

    ## return
    list(cmd=cmd, ext=ext, par=par, rtm=tm)
}

#' Parse GCTA output
#'
#' @param wd the working directory for GCTA
gcta.parse <- function(wd)
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

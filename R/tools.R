.onLoad <- function(libname,pkgname) {
  if(is.na(thr <- as.integer(Sys.getenv('LATENTCOR_THREADS')))) thr <- 1L
  options(latentcor.threads=thr)
}

vectorfun <- function(fun,arity,domain=NULL) {
  force(arity)
  f <- function(x) {
    mc <- match.call(expand.dots=TRUE)
    mc[[1L]] <- quote(list)
    arglist <- eval.parent(mc)
    x <- arglist[[1]]
    if(is.matrix(x) || length(x) == arity) return(do.call(fun,arglist))
    if(arity == 1) {arglist[[1]] <- matrix(x,1); return(do.call(fun,arglist))}
    # try to coerce argument to matrix
    if(length(x) %% arity == 0) {
      numvec <- length(x) / arity
      if(numvec > 1)
        warning(sprintf('Coercing vector of length %d into %dx%d matrix',length(x),arity,numvec))
      dim(x) <- c(arity,numvec)
      arglist[[1]] <- x
      do.call(fun,arglist)
    } else
      stop(sprintf('Function should take %d arguments, you supplied a vector of length %d',arity,length(x)))
  }
  formals(f) <- formals(fun)
  structure(f,arity=arity,domain=as.data.frame(domain),
            latentcor.version=utils::packageVersion('latentcor'))
}

defaultblend <- c('linear','cubic','sigmoid','parodic','square','mean')
blenddef <- function(fun,blend=defaultblend) {
  if(is.null(blend)) return(fun)
  blend <- match.arg(blend)
  pos <- match(blend,defaultblend)
  f <- formals(fun)
  if(!('blend' %in% names(f))) return(f)
  f[['blend']] <- as.call(c(list(as.name('c')),c(defaultblend[pos],defaultblend[-pos])))
  formals(fun) <- f
  ffun <- get('fun',environment(fun))
  ff <- formals(ffun)
  if(!('blend' %in% names(ff))) return(fun)
  ff[['blend']] <- f[['blend']]
  formals(ffun) <- ff
  assign('fun',ffun,envir=environment(fun))
  fun
}





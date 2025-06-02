# this script contains the base functions that can be used to efficiently calculate the MDS coordinates based on a pairwise distance matrix.

#' get mds from distance matrix
get_mds_wdmat <- function(dmat, mds_method = c("cmds", "metmds", "isomds", "monomds"), ndim = 2, smacof_itmax = 200, smacof_eps = 1e-6, isomds_maxit = 200, isomds_tol = 1e-6) {
  
  mds_method <- match.arg(mds_method, c("cmds", "metmds", "isomds", "monomds"), several.ok = TRUE)
  
  mds <- NULL
  cmds <- cmdscale2(dmat, k = max(ndim))
  if ("cmds" %in% mds_method) {
    mds <- list(cmds = cmds)
  }
  
  if ("metmds" %in% mds_method) {
    metmds_alldim <- vector("list", length(ndim))
    for (i in 1:length(ndim)) {
      cmds_this <- cmds[, 1:ndim[i]]
      metmds_alldim[[i]] <- smacofSym2(dmat, ndim = ndim[i], init = cmds_this, type = "ratio", itmax = smacof_itmax, eps = smacof_eps, verbose = TRUE)
    }
    if (length(ndim) == 1) metmds_alldim <- metmds_alldim[[1]]
    mds <- c(mds, list(metmds = metmds_alldim))
  }
  
  if ("isomds" %in% mds_method) {
    isomds_alldim <- vector("list", length(ndim))
    for (i in 1:length(ndim)) {
      cmds_this <- cmds[, 1:ndim[i]]
      isomds_alldim[[i]] <- isoMDS2(dmat, y = cmds_this, k = ndim[i], maxit = isomds_maxit, tol = isomds_tol)
    }
    if (length(ndim) == 1) isomds_alldim <- isomds_alldim[[1]]
    mds <- c(mds, list(isomds = isomds_alldim))
  }
  
  if ("monomds" %in% mds_method) {
    monomds_alldim <- vector("list", length(ndim))
    for (i in 1:length(ndim)) {
      cmds_this <- cmds[, 1:ndim[i]]
      monomds_alldim[[i]] <- vegan::monoMDS(dmat, y = cmds_this, k = ndim[i], maxit = isomds_maxit)
      incidx <- which(!names(monomds_alldim[[i]]) %in% c("diss", "iidx", "jidx", "dist", "dhat"))
      monomds_alldim[[i]] <- monomds_alldim[[i]][incidx]
    }
    if (length(ndim) == 1) monomds_alldim <- monomds_alldim[[1]]
    mds <- c(mds, list(monomds = monomds_alldim))
  }
  
  return(mds)
}

#' faster implementation of the cmdscale function
cmdscale2 <- function (d, k = 2)
{
  if (is.null(n <- attr(d, "Size"))) {
    x <- as.matrix(d^2)
    storage.mode(x) <- "double"
    if ((n <- nrow(x)) != ncol(x)) stop("distances must be result of 'dist' or a square matrix")
    rn <- rownames(x)
  } else {
    rn <- attr(d, "Labels")
    x <- matrix(0, n, n)
    x[row(x) > col(x)] <- d^2
    x <- x + t(x)
  }
  rm(d)
  gc()
  gc()
  
  n <- as.integer(n)
  # if (is.na(n) || n > 46340) 
  if (is.na(n) || n > 1e5) stop("invalid value of 'n'")
  if ((k <- as.integer(k)) > n - 1 || k < 1) stop("'k' must be in {1, 2, ..  n - 1}")
  x <- .Call(stats:::C_DoubleCentre, x)
  
  e <- RSpectra::eigs_sym(-x/2, k, which="LA")
  rm(x)
  gc()
  gc()
  
  ev <- e$values
  evec <- e$vectors
  
  k1 <- sum(ev > 0)
  if (k1 < k) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", k1, k), domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
  }
  points <- evec * rep(sqrt(ev), each = n)
  dimnames(points) <- list(rn, NULL)
  points
}


isoMDS2 <- function (d, y = cmdscale(d, k), k = 2, maxit = 50, trace = TRUE, tol = 0.001, p = 2) 
{
  if (any(!is.finite(d)) && missing(y)) stop("an initial configuration must be supplied with NA/Infs in 'd'")
  if (!is.matrix(y)) stop("'y' must be a matrix")
  
  if (is.null(n <- attr(d, "Size"))) {
    x <- as.matrix(d)
    if ((n <- nrow(x)) != ncol(x)) stop("distances must be result of 'dist' or a square matrix")
    rn <- rownames(x)
  } else {
    x <- matrix(0, n, n)
    x[row(x) > col(x)] <- d
    x <- x + t(x)
    rn <- attr(d, "Labels")
  }
  rm(d)
  gc()
  gc()
  
  n <- as.integer(n)
  if (is.na(n)) stop("invalid size")
  ab <- x[row(x) < col(x)] <= 0
  if (any(ab, na.rm = TRUE)) {
    ab <- !is.na(ab) & ab
    aa <- cbind(as.vector(row(x)), as.vector(col(x)))[row(x) < col(x), ]
    aa <- aa[ab, , drop = FALSE]
    stop(gettextf("zero or negative distance between objects %d and %d", aa[1, 1], aa[1, 2]), domain = NA)
  }
  
  if (any(dim(y) != c(n, k))) stop("invalid initial configuration")
  if (any(!is.finite(y))) stop("initial configuration must be complete")
  
  ord <- order(x[row(x) > col(x)])
  nd <- sum(!is.na(ord))
  rm(x)
  gc()
  gc()
  
  on.exit(.C(MASS:::VR_mds_unload))
  .C(MASS:::VR_mds_init_data, as.integer(nd), as.integer(k), n, as.integer(ord - 1), as.integer(order(ord) - 1), as.double(y), as.double(p))
  tmp <- .C(MASS:::VR_mds_dovm, val = double(1), as.integer(maxit), as.integer(trace), y = as.double(y), as.double(tol))
  points <- matrix(tmp$y, , k)
  dimnames(points) <- list(rn, NULL)
  list(points = points, stress = tmp$val)
}


#' faster and more memory efficient implementation of the smacof function
smacofSym2 <- function (delta, ndim = 2, type = "ratio", init = "torgerson", verbose = FALSE, itmax = 1000, eps = 1e-06) 
{
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- as.dist(diss)
    attr(diss, "Labels") <- rownames(delta)
  }
  rm(delta)
  gc()
  gc()
  
  p <- ndim
  n <- attr(diss, "Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
  nn <- n * (n - 1)/2
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)
  
  x <- initConf2(init, diss, n, p)
  dhat <- sqrt(length(diss)/sum(diss^2)) * diss

  itel <- 1
  d <- dist(x)
  lb <- sum(d * dhat)/sum(d^2)
  x <- lb * x
  d <- lb * d
  sold <- sum((dhat - d)^2)/nn
  
  repeat {
    b <- bmat2(dhat, d, n)
    y <-  1/n * (b %*% x)
    e <- dist(y)
    snon <- sum((dhat - e)^2)/nn
    if (verbose) 
      cat("Iteration: ", formatC(itel, width = 3, format = "d"), " Stress (raw):", formatC(c(snon), digits = 8, width = 12, format = "f"), " Difference: ", 
          formatC(sold - snon, digits = 8, width = 12, format = "f"), "\n")
    if (((sold - snon) < eps) || (itel == itmax)) break
    x <- y
    d <- e
    sold <- snon
    itel <- itel + 1
  }
  
  stress <- sqrt(snon)
  colnames(y) <- paste("D", 1:(dim(y)[2]), sep = "")
  rownames(y) <- labels(diss)

  if (itel == itmax) warning("Iteration limit reached! You may want to increase the itmax argument!")
  result <- list(conf = y, stress = stress, ndim = p, model = "Symmetric SMACOF", niter = itel, nobj = n, type = type, call = match.call())
  class(result) <- c("smacofB", "smacof")
  result
}

torgerson2 <- function (delta, p = 2) 
{
  z <- RSpectra::eigs_sym(-.Call(stats:::C_DoubleCentre, delta^2)/2, p, which="LA")
  v <- pmax(z$values, 0)
  if (p == 1) normdiag <- cbind(sqrt(v[1]))
  else normdiag <- diag(sqrt(v[1:p]))
  conf <- z$vectors[, 1:p] %*% normdiag
  rownames(conf) <- rownames(diss)
  return(conf)
}

initConf2 <- function (init, diss, n, p, inddiff = FALSE) 
{
  if (inddiff) diss <- as.dist(apply(simplify2array(lapply(diss, as.matrix)), c(1, 2), sum, na.rm = TRUE))
  
  if (length(init) == 1) {
    if (init == "torgerson") {
      meandiss <- mean(diss, na.rm = TRUE)
      diss1 <- as.matrix(diss)
      diss1[is.na(diss1)] <- meandiss
      x <- torgerson2(diss1, p = p)
      init <- "dummy"
    }
    if (init == "random") {
      x <- matrix(runif(n * p, min = -1), ncol = p)
    }
  }
  if (is.data.frame(init)) 
    init <- as.matrix(init)
  if (is.matrix(init)) {
    x <- as.matrix(init)
    if (any(dim(x) != c(n, p))) stop(paste0("Dimension of the starting configuration matrix needs to be ", n, " times ", p, "!"))
  }
  return(x)
}

bmat2 <- function (diss, d, n, eps = 1e-12)
{
  z <- d < eps
  df <- as.matrix(diss * (!z)/(d + z))
  r <- rowSums(df)
  return(.Internal(diag(r, n, n)) - df)
}

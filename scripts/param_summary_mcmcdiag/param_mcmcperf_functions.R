# this script contains functions that summarize parameters and compute MCMC diagnostics
library(coda)
library(rstan)
library(convenience)

# check and get the unique value from a vector with all identical elements
get_unique <- function(x) {
  x <- unique(x)
  if (length(x) != 1) {
    stop("different chain lengths")
  }
  return(x)
}

# generate result file output path
get_outpath <- function(path_str, ngen_freq, subsample_freq = 0, add_str = "") {
  if (subsample_freq > 1) {
    ngen_freq <- subsample_freq
  }
  out_path <- paste0(path_str, as.integer(ngen_freq), "sampleFreq", add_str)
  
  outdir_path <- dirname(out_path)
  out_name <- basename(out_path)
  while (nchar(out_name) > 255 && grepl("and[0-9]+(?=pruned)", out_name, perl = TRUE)) {
    out_name <- sub("and[0-9]+(?=pruned)", "", out_name, perl = TRUE)
  }
  
  out_path <- paste0(outdir_path, "/", out_name)
  return(out_path)
}

# calculate PSRF
calc_psrfval <- function(param_list, method = c("coda", "rstan")) {
  method <- match.arg(method, c("coda", "rstan"))
  if (length(param_list) <= 1) return(rep(NA_real_, switch(method, coda = 2, rstan = 1)))
  nsamples_min <- min(lengths(param_list))
  param_list <- lapply(param_list, function(x) x[1:nsamples_min])
  if (method == "coda") {
    param_list <- coda::mcmc.list(lapply(param_list, as.mcmc))
    res <- coda::gelman.diag(param_list, autoburnin = FALSE, multivariate = FALSE)$psrf[1, ]
  } else {
    res <- rstan::Rhat(do.call(cbind, param_list))
  }
  return(res)
}

# prepare log files to compute PSRF
comp_psrf <- function(rep_pbs, method = c("coda", "rstan"), pmin = 0.01) {
  method <- match.arg(method, c("coda", "rstan"))
  
  nreps <- length(rep_pbs)
  nparams <- ncol(rep_pbs[[1]])
  nrow_logs <- sapply(rep_pbs, nrow)
  
  if (length(unique(nrow_logs)) > 1) {
    nrow_max <- max(nrow_logs)
    nrow_min <- min(nrow_logs)
    rep_idx_max <- which(nrow_logs == nrow_max)
    rep_idx_min <- which(nrow_logs == nrow_min)
    if (method == "coda") {
      psrf_rowname <- c(paste0("coda_", c("point", "upper"), "_long"), paste0("coda_", c("point", "upper"), "_short"), paste0("coda_", c("point", "upper"), "_trunc"))
    } else {
      psrf_rowname <- paste0("rstan_", c("long", "short", "trunc"))
    }
    
    nsamp_min1 <- ceiling(nrow_max * pmin)
    nsamp_min2 <- ceiling(nrow_min * pmin)
  } else {
    if (method == "coda") {
      psrf_rowname <- paste0("coda_", c("point", "upper"))
    } else {
      psrf_rowname <- "rstan"
    }
    
    nsamp_min <- ceiling(unique(nrow_logs) * pmin)
  }
  
  psrf_nrow <- length(psrf_rowname)
  psrf_mat <- matrix(nrow = psrf_nrow, ncol = nparams)
  rownames(psrf_mat) <- psrf_rowname
  colnames(psrf_mat) <- colnames(rep_pbs[[1]])
  
  for (j in 1:nparams) {
    param_all <- sapply(rep_pbs, "[[", j)
    
    if (length(unique(nrow_logs)) > 1) {
      param_all_finite <- lapply(param_all, function(x) x[is.finite(x)])
      param_all_trunc <- lapply(param_all, function(x) x[1:nrow_min])
      param_all_trunc_finite <- lapply(param_all_trunc, function(x) x[is.finite(x)])
      
      param_all_longruns <- param_all_finite[rep_idx_max]
      param_all_longruns <- param_all_longruns[which(lengths(param_all_longruns) >= nsamp_min1)]
      param_all_shortruns <- param_all_finite[rep_idx_min]
      param_all_shortruns <- param_all_shortruns[which(lengths(param_all_longruns) >= nsamp_min2)]
      param_all_allrunstrunc <- param_all_trunc_finite[which(lengths(param_all_trunc_finite) >= nsamp_min2)]
      psrf_mat[, j] <- unlist(lapply(list(param_all_longruns, param_all_shortruns, param_all_allrunstrunc), calc_psrfval, method = method))
    } else {
      param_all <- lapply(rep_pbs, function(x) x[[j]][is.finite(x[[j]])])
      param_all <- param_all[which(lengths(param_all) >= nsamp_min)]
      psrf_mat[, j] <- calc_psrfval(param_all, method = method)
    }
  }
  
  return(psrf_mat)
}

# call PSRF computation functions to calculate various types of PSRF
get_psrf <- function(rep_pbs, ngen_freq, subsample_freq = 0) {
  
  if (subsample_freq != ngen_freq && subsample_freq > 1) {
    if (subsample_freq < ngen_freq) {
      cat("subsampling frequency is not greater than the current sampling frequncy, so do nothing\n")
      return()
    } else if (subsample_freq %% ngen_freq != 0) {
      warning("subsampling frequency not multiple of sampled generation frequency")
      return()
    }
    
    by <- as.integer(subsample_freq / ngen_freq)
    rep_pbs <- lapply(rep_pbs, function(x) x[seq(from = by, to = nrow(x), by = by), , drop = FALSE])
  }
  
  psrf_mat <- do.call(rbind, lapply(c("coda", "rstan"), function(x) comp_psrf(rep_pbs, method = x)))
  return(psrf_mat)
}

# read in log files and call downstream PSRF functions
generate_psrf_all <- function(rep_paths, burnin = 2e8, ngen_trunc = 0, overwrite = FALSE) {
  
  if (burnin <= 0) {
    burnin_str <- ""
  } else {
    if (burnin < 1) {
      burnin_str <- paste0("_", burnin * 100, "burninPerc")
    } else {
      burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
    }
  }
  
  outdir_path <- gsub("/analyses", "/intermediate", paste0(getwd(), "/", gsub("^\\.", "", dirname(rep_paths[1]))))
  if (!dir.exists(outdir_path)) dir.create(outdir_path, recursive = TRUE)
  name_str <- gsub("_run(\\d+).*", "", basename(rep_paths[1]))
  path_str <- paste0(outdir_path, "/", "psrf_", name_str, burnin_str, "_")
  
  subsample_freq <- as.integer(c(1000, 10000, 100000))
  rep_2rows <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F, nrows = 2))
  ngen_freq <- get_unique(sapply(rep_2rows, function(x) diff(x[1:2, 1])))
  subsample_freq <- subsample_freq[subsample_freq >= ngen_freq]
  if (length(subsample_freq) == 0) return()
  
  ngen_truncstr <- ifelse(ngen_trunc < ngen_freq, "", paste0(as.integer(ngen_trunc), "keepNgen", "_"))
  path_str <- paste0(path_str, ngen_truncstr)
  
  if (overwrite == FALSE) {
    out_paths <- sapply(subsample_freq, function(x) get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = x, add_str = ".tsv"))
    subsample_freq <- subsample_freq[which(file.exists(out_paths) == FALSE)]
  }
  nsets <- length(subsample_freq)
  if (nsets == 0) return()
  
  nreps <- length(rep_paths)
  rep_dfs <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F))
  
  nrow_logs <- sapply(rep_dfs, nrow)
  ngen_totals <- sapply(1:nreps, function(i) rep_dfs[[i]][nrow_logs[i], 1])
  
  if (burnin <= 0) {
    rep_pbs <- rep_dfs
  } else {
    if (burnin < 1) {
      ngen_burnins <- floor(ngen_totals * burnin) + 1L
    } else {
      ngen_burnins <- rep(as.integer(burnin), nreps)
    }
    
    if (ngen_trunc < ngen_freq) {
      ngen_truncs <- ngen_totals
    } else {
      ngen_truncs <- ngen_burnins + ngen_trunc
      if (any(ngen_truncs > ngen_totals)) stop("The specified number of generations to truncate is longer than the total number of generations.")
    }
    
    rep_pbs <- lapply(1:nreps, function(i) rep_dfs[[i]][rep_dfs[[i]][, 1] > ngen_burnins[i] & rep_dfs[[i]][, 1] <= ngen_truncs[i], -1, drop = FALSE])
  }
  
  rm(rep_dfs)
  gc()
  
  for (j in 1:nsets) {
    out_path <- get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = subsample_freq[j], add_str = ".tsv")
    if (file.exists(out_path) == FALSE || overwrite) {
      cat("working on ", out_path, "\n", sep = "")
      psrf_mat <- get_psrf(rep_pbs, ngen_freq = ngen_freq, subsample_freq = subsample_freq[j])
      if (!is.null(psrf_mat)) {
        write.table(psrf_mat, file = out_path, quote = FALSE, sep = "\t")
      }
    }
  }
}


# ess
comp_ess <- function(rep_pbs, method = c("coda", "rstan", "tracer"), max_lag = 2000L, pmin = 0.01) {
  method <- match.arg(method, c("coda", "rstan", "tracer"))
  
  nreps <- length(rep_pbs)
  nparams <- ncol(rep_pbs[[1]])
  ess_mat <- matrix(nrow = nreps, ncol = nparams)
  rownames(ess_mat) <- paste0("rep", 1:nreps)
  colnames(ess_mat) <- colnames(rep_pbs[[1]])
  
  for (j in 1:nparams) {
    param_all <- lapply(rep_pbs, "[[", j)
    for (i in 1:nreps) {
      param <- param_all[[i]]
      nsamp_min <- ceiling(length(param) * pmin)
      param <- param[is.finite(param)]
      if (length(param) >= nsamp_min) {
        ess_mat[i, j] <- switch(method, coda = coda::effectiveSize(param), rstan = rstan::ess_bulk(param), tracer = convenience:::essTracerC(param, max_lag = max_lag))
      }
    }
  }
  
  nrow_pbs <- sapply(rep_pbs, nrow)
  if (length(unique(nrow_pbs)) > 1) {
    nrow_max <- max(nrow_pbs)
    nrow_min <- min(nrow_pbs)
    rep_idx_max <- which(nrow_pbs == nrow_max)
    nreps_max <- length(rep_idx_max)
    ess_mat_trunc <- matrix(nrow = nreps_max, ncol = nparams)
    rownames(ess_mat_trunc) <- paste0("rep", rep_idx_max, "trunc")
    colnames(ess_mat_trunc) <- colnames(rep_pbs[[1]])
    
    for (j in 1:nparams) {
      param_all <- lapply(rep_pbs[rep_idx_max], "[[", j)
      param_all_trunc <- lapply(param_all, function(x) x[1:nrow_min])
      for (i in 1:nreps_max) {
        param <- param_all_trunc[[i]]
        nsamp_min <- ceiling(length(param) * pmin)
        param <- param[is.finite(param)]
        if (length(param) >= nsamp_min) {
          ess_mat_trunc[i, j] <- switch(method, coda = coda::effectiveSize(param), rstan = rstan::ess_bulk(param), tracer = convenience:::essTracerC(param, max_lag = max_lag))
        }
      }
    }
    
    ess_mat <- rbind(ess_mat, ess_mat_trunc)
  }
  
  return(ess_mat)
}


get_ess <- function(rep_pbs, ngen_freq, method = c("coda", "rstan", "tracer"), max_lag = 2000L, subsample_freq = 0) {
  method <- match.arg(method, c("coda", "rstan", "tracer"))
  
  nrow_pbs <- sapply(rep_pbs, nrow)
  if (subsample_freq != ngen_freq && subsample_freq > 1) {
    if (subsample_freq < ngen_freq) {
      cat("subsampling frequency is not greater than the current sampling frequncy, so do nothing\n")
      return()
    } else if (subsample_freq %% ngen_freq != 0) {
      warning("subsampling frequency not multiple of sampled generation frequency")
      return()
    }
    
    by <- as.integer(subsample_freq / ngen_freq)
    rep_pbs <- lapply(rep_pbs, function(x) x[seq(from = by, to = nrow(x), by = by), , drop = FALSE])
    nrow_pbs <- sapply(rep_pbs, nrow)
  }
  
  if (method == "tracer" && max_lag > 2000 && min(nrow_pbs) < max_lag) {
    cat("max lag is greater than Tracer default and greater than the number of samples of the shortest run\n")
    return()
  }
  
  if (method %in% c("coda", "rstan")) {
    ess_mat <- comp_ess(rep_pbs, method = method)
  } else {
    ess_mat <- comp_ess(rep_pbs, method = method, max_lag = max_lag)
  }
  
  return(ess_mat)
}


generate_ess_all <- function(rep_paths, burnin = 2e8, ngen_trunc = 0, overwrite = FALSE) {
  
  if (burnin == 0) {
    burnin_str <- ""
  } else {
    if (burnin < 1) {
      burnin_str <- paste0("_", burnin * 100, "burninPerc")
    } else {
      burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
    }
  }
  
  outdir_path <- gsub("/analyses", "/intermediate", paste0(getwd(), "/", gsub("^\\.", "", dirname(rep_paths[1]))))
  if (!dir.exists(outdir_path)) dir.create(outdir_path, recursive = TRUE)
  name_str <- gsub("_run(\\d+).*", "", basename(rep_paths[1]))
  path_str <- paste0(outdir_path, "/", "ess_", name_str, burnin_str, "_")
  
  subsample_freq <- as.integer(c(1000, 10000, 100000))
  rep_2rows <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F, nrows = 2))
  ngen_freq <- get_unique(sapply(rep_2rows, function(x) diff(x[1:2, 1])))
  subsample_freq <- subsample_freq[subsample_freq >= ngen_freq]
  if (length(subsample_freq) == 0) return()
  
  ngen_truncstr <- ifelse(ngen_trunc < ngen_freq, "", paste0(as.integer(ngen_trunc), "keepNgen", "_"))
  path_str <- paste0(path_str, ngen_truncstr)
  
  max_lag <- as.integer(c(2000, 10000))
  set_df <- expand.grid(method = c("coda", "rstan"), subsample_freq = subsample_freq, max_lag = 0L, stringsAsFactors = FALSE)
  set_df <- rbind(set_df, expand.grid(method = c("tracer"), subsample_freq = subsample_freq, max_lag = max_lag[1], stringsAsFactors = FALSE))
  set_df <- set_df[order(set_df$subsample_freq, decreasing = TRUE), ]
  set_df <- rbind(set_df, expand.grid(method = c("tracer"), subsample_freq = subsample_freq, max_lag = max_lag[2], stringsAsFactors = FALSE))
  nsets <- nrow(set_df)
  
  if (overwrite == FALSE) {
    out_paths <- character(nsets)
    for (j in 1:nsets) {
      out_paths[j] <- get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = set_df$subsample_freq[j],
                                  add_str = paste0("_", set_df$method[j], ifelse(set_df$method[j] == "tracer", paste0("_", set_df$max_lag[j], "maxLag"), ""), ".tsv"))
    }
    set_df <- set_df[which(file.exists(out_paths) == FALSE), , drop = FALSE]
    nsets <- nrow(set_df)
    if (nsets == 0) return()
  }
  
  nreps <- length(rep_paths)
  rep_dfs <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F))
  
  nrow_logs <- sapply(rep_dfs, nrow)
  ngen_totals <- sapply(1:nreps, function(i) rep_dfs[[i]][nrow_logs[i], 1])
  
  if (burnin <= 0) {
    rep_pbs <- rep_dfs
  } else {
    if (burnin < 1) {
      ngen_burnins <- floor(ngen_totals * burnin) + 1L
    } else {
      ngen_burnins <- rep(as.integer(burnin), nreps)
    }
    
    if (ngen_trunc < ngen_freq) {
      ngen_truncs <- ngen_totals
    } else {
      ngen_truncs <- ngen_burnins + ngen_trunc
      if (any(ngen_truncs > ngen_totals)) stop("The specified number of generations to truncate is longer than the total number of generations.")
    }
    
    rep_pbs <- lapply(1:nreps, function(i) rep_dfs[[i]][rep_dfs[[i]][, 1] > ngen_burnins[i] & rep_dfs[[i]][, 1] <= ngen_truncs[i], -1, drop = FALSE])
  }
  
  rm(rep_dfs)
  gc()
  
  for (j in 1:nsets) {
    out_path <- get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = set_df$subsample_freq[j],
                            add_str = paste0("_", set_df$method[j], ifelse(set_df$method[j] == "tracer", paste0("_", set_df$max_lag[j], "maxLag"), ""), ".tsv"))
    if (file.exists(out_path) == FALSE || overwrite) {
      cat("working on ", out_path, "\n", sep = "")
      ess_mat <- get_ess(rep_pbs, ngen_freq = ngen_freq, method = set_df$method[j], max_lag = set_df$max_lag[j], subsample_freq = set_df$subsample_freq[j])
      if (!is.null(ess_mat)) {
        write.table(ess_mat, file = out_path, quote = FALSE, sep = "\t")
      }
    }
  }
}


get_mean <- function(rep_pbs, ngen_freq, subsample_freq = 0) {
  if (subsample_freq != ngen_freq && subsample_freq > 1) {
    if (subsample_freq < ngen_freq) {
      cat("subsampling frequency is not greater than the current sampling frequncy, so do nothing\n")
      return()
    } else if (subsample_freq %% ngen_freq != 0) {
      warning("subsampling frequency not multiple of sampled generation frequency")
      return()
    }
    
    by <- as.integer(subsample_freq / ngen_freq)
    rep_pbs <- lapply(rep_pbs, function(x) x[seq(from = by, to = nrow(x), by = by), , drop = FALSE])
  }
  
  nreps <- length(rep_pbs)
  mean_mat <- do.call(rbind, lapply(rep_pbs, colMeans, na.rm = TRUE))
  rownames(mean_mat) <- paste0("rep", 1:nreps)
  
  nrow_pbs <- sapply(rep_pbs, nrow)
  if (length(unique(nrow_pbs)) > 1) {
    nrow_max <- max(nrow_pbs)
    nrow_min <- min(nrow_pbs)
    rep_idx_max <- which(nrow_pbs == nrow_max)
    nreps_max <- length(rep_idx_max)
    mean_mat_trunc <- do.call(rbind, lapply(rep_pbs[rep_idx_max], function(x) colMeans(x[1:nrow_min, , drop = FALSE], na.rm = TRUE)))
    rownames(mean_mat_trunc) <- paste0("rep", rep_idx_max, "trunc")
    mean_mat <- rbind(mean_mat, mean_mat_trunc)
  }
  
  return(mean_mat)
}

generate_mean_all <- function(rep_paths, burnin = 2e8, ngen_trunc = 0, overwrite = FALSE) {
  
  if (burnin == 0) {
    burnin_str <- ""
  } else {
    if (burnin < 1) {
      burnin_str <- paste0("_", burnin * 100, "burninPerc")
    } else {
      burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
    }
  }
  
  outdir_path <- gsub("/analyses", "/intermediate", paste0(getwd(), "/", gsub("^\\.", "", dirname(rep_paths[1]))))
  if (!dir.exists(outdir_path)) dir.create(outdir_path, recursive = TRUE)
  name_str <- gsub("_run(\\d+).*", "", basename(rep_paths[1]))
  path_str <- paste0(outdir_path, "/", "mean_", name_str, burnin_str, "_")
  
  subsample_freq <- as.integer(c(1000, 10000, 100000))
  rep_2rows <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F, nrows = 2))
  ngen_freq <- get_unique(sapply(rep_2rows, function(x) diff(x[1:2, 1])))
  subsample_freq <- subsample_freq[subsample_freq >= ngen_freq]
  if (length(subsample_freq) == 0) return()
  
  ngen_truncstr <- ifelse(ngen_trunc < ngen_freq, "", paste0(as.integer(ngen_trunc), "keepNgen", "_"))
  path_str <- paste0(path_str, ngen_truncstr)
  
  if (overwrite == FALSE) {
    out_paths <- sapply(subsample_freq, function(x) get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = x, add_str = ".tsv"))
    subsample_freq <- subsample_freq[which(file.exists(out_paths) == FALSE)]
  }
  nsets <- length(subsample_freq)
  if (nsets == 0) return()
  
  nreps <- length(rep_paths)
  rep_dfs <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F))
  
  nrow_logs <- sapply(rep_dfs, nrow)
  ngen_totals <- sapply(1:nreps, function(i) rep_dfs[[i]][nrow_logs[i], 1])
  
  if (burnin <= 0) {
    rep_pbs <- rep_dfs
  } else {
    if (burnin < 1) {
      ngen_burnins <- floor(ngen_totals * burnin) + 1L
    } else {
      ngen_burnins <- rep(as.integer(burnin), nreps)
    }
    
    if (ngen_trunc < ngen_freq) {
      ngen_truncs <- ngen_totals
    } else {
      ngen_truncs <- ngen_burnins + ngen_trunc
      if (any(ngen_truncs > ngen_totals)) stop("The specified number of generations to truncate is longer than the total number of generations.")
    }
    
    rep_pbs <- lapply(1:nreps, function(i) rep_dfs[[i]][rep_dfs[[i]][, 1] > ngen_burnins[i] & rep_dfs[[i]][, 1] <= ngen_truncs[i], -1, drop = FALSE])
  }
  
  rm(rep_dfs)
  gc()
  
  for (j in 1:nsets) {
    out_path <- get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = subsample_freq[j], add_str = ".tsv")
    if (file.exists(out_path) == FALSE || overwrite) {
      cat("working on ", out_path, "\n", sep = "")
      mean_mat <- get_mean(rep_pbs, ngen_freq = ngen_freq, subsample_freq = subsample_freq[j])
      if (!is.null(mean_mat)) {
        write.table(mean_mat, file = out_path, quote = FALSE, sep = "\t")
      }
    }
  }
}


get_quantiles <- function(rep_pbs, ngen_freq, subsample_freq = 0, qprobs = c(0.005, 0.025, 0.05, 0.125, 0.25, 0.5, 0.75, 0.875, 0.95, 0.975, 0.995)) {
  if (subsample_freq != ngen_freq && subsample_freq > 1) {
    if (subsample_freq < ngen_freq) {
      cat("subsampling frequency is not greater than the current sampling frequncy, so do nothing\n")
      return()
    } else if (subsample_freq %% ngen_freq != 0) {
      warning("subsampling frequency not multiple of sampled generation frequency")
      return()
    }
    
    by <- as.integer(subsample_freq / ngen_freq)
    rep_pbs <- lapply(rep_pbs, function(x) x[seq(from = by, to = nrow(x), by = by), , drop = FALSE])
  }
  
  nreps <- length(rep_pbs)
  nparams <- ncol(rep_pbs[[1]])
  param_stats <- vector("list", nreps + 1L)
  names(param_stats) <- c(paste0("rep", 1:nreps), "combined")
  for (i in 1:nreps) {
    param_stat <- vector("list", nparams)
    names(param_stat) <- names(rep_pbs[[i]])
    for (j in 1:nparams) {
      vec <- c(mean(rep_pbs[[i]][[j]]), quantile(rep_pbs[[i]][[j]], probs = qprobs))
      names(vec)[1] <- "mean"
      param_stat[[j]] <- vec
    }
    param_stats[[i]] <- param_stat
  }
  
  rep_combined <- do.call(rbind, rep_pbs)
  param_stat <- vector("list", nparams)
  names(param_stat) <- names(rep_combined)
  for (j in 1:nparams) {
    vec <- c(mean(rep_combined[[j]]), quantile(rep_combined[[j]], probs = qprobs))
    names(vec)[1] <- "mean"
    param_stat[[j]] <- vec
  }
  param_stats[[nreps + 1L]] <- param_stat
  
  return(param_stats)
}

generate_quantiles_all <- function(rep_paths, burnin = 2e8, ngen_trunc = 0, overwrite = FALSE) {
  
  if (burnin == 0) {
    burnin_str <- ""
  } else {
    if (burnin < 1) {
      burnin_str <- paste0("_", burnin * 100, "burninPerc")
    } else {
      burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
    }
  }
  
  outdir_path <- gsub("/analyses", "/intermediate", paste0(getwd(), "/", gsub("^\\.", "", dirname(rep_paths[1]))))
  if (!dir.exists(outdir_path)) dir.create(outdir_path, recursive = TRUE)
  name_str <- gsub("_run(\\d+).*", "", basename(rep_paths[1]))
  path_str <- paste0(outdir_path, "/", "quantiles_", name_str, burnin_str, "_")
  
  subsample_freq <- as.integer(c(1000, 10000, 100000))
  rep_2rows <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F, nrows = 2))
  ngen_freq <- get_unique(sapply(rep_2rows, function(x) diff(x[1:2, 1])))
  subsample_freq <- subsample_freq[subsample_freq >= ngen_freq]
  if (length(subsample_freq) == 0) return()
  
  ngen_truncstr <- ifelse(ngen_trunc < ngen_freq, "", paste0(as.integer(ngen_trunc), "keepNgen", "_"))
  path_str <- paste0(path_str, ngen_truncstr)
  
  if (overwrite == FALSE) {
    out_paths <- sapply(subsample_freq, function(x) get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = x, add_str = ".rds"))
    subsample_freq <- subsample_freq[which(file.exists(out_paths) == FALSE)]
  }
  nsets <- length(subsample_freq)
  if (nsets == 0) return()
  
  nreps <- length(rep_paths)
  rep_dfs <- lapply(rep_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F))
  
  nrow_logs <- sapply(rep_dfs, nrow)
  ngen_totals <- sapply(1:nreps, function(i) rep_dfs[[i]][nrow_logs[i], 1])
  
  if (burnin <= 0) {
    rep_pbs <- rep_dfs
  } else {
    if (burnin < 1) {
      ngen_burnins <- floor(ngen_totals * burnin) + 1L
    } else {
      ngen_burnins <- rep(as.integer(burnin), nreps)
    }
    
    if (ngen_trunc < ngen_freq) {
      ngen_truncs <- ngen_totals
    } else {
      ngen_truncs <- ngen_burnins + ngen_trunc
      if (any(ngen_truncs > ngen_totals)) stop("The specified number of generations to truncate is longer than the total number of generations.")
    }
    
    rep_pbs <- lapply(1:nreps, function(i) rep_dfs[[i]][rep_dfs[[i]][, 1] > ngen_burnins[i] & rep_dfs[[i]][, 1] <= ngen_truncs[i], -1, drop = FALSE])
  }
  
  rm(rep_dfs)
  gc()
  
  for (j in 1:nsets) {
    out_path <- get_outpath(path_str = path_str, ngen_freq = ngen_freq, subsample_freq = subsample_freq[j], add_str = ".rds")
    if (file.exists(out_path) == FALSE || overwrite) {
      cat("working on ", out_path, "\n", sep = "")
      param_stats <- get_quantiles(rep_pbs, ngen_freq = ngen_freq, subsample_freq = subsample_freq[j])
      if (!is.null(param_stats)) {
        saveRDS(param_stats, file = out_path)
      }
    }
  }
}


process_dir <- function(dirpath, filestr, stat = c("ess", "psrf", "mean", "quantiles"), burnin = 2e8, ngen_trunc = 0) {
  stat <- match.arg(stat, c("ess", "psrf", "mean", "quantiles"))
  
  file_paths <- list.files(dirpath, full.name = T, pattern = paste0("^", filestr, "_run\\d+.*\\.log$"))
  if (length(file_paths) == 0) return()
  run_ids_all <- as.integer(gsub(".*_run(\\d+).*", "\\1", file_paths))
  run_ids <- sort(unique(run_ids_all))
  file_strs <- sort(unique(gsub("_run(\\d+).*", "_run", file_paths)))
  nfs <- length(file_strs)
  if (nfs == 0) return()
  
  for (i in 1:nfs) {
    rep_paths <- grep(file_strs[i], file_paths, value = TRUE)
    cat("processing ", file_strs[i], "\n", sep = "")
    if (stat == "ess") {
      generate_ess_all(rep_paths, burnin = burnin, ngen_trunc = ngen_trunc)
    } else if (stat == "psrf") {
      generate_psrf_all(rep_paths, burnin = burnin, ngen_trunc = ngen_trunc)
    } else if (stat == "mean") {
      generate_mean_all(rep_paths, burnin = burnin, ngen_trunc = ngen_trunc)
    } else {
      generate_quantiles_all(rep_paths, burnin = burnin, ngen_trunc = ngen_trunc)
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)
dirpath <- args[1]
filestr <- args[2]
stat <- args[3]
ngen_trunc <- 0
if (length(args) >= 4 && is.finite(as.integer(args[4]))) ngen_trunc <- as.integer(args[4])
process_dir(dirpath = dirpath, filestr = filestr, stat = stat, ngen_trunc = ngen_trunc)

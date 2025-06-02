# this script contains functions that compute tree ESSs and tree PSRFs

get_dist_mat <- function(file_path) {
  dist_df <- read.table(file_path, header = T, sep = ",", check.names = F, stringsAsFactors = F)
  dist_mat <- as.matrix(dist_df[, -1])
  rownames(dist_mat) <- dist_df[, 1]
  rm(dist_df)
  gc()
  gc()
  if (nrow(dist_mat) != ncol(dist_mat)) stop("column and row numbers do not match")
  
  ntrees <- nrow(dist_mat)
  diag(dist_mat) <- 0L
  for (i in 1:ntrees) {
    idx <- which(is.na(dist_mat[i, ]))
    if (length(idx) > 0) {
      dist_mat[i, idx] <- dist_mat[idx, i]
    }
  }
  storage.mode(dist_mat) <- "integer"
  
  return(dist_mat)
}

get_pb_mats <- function(dist_mat, nsamp_pb, nsamp_burnin) {
  nsamp_total <- nrow(dist_mat)
  nreps <- max(length(nsamp_pb), length(nsamp_burnin))
  if (nreps > 1) {
    if (length(nsamp_burnin) == 1) nsamp_burnin <- rep(nsamp_burnin, nreps)
    if (length(nsamp_pb) == 1) nsamp_pb <- rep(nsamp_pb, nreps)
  }
  
  nsamp_end <- 0L
  mats <- vector("list", nreps)
  for (i in 1:nreps) {
    nsamp_start <- nsamp_end + nsamp_burnin[i] + 1L
    nsamp_end <- nsamp_start + nsamp_pb[i] - 1L
    if (nsamp_end > nsamp_total) nsamp_end <- nsamp_total
    mats[[i]] <- dist_mat[nsamp_start:nsamp_end, nsamp_start:nsamp_end]
  }
  
  return(mats)
}

get_pb_mat <- function(dist_mat, nsamp_pb, nsamp_burnin, inc_idx = NULL, nsamp_trunc = 0L) {
  nsamp_total <- nrow(dist_mat)
  nreps <- max(length(nsamp_pb), length(nsamp_burnin))
  if (nreps > 1) {
    if (length(nsamp_burnin) == 1) nsamp_burnin <- rep(nsamp_burnin, nreps)
    if (length(nsamp_pb) == 1) nsamp_pb <- rep(nsamp_pb, nreps)
  }
  nsamp_trunc <- min(nsamp_trunc, nsamp_pb)
  
  nsamp_end <- 0L
  nsamps_inc <- integer()
  for (i in 1:nreps) {
    nsamp_start <- nsamp_end + nsamp_burnin[i] + 1L
    nsamp_end <- nsamp_start + nsamp_pb[i] - 1L
    if (nsamp_trunc > 0) nsamp_end <- nsamp_start + nsamp_trunc - 1L
    if (nsamp_end > nsamp_total) nsamp_end <- nsamp_total
    if (nsamp_start <= nsamp_end && (is.null(inc_idx) || i %in% inc_idx)) nsamps_inc <- c(nsamps_inc, nsamp_start:nsamp_end)
    nsamp_end <- nsamp_start + nsamp_pb[i] - 1L
  }
  
  pb_mat <- dist_mat[nsamps_inc, nsamps_inc]
  return(pb_mat)
}


generate_treepsrf <- function(dir_path, dist_str = "RF_combined", ngen_freq = 100000L, ngen_total = 500000000L, burnin = 0.5, overwrite = FALSE) {
  
  if (burnin <= 0) {
    burnin_str <- ""
    burnin_str2 <- ""
    nsamp_burnin <- 0
  } else if (burnin < 1) {
    burnin_str <- paste0("_", burnin * 100, "burninPerc")
    burnin_str2 <- paste0("_", "burninPerc", burnin * 100)
    nsamp_burnin <- floor(ngen_total * burnin / ngen_freq) + 1L
  } else {
    burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
    burnin_str2 <- paste0("_", "burninNgen", as.integer(burnin))
    nsamp_burnin <- as.integer(burnin / ngen_freq + 1L)
  }
  
  nreps <- length(ngen_total)
  nsamp_rep <- floor(ngen_total / ngen_freq) + 1L
  nsamp_pb <- nsamp_rep - nsamp_burnin
  
  file_path <- list.files(dir_path, full.name = TRUE, pattern = paste0("pairwise", dist_str, "_resampfreq", as.integer(ngen_freq), burnin_str2, "\\.log$"))
  if (length(file_path) > 1) {
    stop("multiple files found")
  } else if (length(file_path) == 1) {
    if (burnin > 0) {
      nsamp_rep <- nsamp_pb
    }
    nsamp_burnin <- 0
  } else {
    file_path <- list.files(dir_path, full.name = TRUE, pattern = paste0("pairwise", dist_str, "_resampfreq", as.integer(ngen_freq), "\\.log$"))
    if (length(file_path) > 1) {
      stop("multiple files found")
    } else if (length(file_path) == 0) {
      return()
    }
  }
  
  cat("processing ", file_path, "\n", sep = "")
  outdir_path <- gsub("/analyses", "/intermediate", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  if (dir.exists(outdir_path) == FALSE) dir.create(outdir_path, recursive = TRUE, showWarnings = FALSE)
  
  out_path <- paste0(outdir_path, "/", "treepsrf_", dist_str, "_resampfreq", as.integer(ngen_freq), burnin_str, ".tsv")
  if (file.exists(out_path) && overwrite == FALSE) return()
  
  dist_mat <- get_dist_mat(file_path)
  pb_mats <- get_pb_mats(dist_mat, nsamp_pb, nsamp_burnin)
  
  if (length(unique(nsamp_pb)) > 1) {
    nsamp_max <- max(nsamp_pb)
    nsamp_min <- min(nsamp_pb)
    rep_idx_max <- which(nsamp_pb == nsamp_max)
    rep_idx_min <- which(nsamp_pb == nsamp_min)
    
    pb_mats_longruns <- pb_mats[rep_idx_max]
    pb_mat_longruns <- get_pb_mat(dist_mat, nsamp_pb, nsamp_burnin, inc_idx = rep_idx_max)
    
    pb_mats_shortruns <- pb_mats[rep_idx_min]
    pb_mat_shortruns <- get_pb_mat(dist_mat, nsamp_pb, nsamp_burnin, inc_idx = rep_idx_min)
    
    pb_mats_allruns <- lapply(pb_mats, function(mat) mat[1:nsamp_min, 1:nsamp_min])
    pb_mat_allruns <- get_pb_mat(dist_mat, nsamp_pb, nsamp_burnin, nsamp_trunc = nsamp_min)
    
    rm(dist_mat)
    gc()
    gc()
    
    W_longruns <- mean(sapply(pb_mats_longruns, function(mat) sum(mat^2)) / (nsamp_max * (nsamp_max - 1)))
    nm_longruns <- nrow(pb_mat_longruns)
    V_longruns <- sum(pb_mat_longruns^2) / (nm_longruns * (nm_longruns - 1))
    rhat_longruns <- sqrt(V_longruns / W_longruns)
    
    W_shortruns <- mean(sapply(pb_mats_shortruns, function(mat) sum(mat^2)) / (nsamp_min * (nsamp_min - 1)))
    nm_shortruns <- nrow(pb_mat_shortruns)
    V_shortruns <- sum(pb_mat_shortruns^2) / (nm_shortruns * (nm_shortruns - 1))
    rhat_shortruns <- sqrt(V_shortruns / W_shortruns)
    
    W_allruns <- mean(sapply(pb_mats_allruns, function(mat) sum(mat^2)) / (nsamp_min * (nsamp_min - 1)))
    nm_allruns <- nrow(pb_mat_allruns)
    V_allruns <- sum(pb_mat_allruns^2) / (nm_allruns * (nm_allruns - 1))
    rhat_allruns <- sqrt(V_allruns / W_allruns)
    
    rhat_mat <- matrix(c(rhat_longruns, W_longruns, V_longruns, rhat_shortruns, W_shortruns, V_shortruns, rhat_allruns, W_allruns, V_allruns), nrow = 1L)
    colnames(rhat_mat) <- c("Rhat_long", "W_long", "V_long", "Rhat_short", "W_short", "V_short", "Rhat_trunc", "W_trunc", "V_trunc")
    
  } else {
    pb_mat <- get_pb_mat(dist_mat, nsamp_pb, nsamp_burnin)
    rm(dist_mat)
    gc()
    gc()
    
    nm <- nrow(pb_mat)
    W <- mean(sapply(pb_mats, function(mat) sum(mat^2)) / (nsamp_pb[1] * (nsamp_pb[1] - 1)))
    V <- sum(pb_mat^2) / (nm * (nm - 1))
    rhat <- sqrt(V / W)
    
    rhat_mat <- matrix(c(rhat, W, V), nrow = 1L)
    colnames(rhat_mat) <- c("Rhat", "W", "V")
  }
  
  write.table(rhat_mat, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE)
}


get_mds_mats <- function(mat, nsamp_pb, nsamp_burnin) {
  nsamp_total <- nrow(mat)
  nreps <- max(length(nsamp_pb), length(nsamp_burnin))
  if (nreps > 1) {
    if (length(nsamp_burnin) == 1) nsamp_burnin <- rep(nsamp_burnin, nreps)
    if (length(nsamp_pb) == 1) nsamp_pb <- rep(nsamp_pb, nreps)
  }
  
  nsamp_end <- 0L
  mats <- vector("list", nreps)
  for (i in 1:nreps) {
    nsamp_start <- nsamp_end + nsamp_burnin[i] + 1L
    nsamp_end <- nsamp_start + nsamp_pb[i] - 1L
    if (nsamp_end > nsamp_total) nsamp_end <- nsamp_total
    mats[[i]] <- mat[nsamp_start:nsamp_end, , drop = FALSE]
  }
  
  return(mats)
}

get_mds_mat <- function(mds_path, mds_method = "cmds") {
  mds <- readRDS(mds_path)
  if (mds_method %in% c("isomds", "monomds")) {
    mds_mat <- mds[[mds_method]]$points
  } else {
    if (length(mds[[1]]) <= 3 && mds_method %in% names(mds[[1]])) {
      mds_mat <- mds[[1]][[mds_method]]
      if (mds_method == "metmds") mds_mat <- mds_mat$conf
    } else {
      mds_mat <- mds[[mds_method]]
      if (mds_method == "metmds") mds_mat <- mds_mat$conf
    }
  }
  return(mds_mat)
}


calc_ess_mat <- function(pb_mat) {
  fs <- c(treess:::frechetCorrelationESS, treess:::minPseudoESS, treess:::medianPseudoESS, treess:::approximateESS)
  esss <- sapply(fs, function(f) f(pb_mat), USE.NAMES = FALSE)
  ess_mat <- matrix(esss, nrow = 1L, ncol = 4L)
  colnames(ess_mat) <- c("frechetCorrelationESS", "minPseudoESS", "medianPseudoESS", "approximateESS")
  return(ess_mat)
}

calc_mdsess_mat <- function(mds_mat) {
  ndim <- ncol(mds_mat)
  mcmcse_multiess <- mcmcse::multiESS(mds_mat)
  mcmcse_ess <- sapply(1:ndim, function(j) mcmcse::ess(mds_mat[, j]), USE.NAMES = FALSE)
  coda_ess <- as.numeric(sapply(1:ndim, function(j) coda::effectiveSize(mds_mat[, j]), USE.NAMES = FALSE))
  tracer_ess <- sapply(1:ndim, function(j) convenience:::essTracerC(mds_mat[, j], max_lag = 2000L), USE.NAMES = FALSE)
  rstan_ess <- sapply(1:ndim, function(j) rstan::ess_bulk(mds_mat[, j]), USE.NAMES = FALSE)
  esss <- c(mcmcse_multiess, mcmcse_ess, coda_ess, tracer_ess, rstan_ess)
  mdsess_mat <- matrix(esss, nrow = 1L, ncol = length(esss))
  colnames(mdsess_mat) <- c(paste0("mds", ndim, "DMcmcseMultiESS"), paste0("mdsD", 1:ndim, "McmcseESS"), paste0("mdsD", 1:ndim, "CodaESS"), 
                            paste0("mdsD", 1:ndim, "TracerLag2000ESS"), paste0("mdsD", 1:ndim, "RstanESS"))
  return(mdsess_mat)
}


generate_treess <- function(dir_path, dist_str = "RF_combined", ngen_freq = 100000L, ngen_total = 500000000L, burnin = 0.5, overwrite = FALSE) {
  
  if (burnin <= 0) {
    burnin_str <- ""
    burnin_str2 <- ""
    nsamp_burnin <- 0
  } else if (burnin < 1) {
    burnin_str <- paste0("_", burnin * 100, "burninPerc")
    burnin_str2 <- paste0("_", "burninPerc", burnin * 100)
    nsamp_burnin <- floor(ngen_total * burnin / ngen_freq) + 1L
  } else {
    burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
    burnin_str2 <- paste0("_", "burninNgen", as.integer(burnin))
    nsamp_burnin <- as.integer(burnin / ngen_freq + 1L)
  }
  
  nreps <- length(ngen_total)
  nsamp_rep <- floor(ngen_total / ngen_freq) + 1L
  nsamp_pb <- nsamp_rep - nsamp_burnin
  
  file_path <- list.files(dir_path, full.name = T, pattern = paste0("pairwise", dist_str, "_resampfreq", as.integer(ngen_freq), burnin_str2, "\\.log$"))
  if (length(file_path) > 1) {
    stop("multiple files found")
  } else if (length(file_path) == 1) {
    if (burnin > 0) {
      nsamp_rep <- nsamp_pb
    }
    nsamp_burnin <- 0
  } else {
    file_path <- list.files(dir_path, full.name = TRUE, pattern = paste0("pairwise", dist_str, "_resampfreq", as.integer(ngen_freq), "\\.log$"))
    if (length(file_path) > 1) {
      stop("multiple files found")
    } else if (length(file_path) == 0) {
      return()
    }
  }
  
  cat("processing ", file_path, "\n", sep = "")
  outdir_path <- gsub("/analyses", "/intermediate", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  if (dir.exists(outdir_path) == FALSE) dir.create(outdir_path, recursive = TRUE, showWarnings = FALSE)
  
  out_path <- paste0(outdir_path, "/", "treess_", dist_str, "_resampfreq", as.integer(ngen_freq), burnin_str, ".tsv")
  if (file.exists(out_path) && (overwrite == FALSE)) return()
  
  dist_mat <- get_dist_mat(file_path)
  pb_mats <- get_pb_mats(dist_mat, nsamp_pb, nsamp_burnin)
  ess_mats <- lapply(pb_mats, calc_ess_mat)
  ess_mat <- do.call(rbind, ess_mats)
  rownames(ess_mat) <- paste0("rep", 1:nreps)
  
  if (length(unique(nsamp_pb)) > 1) {
    nsamp_max <- max(nsamp_pb)
    nsamp_min <- min(nsamp_pb)
    rep_idx_max <- which(nsamp_pb == nsamp_max)
    
    pb_mats_trunc <- lapply(pb_mats[rep_idx_max], function(mat) mat[1:nsamp_min, 1:nsamp_min])
    ess_mats_trunc <- lapply(pb_mats_trunc, calc_ess_mat)
    ess_mat_trunc <- do.call(rbind, ess_mats_trunc)
    rownames(ess_mat_trunc) <- paste0("rep", rep_idx_max, "trunc")
    
    ess_mat <- rbind(ess_mat, ess_mat_trunc)
  }
  
  # mds
  ndims <- c(2, 10)
  mds_method <- "cmds"
  file_strs <- paste0(dist_str, "_resampfreq", as.integer(ngen_freq), burnin_str, "_", mds_method, ndims, "d")
  nfiles <- length(file_strs)
  for (j in 1:nfiles) {
    mds_path <- list.files(outdir_path, full.name = T, pattern = paste0(file_strs[j], "\\.rds$"))
    if (length(mds_path) != 1) {
      message("mds file missing so abort")
      return()
    }
    mds_mat <- get_mds_mat(mds_path, mds_method)
    mds_mats <- get_mds_mats(mds_mat, nsamp_pb, nsamp_burnin)
    
    if (j == 1) {
      mdsess_mat <- do.call(rbind, lapply(mds_mats, calc_mdsess_mat))
    } else if (j == 2) {
      ndim <- ncol(mds_mats[[1]])
      mdsess_ndims <- sapply(mds_mats, function(x) max(sapply(1:ndim, function(n) mcmcse::multiESS(x[, 1:n, drop = FALSE]))), USE.NAMES = FALSE)
      mdsess_mat <- cbind(mdsess_ndims, mdsess_mat)
      colnames(mdsess_mat)[1] <- paste0("mds", ndim, "DMaxMcmcseMultiESS")
    }
  }
  rownames(mdsess_mat) <- paste0("rep", 1:nreps)
  
  if (length(unique(nsamp_pb)) > 1) {
    nsamp_max <- max(nsamp_pb)
    nsamp_min <- min(nsamp_pb)
    rep_idx_max <- which(nsamp_pb == nsamp_max)
    
    for (j in 1:nfiles) {
      mds_path <- list.files(outdir_path, full.name = T, pattern = paste0(file_strs[j], "\\.rds$"))
      if (length(mds_path) != 1) {
        message("mds file missing so abort")
        return()
      }
      mds_mat <- get_mds_mat(mds_path, mds_method)
      mds_mats <- get_mds_mats(mds_mat, nsamp_pb, nsamp_burnin)
      mds_mats_trunc <- lapply(mds_mats[rep_idx_max], function(mat) mat[1:nsamp_min, , drop = FALSE])
      
      if (j == 1) {
        mdsess_mat_trunc <- do.call(rbind, lapply(mds_mats_trunc, calc_mdsess_mat))
      } else if (j == 2) {
        ndim <- ncol(mds_mats_trunc[[1]])
        mdsess_ndims <- sapply(mds_mats_trunc, function(x) max(sapply(1:ndim, function(n) mcmcse::multiESS(x[, 1:n, drop = FALSE]))), USE.NAMES = FALSE)
        mdsess_mat_trunc <- cbind(mdsess_ndims, mdsess_mat_trunc)
        colnames(mdsess_mat_trunc)[1] <- paste0("mds", ndim, "DMaxMcmcseMultiESS")
      }
    }
    
    rownames(mdsess_mat_trunc) <- paste0("rep", rep_idx_max, "trunc")
    mdsess_mat <- rbind(mdsess_mat, mdsess_mat_trunc)
  }
  
  colnames(mdsess_mat) <- gsub("^mds", mds_method, colnames(mdsess_mat))
  mat <- cbind(ess_mat, mdsess_mat)
  colnames(mat) <- paste0(colnames(mat), "_", gsub("_.*", "", dist_str))
  
  write.table(mat, file = out_path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}


process_dir <- function(dir_path, stat = c("ess", "psrf"), dist_str = "RF_combined", ngen_freq = 100000L, ngen_total = NULL, burnin = 0.5, overwrite = FALSE) {
  stat <- match.arg(stat, c("ess", "psrf"))
  
  if (burnin <= 0) {
    burnin_str <- ""
  } else if (burnin < 1) {
    burnin_str <- paste0("_", burnin * 100, "burninPerc")
  } else {
    burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
  }
  
  if (is.null(ngen_total)) {
    trees_paths <- list.files(dir_path, pattern = paste0("^run\\d+", ".trees"))
    if (length(trees_paths) < 1) stop("cannot find the corresponding trees file")
    ngen_str <- "^tree STATE_"
    rpl_str <- "^tree STATE_(\\d+)"
    ngen_total <- sapply(trees_paths, function(trees_path) as.integer(gsub(paste0(rpl_str, ".*"), "\\1", system(paste0('tac ', trees_path, ' | grep -m1 "', ngen_str, '"'), intern = TRUE))))
  }

  cat("processing tree", stat, "for", dir_path, "\n", sep = "")
  if (stat == "ess") {
    generate_treess(dir_path, dist_str = dist_str, ngen_freq = ngen_freq, ngen_total = ngen_total, burnin = burnin, overwrite = overwrite)
  } else {
    generate_treepsrf(dir_path, dist_str = dist_str, ngen_freq = ngen_freq, ngen_total = ngen_total, burnin = burnin, overwrite = overwrite)
  }
}


args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
ngen_freq <- as.integer(args[2])
dist_str <- args[3]
stat <- args[4]

process_dir(dir_path = dir_path, stat = stat, dist_str = dist_str, ngen_freq = ngen_freq, burnin = 2e8, overwrite = FALSE)

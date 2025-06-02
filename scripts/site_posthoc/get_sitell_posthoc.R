# this script contains functions that compute the log likelihood each site in the alignment given the tree and substitution model parameters
library(phangorn)
source("scripts/beasttree_readwrite/tree_proc_functions.R")


calcwrite_phyll <- function(sitell_path) {
  dir_path <- dirname(sitell_path)
  out_name <- gsub("sitell_posthoc_", "phyll_posthoc_", basename(sitell_path))
  out_path <- paste0(dir_path, "/", out_name)
  
  datadir_path <- gsub("/analyses", "/data", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  ali_paths <- list.files(datadir_path, full.names = TRUE, pattern = "\\.fasta$")
  if (length(ali_paths) > 1) {
    ali_paths <- ali_paths[grep("ali\\.fasta$", basename(ali_paths), invert = T)]
  }
  alis <- lapply(ali_paths, function(x) ape::read.dna(x, format = "fasta", as.character = T))
  
  sp_path <- list.files(datadir_path, full.names = TRUE, pattern = "^site_patterns_notuseambi\\.tsv$")
  sp_df <- read.table(sp_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  
  sitell_df <- read.table(sitell_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  lls <- apply(sitell_df[, -1], 1, function(x) sum(x * sp_df$count))
  ll_df <- data.frame(sitell_df[, 1], lls)
  colnames(ll_df) <- c("state", "likelihood")
  write.table(ll_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE)
}

calcwrite_sitell <- function(trees_path, param_path, tipidx_toprune = NULL, ncores = 1L, overwrite = FALSE) {
  
  dir_path <- dirname(trees_path)
  out_name <- gsub("param_", "sitell_posthoc_", basename(param_path))
  if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) out_name <- paste0("tip", paste(tipidx_toprune, collapse = "and"), "pruned", "_", out_name)
  
  while (nchar(out_name) > 255 && grepl("and[0-9]+(?=pruned)", out_name, perl = TRUE)) {
    out_name <- sub("and[0-9]+(?=pruned)", "", out_name, perl = TRUE)
  }
  out_path <- paste0(dir_path, "/", out_name)
  
  ngen_burnin <- 0L
  if (grepl("burninNgen", trees_path)) ngen_burnin <- as.integer(gsub(".*burninNgen(\\d+).*", "\\1", trees_path))
  
  datadir_path <- gsub("/analyses", "/data", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  ali_paths <- list.files(datadir_path, full.names = TRUE, pattern = "\\.fasta$")
  if (length(ali_paths) > 1) {
    ali_paths <- ali_paths[grep("ali\\.fasta$", basename(ali_paths), invert = T)]
  }
  alis <- lapply(ali_paths, function(x) ape::read.dna(x, format = "fasta", as.character = T))
  
  sp_path <- list.files(datadir_path, full.names = TRUE, pattern = "^site_patterns_notuseambi\\.tsv$")
  sp_df <- read.table(sp_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  part_names <- unique(sp_df$partition)
  npart <- length(part_names)
  sps <- vector("list", npart)
  for (j in 1:npart) {
    sps[[j]] <- do.call(cbind, strsplit(sp_df[sp_df$partition == part_names[j], ]$pattern, ""))
    rownames(sps[[j]]) <- rownames(alis[[j]])
    if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) {
      if (any(tipidx_toprune < 0 | tipidx_toprune > nrow(sps[[j]]))) stop("tip indices to prune are misspecified")
      sps[[j]] <- sps[[j]][-tipidx_toprune, , drop = FALSE]
    }
  }
  
  nsps <- sapply(sps, ncol)
  spdats <- lapply(sps, phangorn::phyDat, compress = FALSE)
  for (j in 1:npart) {
    nsp_this <- unique(lengths(spdats[[j]]))
    if (length(nsp_this) != 1 || nsp_this != nsps[j]) stop("wrong site pattern number")
  }
  nsp_total <- sum(nsps)
  ntips <- nrow(sps[[1]])
  tiplabs <- rownames(sps[[1]])
  
  ntrees <- as.integer(system(paste0("grep -c 'tree STATE_' ", trees_path), intern = TRUE))
  ntrees_done <- 0L
  if (file.exists(out_path) && overwrite == FALSE) {
    nlines <- as.integer(system(paste0("wc -l < ", out_path), intern = TRUE))
    if (nlines == ntrees + 1L) return()
    
    ncol_tar <- nsp_total + 1L
    row_last <- system(paste0("tail -n1 ", out_path), intern = TRUE)
    sep_char <- "\t"
    row_last_nseps <- nchar(row_last) - nchar(gsub(sep_char, "", row_last))
    row_last_lastchar <- system(paste0("tail -c 1 ", out_path), intern = TRUE)
    if (row_last_nseps != ncol_tar - 1 || (length(row_last_lastchar) == 1 && row_last_lastchar != "")) {
      out_path_tmp <- paste0(out_path, "tmp")
      system(paste0("head -n -1 ", out_path, " > ", out_path_tmp, " && mv ", out_path_tmp, " ", out_path))
      nlines <- nlines - 1L
    }
    
    ntrees_done <- max(0L, nlines - 1L)
  }

  # read trees
  trees <- read.nexus.beast(trees_path)
  tiplabs_tree <- trees[[1]]$tip.label
  if (ntips != length(tiplabs_tree) - length(tipidx_toprune)) stop("tip numbers do not match")
  if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) {
    tiplabs_tree <- tiplabs_tree[-tipidx_toprune]
  }
  if (identical(tiplabs, tiplabs_tree) == FALSE) stop("tip labels do not match")
  if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) trees <- ape::drop.tip(trees, tipidx_toprune)
  
  ngen_trees <- as.integer(gsub("STATE_", "", names(trees)))
  if (ntrees != length(trees)) stop("tree number do not match")
  if (ntrees_done > 0) {
    trees <- trees[(ntrees_done + 1L):ntrees]
    ngen_trees <- ngen_trees[(ntrees_done + 1L):ntrees]
    ntrees <- length(trees)
    gc()
  }
  ngen_trees <- ngen_trees + ngen_burnin
  
  # read in log files
  firstrow <- read.table(param_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F, nrows = 1)
  
  f1idx <- grep("frequencies1$", colnames(firstrow))
  if (length(f1idx) == 0) {
    stop("cannot find substitution model parameter")
  } else if (length(f1idx) != npart) {
    stop("partition number do not match")
  }
  if (npart > 1 && identical(part_names, gsub("\\..*", "", colnames(firstrow)[f1idx])) == FALSE) stop("partition names do not match")
  
  idx_all <- vector("list", npart)
  for (j in 1:npart) {
    part_str <- ifelse(npart > 1, paste0(part_names[j], "."), "")
    freq_str <- paste0(part_str, "frequencies\\d{1}$")
    freq_idx <- grep(freq_str, colnames(firstrow))
    
    idx <- list(freq_idx = freq_idx)
    
    kappa_str <- paste0(part_str, "kappa")
    kappa_idx <- which(colnames(firstrow) == kappa_str)
    idx <- c(idx, list(kappa_idx = kappa_idx))
    if (length(kappa_idx) == 0) {
      gtr_str <- paste0(part_str, "gtr.rates.rate")
      gtr_idx <- grep(gtr_str, colnames(firstrow))
      if (length(gtr_idx) != 6) stop("cannot find either kappa or gtr")
      idx <- c(idx, list(gtr_idx = gtr_idx))
    }
    
    alpha_str <- paste0(part_str, "alpha")
    alpha_idx <- which(colnames(firstrow) == alpha_str)
    idx <- c(idx, list(alpha_idx = alpha_idx))
    
    pinv_str <- paste0(part_str, "pInv")
    pinv_idx <- which(colnames(firstrow) == pinv_str)
    idx <- c(idx, list(pinv_idx = pinv_idx))
    
    if (npart > 1) {
      mu_str <- paste0(part_str, "mu")
      mu_idx <- which(colnames(firstrow) == mu_str)
      idx <- c(idx, list(mu_idx = mu_idx))
    }
    
    idx_all[[j]] <- idx 
  }
  
  param_df <- read.table(param_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  param_df <- param_df[param_df[, 1] %in% ngen_trees, , drop = FALSE]
  nrow_log <- nrow(param_df)
  dum <- rep(1, nrow_log)
  
  res_all <- vector("list", npart)
  for (j in 1:npart) {
    idx <- idx_all[[j]]
    bfs <- param_df[, idx$freq_idx]
    
    if (length(idx$alpha_idx) > 0) {
      alphas <- param_df[, idx$alpha_idx]
      site.rate <- "gamma"
      ngcat <- 4
    } else {
      alphas <- dum
      site.rate <- "free_rate"
      ngcat <- 1
    }
    
    if (length(idx$pinv_idx) > 0) {
      pinvs <- param_df[, idx$pinv_idx]
    } else {
      pinvs <- rep(0, nrow_log)
    }
    
    if ("mu_idx" %in% names(idx) && length(idx$mu_idx) > 0) {
      mus <- param_df[, idx$mu_idx]
    } else {
      mus <- dum
    }
    
    if (length(idx$kappa_idx) > 0) {
      kappas <- param_df[, idx$kappa_idx]
      ers <- cbind(dum, kappas, dum, dum, kappas, dum)
    } else {
      ers <- param_df[, idx$gtr_idx]
    }
    
    res <- list(bfs = bfs, alphas = alphas, pinvs = pinvs, mus = mus, ers = ers, site.rate = site.rate, ngcat = ngcat)
    res_all[[j]] <- res
  }
  
  if (ncores > 1) {
    
    cl <- parallel::makeCluster(ncores, outfile = "")
    parallel::clusterEvalQ(cl, {library(phangorn)})
    # parallel::clusterExport(cl = cl, varlist = c("trees", "spdats", "res_all", "npart", "pml_custom", "pml.fit_custom"), envir = environment())
    parallel::clusterExport(cl = cl, varlist = c("trees", "spdats", "res_all", "npart"), envir = environment())
    sid <- 1L
    ncores_multi <- 40L
    
    while (sid <= ntrees) {
      eid <- sid + ncores * ncores_multi - 1L
      if (eid > ntrees) eid <- ntrees
      
      sitell_mat <- parallel::parSapply(cl, sid:eid, function(l) {
        tree <- trees[[l]]
        
        sitell_this <- vector("list", npart)
        for (j in 1:npart) {
          spdat <- spdats[[j]]
          bf <- as.numeric(res_all[[j]]$bfs[l, ])
          Q <- as.numeric(res_all[[j]]$ers[l, ])
          inv <- res_all[[j]]$pinvs[l]
          k <- res_all[[j]]$ngcat
          shape <- res_all[[j]]$alphas[l]
          rate <- res_all[[j]]$mus[l]
          site.rate <- res_all[[j]]$site.rate
          # res <- pml_custom(tree = tree, data = spdat, bf = bf, Q = Q, inv = inv, k = k, shape = shape, rate = rate, site.rate = site.rate, geps = 1e-6)
          res <- phangorn::pml(tree = tree, data = spdat, bf = bf, Q = Q, inv = inv, k = k, shape = shape, rate = rate, site.rate = site.rate)
          sitell_this[[j]] <- res$siteLik
        }
        
        sitells <- do.call(c, sitell_this)
        return(sitells)
      })
      
      out_df <- data.frame(ngen_trees[sid:eid], t(sitell_mat))
      out_colnames <- c("state", paste0("siteLikelihood_", 1:nsp_total))
      colnames(out_df) <- out_colnames
      first_write <- ntrees_done + sid == 1
      write.table(out_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
      
      cat("processed gen ", ngen_trees[eid], "\n")
      sid <- eid + 1L
    }
    
    parallel::stopCluster(cl)
  } else {
    
    sitell_colnames <- c("state", paste0("siteLikelihood_", 1:nsp_total))
    niter_out <- 100L
    sitell_mat <- matrix(ncol = nsp_total, nrow = niter_out)
    ntrees_processed <- 1L
    write_now <- FALSE
    for (l in 1:ntrees) {
      tree <- trees[[l]]
      
      sitell_this <- vector("list", npart)
      for (j in 1:npart) {
        spdat <- spdats[[j]]
        bf <- as.numeric(res_all[[j]]$bfs[l, ])
        Q <- as.numeric(res_all[[j]]$ers[l, ])
        inv <- res_all[[j]]$pinvs[l]
        k <- res_all[[j]]$ngcat
        shape <- res_all[[j]]$alphas[l]
        rate <- res_all[[j]]$mus[l]
        site.rate <- res_all[[j]]$site.rate
        # res <- pml_custom(tree = tree, data = spdat, bf = bf, Q = Q, inv = inv, k = k, shape = shape, rate = rate, site.rate = site.rate, geps = 1e-6)
        res <- phangorn::pml(tree = tree, data = spdat, bf = bf, Q = Q, inv = inv, k = k, shape = shape, rate = rate, site.rate = site.rate)
        sitell_this[[j]] <- res$siteLik
      }
      sitells <- do.call(c, sitell_this)
      
      if (l %% niter_out == 0) {
        sitell_mat[niter_out, ] <- sitells
        write_now <- TRUE
      } else {
        sitell_mat[l %% niter_out, ] <- sitells
        if (l == ntrees) {
          sitell_mat <- sitell_mat[1:(l %% niter_out), , drop = FALSE]
          write_now <- TRUE
        }
      }

      if (write_now) {
        sitell_df <- data.frame(ngen_trees[ntrees_processed:l], formatC(sitell_mat, digits = 3, format = "f"))
        colnames(sitell_df) <- sitell_colnames
        first_write <- ntrees_done + l == niter_out
        write.table(sitell_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)

        ntrees_processed <- l + 1
        cat("processed gen ", ngen_trees[l], "\n")
        sitell_mat <- matrix(ncol = nsp_total, nrow = niter_out)
        write_now <- FALSE
      }
    }
    
  }
}

generate_sitell <- function(dir_path, rep_id, tipidx_toprune = NULL, ncores = 1L, overwrite = FALSE) {
  ridstr <- formatC(rep_id, width = 2, format = "d", flag = 0)
  trees_path <- list.files(dir_path, full.name = T, pattern = paste0("^bltree_run", ridstr, "\\.trees$"))
  if (length(trees_path) != 1) stop("cannot find tree file")
  
  param_path <- list.files(dir_path, full.name = T, pattern = paste0("^param_run", ridstr, "\\.log$"))
  if (length(param_path) != 1) stop("cannot find param file")
  
  calcwrite_sitell(trees_path, param_path, tipidx_toprune, ncores, overwrite)
  
  dir_path <- dirname(trees_path)
  out_name <- gsub("param_", "sitell_posthoc_", basename(param_path))
  if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) out_name <- paste0("tip", paste(tipidx_toprune, collapse = "and"), "pruned", "_", out_name)
  
  while (nchar(out_name) > 255 && grepl("and[0-9]+(?=pruned)", out_name, perl = TRUE)) {
    out_name <- sub("and[0-9]+(?=pruned)", "", out_name, perl = TRUE)
  }
  out_path <- paste0(dir_path, "/", out_name)
  
  out_path_ll <- gsub("sitell_", "phyll_", out_path)
  if (file.exists(out_path) && (file.exists(out_path_ll) == FALSE || overwrite)) calcwrite_phyll(out_path)
}

args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
rep_id <- as.integer(args[2])
tipidx_toprune <- NULL
if (length(args) >= 3 && args[3] != "") {
  tipidx_toprune <- as.integer(gsub(" ", "", unlist(strsplit(args[3], "and"))))
  if (any(is.na(tipidx_toprune))) stop("tip indices to prune are misspecified")
  if (length(tipidx_toprune) == 0) tipidx_toprune <- NULL
}
ncores <- 1L
if (length(args) >= 4) {
  ncores <- as.integer(args[4])
}
generate_sitell(dir_path = dir_path, rep_id = rep_id, tipidx_toprune = tipidx_toprune, ncores = ncores, overwrite = FALSE)

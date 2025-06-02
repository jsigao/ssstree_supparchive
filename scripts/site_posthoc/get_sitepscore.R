# this script contains functions that compute the parsimony score of each site in the alignment
library(phangorn)
source('scripts/beasttree_readwrite/tree_proc_functions.R')

calcwrite_sitepscore <- function(trees_path, tipidx_toprune = NULL, ncores = 1L, overwrite = FALSE) {
  
  dir_path <- dirname(trees_path)
  ridstr <- gsub("^run(\\d+)\\.trees$", "\\1", basename(trees_path))
  out_name <- paste0("sitepscore", "_run", ridstr, ".log")
  if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) out_name <- paste0("tip", paste(tipidx_toprune, collapse = "and"), "pruned", "_", out_name)
  
  while (nchar(out_name) > 255 && grepl("and[0-9]+(?=pruned)", out_name, perl = TRUE)) {
    out_name <- sub("and[0-9]+(?=pruned)", "", out_name, perl = TRUE)
  }
  out_path <- paste0(dir_path, "/", out_name)
  
  ngen_burnin <- 0L
  if (grepl("burninNgen", trees_path)) ngen_burnin <- as.integer(gsub(".*burninNgen(\\d+).*", "\\1", trees_path))
  
  datadir_path <- gsub("/analyses", "/data", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  ali_path <- list.files(datadir_path, full.names = TRUE, pattern = "^ali\\.fasta$")
  ali <- ape::read.dna(ali_path, format = "fasta", as.character = T)
  
  sp_path <- list.files(datadir_path, full.names = TRUE, pattern = "^site_patterns_notuseambi\\.tsv$")
  sp_df <- read.table(sp_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  spmat <- do.call(cbind, strsplit(sp_df$pattern, ""))
  rownames(spmat) <- rownames(ali)
  
  if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) {
    if (any(tipidx_toprune < 0 | tipidx_toprune > nrow(spmat))) stop("tip indices to prune are misspecified")
    spmat <- spmat[-tipidx_toprune, , drop = FALSE]
  }
  
  nsp_total <- ncol(spmat)
  ntips <- nrow(spmat)
  tiplabs <- rownames(spmat)
  
  spdat <- phangorn::phyDat(spmat, compress = FALSE)
  nsp_this <- unique(lengths(spdat))
  if (length(nsp_this) != 1 || nsp_this != nsp_total) stop("wrong site pattern number")
  
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
  
  if (ncores > 1) {
    
    cl <- parallel::makeCluster(ncores, outfile = "")
    parallel::clusterEvalQ(cl, {library(phangorn)})
    parallel::clusterExport(cl = cl, varlist = c("trees", "spdat"), envir = environment())
    sid <- 1L
    ncores_multi <- 40L
    
    while (sid <= ntrees) {
      eid <- sid + ncores * ncores_multi - 1L
      if (eid > ntrees) eid <- ntrees
      
      pscore_mat <- parallel::parSapply(cl, sid:eid, function(l) {
        tree <- trees[[l]]
        pscores <- phangorn::parsimony(tree = tree, data = spdat, method = "sankoff", site = "site")
        return(pscores)
      })
      
      storage.mode(pscore_mat) <- "integer"
      out_df <- data.frame(ngen_trees[sid:eid], t(pscore_mat))
      out_colnames <- c("state", paste0("pscore_", 1:nsp_total))
      colnames(out_df) <- out_colnames
      first_write <- ntrees_done + sid == 1
      write.table(out_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
      
      cat("processed gen ", ngen_trees[eid], "\n")
      sid <- eid + 1L
    }
    
    parallel::stopCluster(cl)
  } else {
    pscore_colnames <- c("state", paste0("pscore_", 1:nsp_total))
    niter_out <- 1000L
    pscore_mat <- matrix(ncol = nsp_total, nrow = niter_out)
    ntrees_processed <- 1L
    write_now <- FALSE
    for (l in 1:ntrees) {
      tree <- trees[[l]]
      pscores <- phangorn::parsimony(tree = tree, data = spdat, method = "sankoff", site = "site")
      
      if (l %% niter_out == 0) {
        pscore_mat[niter_out, ] <- pscores
        write_now <- TRUE
      } else {
        pscore_mat[l %% niter_out, ] <- pscores
        if (l == ntrees) {
          pscore_mat <- pscore_mat[1:(l %% niter_out), , drop = FALSE]
          write_now <- TRUE
        }
      }
      
      if (write_now) {
        storage.mode(pscore_mat) <- "integer"
        pscore_df <- data.frame(ngen_trees[ntrees_processed:l], pscore_mat)
        colnames(pscore_df) <- pscore_colnames
        first_write <- ntrees_done + l == niter_out
        write.table(pscore_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
        
        pscore_mat <- matrix(ncol = nsp_total, nrow = niter_out)
        ntrees_processed <- l + 1
        cat("processed gen ", ngen_trees[l], "\n")
        write_now <- FALSE
      }
    }
  }
  
}


calcwrite_alipscore <- function(sitepscore_path) {
  dir_path <- dirname(sitepscore_path)
  out_path <- gsub("sitepscore_", "alipscore_", sitepscore_path)
  
  datadir_path <- gsub("/analyses", "/data", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  sp_path <- list.files(datadir_path, full.names = TRUE, pattern = "^site_patterns_notuseambi\\.tsv$")
  if (length(sp_path) != 1) stop("cannot site pattern file")
  sp_df <- read.table(sp_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  
  sitepscore_df <- read.table(sitepscore_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  pscores <- as.integer(apply(sitepscore_df[, -1], 1, function(x) sum(x * sp_df$count)))
  
  alipscore_df <- data.frame(state = sitepscore_df[, 1], pscore = pscores, stringsAsFactors = FALSE)
  write.table(alipscore_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE)
}


generate_sitepscore <- function(dir_path, rep_id, tipidx_toprune = NULL, ncores = 1L, overwrite = FALSE) {
  ridstr <- formatC(rep_id, width = 2, format = "d", flag = 0)
  trees_path <- list.files(dir_path, full.name = T, pattern = paste0("^run", ridstr, "\\.trees$"))
  if (length(trees_path) != 1) stop("cannot find tree file")
  
  calcwrite_sitepscore(trees_path, tipidx_toprune, ncores, overwrite)
  
  out_name <- paste0("sitepscore", "_run", ridstr, ".log")
  if (is.null(tipidx_toprune) == FALSE && length(tipidx_toprune) > 0) out_name <- paste0("tip", paste(tipidx_toprune, collapse = "and"), "pruned", "_", out_name)
  
  while (nchar(out_name) > 255 && grepl("and[0-9]+(?=pruned)", out_name, perl = TRUE)) {
    out_name <- sub("and[0-9]+(?=pruned)", "", out_name, perl = TRUE)
  }
  out_path <- paste0(dir_path, "/", out_name)
  
  ps_path <- gsub("sitepscore_", "alipscore_", out_path)
  if (file.exists(out_path) && (file.exists(ps_path) == FALSE || overwrite)) calcwrite_alipscore(out_path)
}

args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
rep_id <- as.integer(args[2])
tipidx_toprune <- NULL
if (length(args) >= 3 || args[3] != "") {
  tipidx_toprune <- as.integer(gsub(" ", "", unlist(strsplit(args[3], "and"))))
  if (any(is.na(tipidx_toprune))) stop("tip indices to prune are misspecified")
  if (length(tipidx_toprune) == 0) tipidx_toprune <- NULL
}
ncores <- 1L
if (length(args) >= 4) {
  ncores <- as.integer(args[4])
}
generate_sitepscore(dir_path = dir_path, rep_id = rep_id, tipidx_toprune = tipidx_toprune, ncores = ncores, overwrite = FALSE)

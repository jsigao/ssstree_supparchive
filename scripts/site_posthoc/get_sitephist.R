# this script contains functions that generate the parsimony mutation history of each site in the alignment
library(phangorn)
source("scripts/beasttree_readwrite/tree_proc_functions.R")

nucstates <- c("a", "c", "g", "t")
nucstates_all <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y", "k", "v", "h", "d", "b", "n", "?", "-")
contrast <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1), c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1), 
                     c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)), 
                   18, 4, dimnames = list(nucstates_all, nucstates))

get_anc <- function(tree, ancprob, contrast) {
  
  nstates <- ncol(contrast)
  edge <- tree$edge
  nedges <- nrow(edge)
  nnodes <- nedges + 1L
  root <- unique(edge[!edge[, 1] %in% edge[, 2], 1])
  interidx <- which(edge[, 2] %in% edge[, 1])
  states_idx <- match(colnames(contrast), rownames(contrast))
  
  node_states <- integer(nnodes)
  for (k in 1:nnodes) {
    idx <- which(ancprob[k, ] == 1)
    if (length(idx) == 1) node_states[k] <- idx 
  }
  
  k <- 0L # up pass
  while (k < nedges) {
    k <- k + 2L
    tnode <- edge[k - 1, 1]
    dnodes <- edge[k - (1:0), 2]
    tstate <- node_states[tnode]
    dstates <- node_states[dnodes]
    if (all(c(tstate, dstates) %in% states_idx) || all(c(tstate, dstates) %in% states_idx == FALSE)) next
    
    if (tstate %in% states_idx) {
      editd1 <- dstates[1] %in% states_idx == FALSE && (k - 1) %in% interidx == FALSE
      editd2 <- dstates[2] %in% states_idx == FALSE && k %in% interidx == FALSE
      if (editd1) node_states[dnodes[1]] <- ifelse(ancprob[dnodes[1], tstate] > 0, tstate, sample.int(nstates, size = 1, prob = ancprob[dnodes[1], ]))
      if (editd2) node_states[dnodes[2]] <- ifelse(ancprob[dnodes[2], tstate] > 0, tstate, sample.int(nstates, size = 1, prob = ancprob[dnodes[2], ]))
    } else { # this internal node ambiguous
      if (tnode != root) {
        anode <- edge[which(tnode == edge[, 2]), 1]
        astate <- node_states[anode]
      }
      
      if (all(dstates %in% states_idx)) {
        if (length(unique(dstates)) == 1) {
          node_states[tnode] <- dstates[1]
        } else {
          if (tnode != root) {
            if (astate %in% dstates) {
              node_states[tnode] <- astate
            } else if (astate %in% states_idx) {
              tp <- ancprob[tnode, c(astate, dstates)]
              if (any(tp > 0)) {
                node_states[tnode] <- c(astate, dstates)[sample.int(3, size = 1, prob = tp)]
              } else {
                node_states[tnode] <- sample.int(nstates, size = 1, prob = ancprob[tnode, ])
              }
            }
          } else {
            tp <- ancprob[tnode, dstates]
            if (any(tp > 0)) {
              node_states[tnode] <- dstates[sample.int(2, size = 1, prob = tp)]
            } else {
              node_states[tnode] <- sample.int(nstates, size = 1, prob = ancprob[tnode, ])
            }
          }
        }
      } else {
        if (tnode != root) {
          if (astate %in% states_idx) {
            if (astate %in% dstates) {
              node_states[tnode] <- astate
            } else {
              if (dstates[1] %in% states_idx) {
                tp <- ancprob[tnode, ] * ancprob[dnodes[2], ]
                tp <- tp[c(astate, dstates[1])]
                if (any(tp > 0)) {
                  node_states[tnode] <- c(astate, dstates[1])[sample.int(2, size = 1, prob = tp)]
                }
              } else if (dstates[2] %in% states_idx) {
                tp <- ancprob[tnode, ] * ancprob[dnodes[1], ]
                tp <- tp[c(astate, dstates[2])]
                if (any(tp > 0)) {
                  node_states[tnode] <- c(astate, dstates[2])[sample.int(2, size = 1, prob = tp)]
                }
              }
            }
          }
        } else {
          if (dstates[1] %in% states_idx) {
            node_states[tnode] <- ifelse(ancprob[tnode, dstates[1]] > 0, dstates[1], sample.int(nstates, size = 1, prob = ancprob[tnode, ]))
          } else if (dstates[2] %in% states_idx) {
            node_states[tnode] <- ifelse(ancprob[tnode, dstates[2]] > 0, dstates[2], sample.int(nstates, size = 1, prob = ancprob[tnode, ]))
          }
        }
      }
    }
  }
  
  if (!node_states[root] %in% states_idx) node_states[root] <- sample.int(nstates, size = 1, prob = ancprob[root, ])
  
  edge <- edge[nedges:1, ]
  interidx <- which(edge[, 2] %in% edge[, 1])
  root <- unique(edge[!edge[, 1] %in% edge[, 2], 1])
  k <- 0L # down pass
  for (k in 1:nedges) {
    tnode <- edge[k, 2]
    tstate <- node_states[tnode]
    if (tstate %in% states_idx) next
    
    anode <- edge[k, 1]
    astate <- node_states[anode]
    node_states[tnode] <- ifelse(ancprob[tnode, astate] > 0, astate, sample.int(nstates, size = 1, prob = ancprob[tnode, ]))
  }
  
  return(node_states)
}


process_amb <- function(tree, ancstates, contrast) {
  nstates <- ncol(contrast)
  
  edge <- tree$edge
  nedges <- nrow(edge)
  root <- unique(edge[!edge[, 1] %in% edge[, 2], 1])
  interidx <- which(edge[, 2] %in% edge[, 1])
  eidx <- integer()
  
  for (k in 1:nedges) {
    statepair <- ancstates[edge[k, ]]
    if (length(unique(statepair)) == 1) next
    if (all(statepair %in% states_idx)) next
    if ((statepair[1] %in% states_idx) && (statepair[2] %in% states_idx == FALSE)) {
      if (contrast[statepair[2], statepair[1]] == 1) {
        ancstates[edge[k, 2]] <- ancstates[edge[k, 1]]
        if (k %in% interidx) eidx <- c(eidx, k)
      }
    } else if ((statepair[1] %in% states_idx == FALSE) && (statepair[2] %in% states_idx == FALSE)) {
      idx_start <- contrast[statepair[1], ]
      idx_end <- contrast[statepair[2], ]
      idx_intersect <- as.integer(idx_start * idx_end)
      if (any(idx_intersect > 0)) {
        id <- as.integer(which(apply(contrast, 1, function(x) identical(as.integer(x), idx_intersect))))
        if (length(id) > 1) stop("unexpected")
        ancstates[edge[k, ]] <- id
        if (k %in% interidx) eidx <- c(eidx, k)
      }
    } else if ((statepair[1] %in% states_idx == FALSE) && (statepair[2] %in% states_idx)) {
      if (contrast[statepair[1], statepair[2]] == 1) {
        ancstates[edge[k, 1]] <- ancstates[edge[k, 2]]
        if (edge[k, 1] != root) eidx <- c(eidx, which(edge[, 2] == edge[k, 1]))
      }
    }
  }
  
  if (length(eidx) > 0) {
    eidx <- eidx[order(eidx, decreasing = TRUE)]
    for (k in eidx) {
      node <- edge[k, 2]
      astate <- ancstates[edge[k, 1]]
      dstates <- ancstates[edge[edge[, 1] == node, 2]]
      if (dstates[1] != dstates[2] && astate %in% dstates) ancstates[node] <- astate
    }
  }
  
  edge <- edge[nedges:1, ]
  interidx <- which(edge[, 2] %in% edge[, 1])
  for (k in 1:nedges) {
    node <- edge[k, 2]
    tstate <- ancstates[node]
    if (tstate %in% states_idx) next
    astate <- ancstates[edge[k, 1]]
    if (tstate == astate) next
    idx_start <- contrast[astate, ]
    
    if (!k %in% interidx) {
      idx_end <- contrast[tstate, ]
      idx_intersect <- as.integer(idx_start * idx_end)
      if (all(idx_intersect == 0) || any(idx_end - idx_start < 0)) {
        next
      }
      ancstates[node] <- astate
    } else {
      didx <- which(edge[, 1] == node)
      dstates <- ancstates[edge[didx, 2]]
      if (all(dstates %in% states_idx)) next
      
      tpstates <- unique(c(tstate, dstates))
      modify <- TRUE
      for (l in 1:length(tpstates)) {
        idx_end <- contrast[tpstates[l], ]
        idx_intersect <- as.integer(idx_start * idx_end)
        if (all(idx_intersect == 0) || any(idx_end - idx_start < 0)) {
          modify <- FALSE
          break
        }
      }
      
      if (modify) ancstates[node] <- astate
    }
  }
  
  return(ancstates)
}

calcwrite_sitephist <- function(trees_path, ncores = 1L, overwrite = FALSE) {
  
  dir_path <- dirname(trees_path)
  out_path <- paste0(dir_path, "/", gsub("\\.trees$", ".log", gsub(".*run", "sitephist_run", basename(trees_path))))
  ngen_burnin <- 0L
  if (grepl("burninNgen", trees_path)) ngen_burnin <- as.integer(gsub(".*burninNgen(\\d+).*", "\\1", trees_path))
  
  ntrees <- as.integer(system(paste0("grep -c 'tree STATE_' ", trees_path), intern = TRUE))
  ntrees_done <- 0L
  if (file.exists(out_path) && overwrite == FALSE) {
    nlines <- as.integer(system(paste0("wc -l < ", out_path), intern = TRUE))
    if (nlines == ntrees + 1L) return()
    
    ncol_tar <- 2L
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
  cat("processing ", out_path, "\n")
  
  datadir_path <- gsub("/analyses", "/data", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  ali_path <- list.files(datadir_path, full.names = TRUE, pattern = "^ali\\.fasta$")
  ali <- ape::read.dna(ali_path, format = "fasta", as.character = T)
  
  sp_path <- list.files(datadir_path, full.names = TRUE, pattern = "^site_patterns_notuseambi\\.tsv$")
  sp_df <- read.table(sp_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  spmat <- do.call(cbind, strsplit(sp_df$pattern, ""))
  rownames(spmat) <- rownames(ali)
  rm(ali)
  gc()
  
  nsp_total <- ncol(spmat)
  ntips <- nrow(spmat)
  tiplabs <- rownames(spmat)
  
  spdat <- phangorn::phyDat(spmat, compress = FALSE)
  nsp_this <- unique(lengths(spdat))
  if (length(nsp_this) != 1 || nsp_this != nsp_total) stop("wrong site pattern number")
  
  # read trees
  trees <- ape::read.nexus(trees_path)
  if (ntips != length(trees[[1]]$tip.label)) stop("tip numbers do not match")
  if (identical(tiplabs, trees[[1]]$tip.label) == FALSE) {
    trees <- read.nexus.beast(trees_path)
    if (identical(tiplabs, trees[[1]]$tip.label) == FALSE) stop("tip labels do not match")
  }
  ngen_trees <- as.integer(gsub("STATE_", "", names(trees)))
  if (ntrees != length(trees)) stop("tree number do not match")
  if (ntrees_done > 0) {
    trees <- trees[(ntrees_done + 1L):ntrees]
    ngen_trees <- ngen_trees[(ntrees_done + 1L):ntrees]
    ntrees <- length(trees)
    gc()
  }
  ngen_trees <- ngen_trees + ngen_burnin
  
  out_colnames <- c("state", "history")
  
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores, outfile = "")
    parallel::clusterEvalQ(cl, {library(phangorn)})
    # parallel::clusterExport(cl = cl, varlist = c("trees", "spmat", "spdat", "contrast", "ngen_trees", "process_amb", "get_anc"), envir = environment())
    parallel::clusterExport(cl = cl, varlist = c("trees", "spmat", "spdat", "contrast", "ngen_trees", "get_anc"), envir = environment())
    sid <- 1L
    ncores_multi <- 5L
    
    while (sid <= ntrees) {
      eid <- sid + ncores * ncores_multi - 1L
      if (eid > ntrees) eid <- ntrees
      
      mutstrs <- parallel::parSapply(cl, sid:eid, function(l) {
        tree <- trees[[l]]
        tree_postord <- ape::reorder.phylo(tree, order = "postorder")
        
        pscores <- phangorn::parsimony(tree = tree, data = spdat, method = "sankoff", site = "site")
        pscore_total <- sum(pscores)
        
        sidx <- which(pscores > 0)
        nsp_this <- length(sidx)
        spdat_this <- phangorn::phyDat(spmat[, sidx], compress = FALSE)
        # ancrec <- phangorn::ancestral.pars(tree = tree, data = spdat_this, type = "ACCTRAN", return = "phyDat")
        ancprobs <- phangorn::ancestral.pars(tree = tree, data = spdat_this, type = "ACCTRAN", return = "prob")
        
        eidx_all <- vector("list", nsp_this)
        node.states_all <- vector("list", nsp_this)
        for (k in 1:nsp_this) {
          # ancrec_this <- as.integer(sapply(ancrec, "[[", k))
          # if (any(!ancrec_this %in% states_idx)) {
          #   ancrec_this <- process_amb(tree_postord, ancrec_this, contrast)
          # }
          ancprob_this <- do.call(rbind, lapply(ancprobs, function(x) x[k, ]))
          ancrec_this <- get_anc(tree_postord, ancprob_this, contrast)
          node.states <- matrix(ancrec_this[tree$edge], ncol = 2)
          eidx_all[[k]] <- which(node.states[, 1] != node.states[, 2])
          node.states_all[[k]] <- node.states
        }
        
        nevents_total <- sum(lengths(eidx_all))
        if (nevents_total != pscore_total) message("generation ", ngen_trees[l], ": ", "parsimony score is ", pscore_total, " while reconstruction contains ", nevents_total, " events")
        
        cid <- 0L
        mutmat <- matrix(nrow = nevents_total, ncol = 4)
        colnames(mutmat) <- c("nodeid", "sitepatid", "from", "to")
        for (k in 1:nsp_this) {
          node.states <- node.states_all[[k]]
          eidx <- eidx_all[[k]]
          nevent <- length(eidx)
          mutmat[1:nevent + cid, ] <- cbind(tree$edge[eidx, 2], rep(sidx[k], nevent), node.states[eidx, , drop = FALSE])
          cid <- cid + nevent
        }
        
        mutdf <- data.frame(mutmat)
        mutdf_grbynode <- split(mutdf[, -1], mutdf[, 1])
        nnodes <- length(mutdf_grbynode)
        nodeidx <- names(mutdf_grbynode)
        mutstr <- paste(sapply(1:nnodes, function(x) paste0(nodeidx[x], "[", "{", paste(apply(mutdf_grbynode[[x]], 1, function(x) paste(x, collapse = ",")), collapse = "},{"), "}", "]")), collapse = ",")
        
        return(mutstr)
      })
      
      out_df <- data.frame(ngen_trees[sid:eid], mutstrs)
      colnames(out_df) <- out_colnames
      first_write <- ntrees_done + sid == 1
      write.table(out_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
      
      cat("processed gen ", ngen_trees[eid], "\n")
      sid <- eid + 1L
    }
    
    parallel::stopCluster(cl)
  } else {
    
    niter_out <- 100L
    mutstrs <- character(niter_out)
    ntrees_processed <- 1L
    write_now <- FALSE
    
    for (l in 1:ntrees) {
      tree <- trees[[l]]
      tree_postord <- ape::reorder.phylo(tree, order = "postorder")
      
      pscores <- phangorn::parsimony(tree = tree, data = spdat, method = "sankoff", site = "site")
      pscore_total <- sum(pscores)
      
      sidx <- which(pscores > 0)
      nsp_this <- length(sidx)
      spdat_this <- phangorn::phyDat(spmat[, sidx], compress = FALSE)
      ancprobs <- phangorn::ancestral.pars(tree = tree, data = spdat_this, type = "ACCTRAN", return = "prob")
      
      eidx_all <- vector("list", nsp_this)
      node.states_all <- vector("list", nsp_this)
      for (k in 1:nsp_this) {
        ancprob_this <- do.call(rbind, lapply(ancprobs, function(x) x[k, ]))
        ancrec_this <- get_anc(tree_postord, ancprob_this, contrast)
        node.states <- matrix(ancrec_this[tree$edge], ncol = 2)
        eidx_all[[k]] <- which(node.states[, 1] != node.states[, 2])
        node.states_all[[k]] <- node.states
      }
      
      nevents_total <- sum(lengths(eidx_all))
      if (nevents_total != pscore_total) message("generation ", ngen_trees[l], ": ", "parsimony score is ", pscore_total, " while reconstruction contains ", nevents_total, " events")
      
      cid <- 0L
      mutmat <- matrix(nrow = nevents_total, ncol = 4)
      colnames(mutmat) <- c("nodeid", "sitepatid", "from", "to")
      for (k in 1:nsp_this) {
        node.states <- node.states_all[[k]]
        eidx <- eidx_all[[k]]
        nevent <- length(eidx)
        mutmat[1:nevent + cid, ] <- cbind(tree$edge[eidx, 2], rep(sidx[k], nevent), node.states[eidx, , drop = FALSE])
        cid <- cid + nevent
      }
      
      mutdf <- data.frame(mutmat)
      mutdf_grbynode <- split(mutdf[, -1], mutdf[, 1])
      nnodes <- length(mutdf_grbynode)
      nodeidx <- names(mutdf_grbynode)
      mutstr <- paste(sapply(1:nnodes, function(x) paste0(nodeidx[x], "[", "{", paste(apply(mutdf_grbynode[[x]], 1, function(x) paste(x, collapse = ",")), collapse = "},{"), "}", "]")), collapse = ",")
      
      # add state id then write out
      if (l %% niter_out == 0) {
        mutstrs[niter_out] <- mutstr
        write_now <- TRUE
      } else {
        mutstrs[l %% niter_out] <- mutstr
        if (l == ntrees) {
          mutstrs <- mutstrs[1:(l %% niter_out)]
          write_now <- TRUE
        }
      }
      
      if (write_now) {
        out_df <- data.frame(ngen_trees[ntrees_processed:l], mutstrs)
        colnames(out_df) <- out_colnames
        
        first_write <- (ntrees_done + l == niter_out) || (ntrees_done + l < niter_out && l == ntrees)
        write.table(out_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
        
        ntrees_processed <- l + 1
        cat("processed gen ", ngen_trees[l], "\n")
        mutstrs <- character(niter_out)
        write_now <- FALSE
      }
    }
  }
  
}


generate_sitephist <- function(dir_path, rep_id, ncores = 1L, overwrite = FALSE) {
  ridstr <- formatC(rep_id, width = 2, format = "d", flag = 0)
  # trees_path <- list.files(dir_path, full.name = T, pattern = paste0("^run", ridstr, "\\.trees$"))
  # trees_path <- list.files(dir_path, full.name = T, pattern = paste0("^topo_run", ridstr, "_resampfreq100000_burninNgen200000000\\.trees$"))
  trees_path <- list.files(dir_path, full.name = T, pattern = paste0("^run", ridstr, "_resampfreq100000_burninNgen200000000\\.trees$"))
  if (length(trees_path) != 1) stop("cannot find tree file")
  
  calcwrite_sitephist(trees_path, ncores = ncores, overwrite = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
rep_id <- as.integer(args[2])
ncores <- 1L
if (length(args) >= 3) {
  ncores <- as.integer(args[3])
}
generate_sitephist(dir_path = dir_path, rep_id = rep_id, ncores = ncores)

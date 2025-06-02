# this script contains functions that compute the parsimony number of geographic dispersal events among all areas and between each pair of areas
library(phangorn)
source("scripts/beasttree_readwrite/tree_proc_functions.R")

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

get_pairwisemat <- function(node.states, contrast) {
  # count the number of events between each pair
  nstates <- ncol(contrast)
  nchanges_pairwise <- matrix(0, nrow = nstates, ncol = nstates, byrow = T)
  
  node.states_diff <- node.states[node.states[, 1] != node.states[, 2], , drop = F]
  if (nrow(node.states_diff) == 0) return(nchanges_pairwise)
  
  statepairs_tab <- table(apply(node.states_diff, 1, function(x) paste(x, collapse = ", ")))
  statepairs_nchanges <- as.integer(statepairs_tab)
  statepairs_mat <- matrix(unlist(strsplit(names(statepairs_tab), ", ")), ncol = 2, byrow = T)
  mode(statepairs_mat) <- "integer"
  nchanges_pairwise[statepairs_mat] <- statepairs_nchanges
  return(nchanges_pairwise)
}

get_mutstr <- function(node.states, tree) {
  eidx <- which(node.states[, 1] != node.states[, 2])
  nevents <- length(eidx)
  node.states_diff <- node.states[eidx, , drop = FALSE]
  
  mutmat <- cbind(tree$edge[eidx, 2], node.states_diff)
  colnames(mutmat) <- c("nodeid", "from", "to")
  mutmat <- mutmat[order(mutmat[, 1]), ]
  
  mutstr <- paste0("{", paste(apply(mutmat, 1, function(x) paste(x, collapse = ",")), collapse = "},{"), "}")
  return(mutstr)
}

get_pair_name <- function(mat) {
  from <- c(col(mat)[lower.tri(mat)], row(mat)[lower.tri(mat)])
  to <- c(row(mat)[lower.tri(mat)], col(mat)[lower.tri(mat)])
  pairs <- apply(cbind(from, to), 1, function(x) paste(x, collapse = ","))
  return(pairs)
}

remove_lastnlines <- function(out_path, n = 1L) {
  out_path_tmp <- paste0(out_path, "tmp")
  system(paste0("head -n -", n, " ", out_path, " > ", out_path_tmp, " && mv ", out_path_tmp, " ", out_path))
}

get_ntrees_done <- function(out_path, ncol_tar) {
  nlines <- as.integer(system(paste0("wc -l < ", out_path), intern = TRUE))
  
  row_last <- system(paste0("tail -n1 ", out_path), intern = TRUE)
  sep_char <- "\t"
  row_last_nseps <- nchar(row_last) - nchar(gsub(sep_char, "", row_last))
  row_last_lastchar <- system(paste0("tail -c 1 ", out_path), intern = TRUE)
  if (row_last_nseps != ncol_tar - 1 || (length(row_last_lastchar) == 1 && row_last_lastchar != "")) {
    remove_lastnlines(out_path)
    nlines <- nlines - 1L
  }
  
  return(nlines)
}

get_pscore_pairwise <- function(file_path, trees_path, tnevent_outpath, pnevent_outpath, hist_outpath, pars_type = "ACCTRAN", ncores = 1L, overwrite = FALSE) {
  
  if (file.exists(file_path) == FALSE || file.exists(trees_path) == FALSE) return()
  
  # read in geographic data
  geo_df <- read.table(file_path, header = T, sep = ",", check.names = F, stringsAsFactors = F)
  states <- unique(geo_df$geography)
  states <- sort(unique(unlist(strsplit(states, "\\|"))))
  states <- states[states != "?"]
  
  states_amb <- sort(unique(geo_df$geography))
  states_amb <- c(states, states_amb[!states_amb %in% states])
  nstates <- length(states)
  nstates_amb <- length(states_amb)
  contrast <- matrix(0, ncol = nstates, nrow = nstates_amb)
  rownames(contrast) <- states_amb
  colnames(contrast) <- states
  diag(contrast) <- 1L
  if (nstates_amb > nstates) {
    for (j in (nstates + 1L):nstates_amb) {
      if (states_amb[j] == "?") {
        contrast[j, ] <- 1L
      } else if (grepl("\\|", states_amb[j])) {
        states_this <- sort(unlist(strsplit(states_amb[j], "\\|")))
        contrast[j, match(states_this, states)] <- 1L
      } else {
        stop("not supported yet")
      }
    }
  }
  
  tipstates <- geo_df$geography
  names(tipstates) <- geo_df$taxon_id
  data_this <- phangorn::phyDat(t(t(tipstates)), type = "USER", contrast = contrast)
  
  # read in trees
  ngen_burnin <- 0L
  if (grepl("burninNgen", trees_path)) ngen_burnin <- as.integer(gsub(".*burninNgen(\\d+).*", "\\1", trees_path))
  
  ntrees <- as.integer(system(paste0("grep -c 'tree STATE_' ", trees_path), intern = TRUE))
  ntrees_done <- 0L
  if ((file.exists(pnevent_outpath) || file.exists(hist_outpath)) && overwrite == FALSE) {
    nlines1 <- 0L
    if (file.exists(pnevent_outpath)) nlines1 <- get_ntrees_done(pnevent_outpath, ncol_tar = nstates * (nstates - 1) + 1L)
    nlines2 <- 0L
    if (file.exists(hist_outpath)) nlines2 <- get_ntrees_done(hist_outpath, ncol_tar = 2L)
    
    if (nlines1 == ntrees + 1L && nlines2 == ntrees + 1L) return()
    if (nlines1 != nlines2) {
      if (nlines1 > nlines2) {
        remove_lastnlines(pnevent_outpath, nlines1 - nlines2)
      } else {
        remove_lastnlines(hist_outpath, nlines2 - nlines1)
      }
    }
    
    ntrees_done <- max(0L, min(nlines1, nlines2) - 1L)
  }
  cat("processing ", pnevent_outpath, "\n")
  
  trees <- ape::read.nexus(trees_path)
  if (length(tipstates) != length(trees[[1]]$tip.label)) stop("tip numbers do not match")
  if (identical(sort(names(tipstates)), sort(trees[[1]]$tip.label)) == FALSE) {
    trees <- read.nexus.beast(trees_path)
    if (identical(sort(names(tipstates)), sort(trees[[1]]$tip.label)) == FALSE) stop("tip labels do not match")
  }
  if (ntrees != length(trees)) stop("tree number do not match")
  
  ngen_trees <- as.integer(gsub("STATE_", "", names(trees)))
  ngen_trees <- ngen_trees + ngen_burnin
  
  pscores <- as.integer(vapply(trees, function(tree) phangorn::parsimony(tree = tree, data = data_this, method = "sankoff"), FUN.VALUE = double(1L), USE.NAMES = FALSE))
  if (file.exists(tnevent_outpath) == FALSE || overwrite == TRUE) {
    out_df <- data.frame(ngen_trees, pscores)
    colnames(out_df) <- c("state", "total")
    write.table(out_df, file = tnevent_outpath, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  rm(pscores)
  gc()
  
  if (ntrees_done > 0) {
    trees <- trees[(ntrees_done + 1L):ntrees]
    ngen_trees <- ngen_trees[(ntrees_done + 1L):ntrees]
    ntrees <- length(trees)
    gc()
  }
  
  amat <- matrix(ncol = nstates, nrow = nstates)
  out_colnames <- c("state", get_pair_name(amat))
  
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores, outfile = "")
    parallel::clusterEvalQ(cl, {library(phangorn)})
    parallel::clusterExport(cl = cl, varlist = c("trees", "data_this", "contrast", "get_anc", "pars_type"), envir = environment())
    sid <- 1L
    ncores_multi <- 100L
    
    while (sid <= ntrees) {
      eid <- sid + ncores * ncores_multi - 1L
      if (eid > ntrees) eid <- ntrees
      treeidx <- sid:eid
      ntrees_this <- length(treeidx)
      
      ps <- as.integer(parallel::parSapply(cl, treeidx, function(l) {
        phangorn::parsimony(tree = trees[[l]], data = data_this, method = "sankoff")
      }, USE.NAMES = FALSE))
      
      ancrecs <- parallel::parLapply(cl, treeidx, function(l) {
        tree <- trees[[l]]
        tree_postord <- ape::reorder.phylo(tree, order = "postorder")
        
        ancprob <- phangorn::ancestral.pars(tree = tree, data = data_this, type = pars_type, return = "prob")
        ancprob_plain <- do.call(rbind, lapply(ancprob, function(x) x[1, ]))
        ancrec <- get_anc(tree_postord, ancprob_plain, contrast)
        return(ancrec)
      })
      
      mutstrs <- character(ntrees_this)
      pevent_mat <- matrix(ncol = nstates * (nstates - 1L), nrow = ntrees_this)
      for (l in 1:ntrees_this) {
        ancrec <- ancrecs[[l]]
        tree <- trees[[treeidx[l]]]
        node.states <- matrix(ancrec[tree$edge], ncol = 2L)
        mutstrs[l] <- get_mutstr(node.states, tree)
        psmat <- get_pairwisemat(node.states, contrast)
        pevent_mat[l, ] <- c(t(psmat)[lower.tri(t(psmat))], psmat[lower.tri(psmat)])
        
        nevents <- sum(psmat)
        if (nevents != ps[l]) message("generation ", ngen_trees[treeidx[l]], ": ", "parsimony score is ", ps[l], " while reconstruction contains ", nevents, " events")
      }
      
      first_write <- ntrees_done + sid == 1
      
      out_df <- data.frame(ngen_trees[treeidx], pevent_mat)
      colnames(out_df) <- out_colnames
      write.table(out_df, file = pnevent_outpath, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
      
      out_df <- data.frame(ngen_trees[treeidx], mutstrs)
      colnames(out_df) <- c("state", "history")
      write.table(out_df, file = hist_outpath, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
      
      cat("processed gen ", ngen_trees[eid], "\n")
      sid <- eid + 1L
    }
    
    parallel::stopCluster(cl)
  } else {
    
    ps <- as.integer(vapply(trees, function(tree) phangorn::parsimony(tree = tree, data = data_this, method = "sankoff"), FUN.VALUE = double(1L), USE.NAMES = FALSE))
    
    niter_out <- 100L
    mutstrs <- character(niter_out)
    pevent_mat <- matrix(ncol = nstates * (nstates - 1L), nrow = niter_out)
    mode(pevent_mat) <- "integer"
    ntrees_processed <- 1L
    write_now <- FALSE
    
    for (l in 1:ntrees) {
      tree <- trees[[l]]
      tree_postord <- ape::reorder.phylo(tree, order = "postorder")
      
      ancprob <- phangorn::ancestral.pars(tree = tree, data = data_this, type = pars_type, return = "prob")
      ancprob_plain <- do.call(rbind, lapply(ancprob, function(x) x[1, ]))
      ancrec <- get_anc(tree_postord, ancprob_plain, contrast)
      node.states <- matrix(ancrec[tree$edge], ncol = 2L)
      
      mutstr <- get_mutstr(node.states, tree)
      psmat <- get_pairwisemat(node.states, contrast)
      pevent_this <- c(t(psmat)[lower.tri(t(psmat))], psmat[lower.tri(psmat)])
      
      nevents <- sum(psmat)
      if (nevents != ps[l]) message("generation ", ngen_trees[l], ": ", "parsimony score is ", ps[l], " while reconstruction contains ", nevents, " events")
      
      # add state id then write out
      if (l %% niter_out == 0) {
        mutstrs[niter_out] <- mutstr
        pevent_mat[niter_out, ] <- pevent_this
        write_now <- TRUE
      } else {
        mutstrs[l %% niter_out] <- mutstr
        pevent_mat[l %% niter_out, ] <- pevent_this
        if (l == ntrees) {
          mutstrs <- mutstrs[1:(l %% niter_out)]
          pevent_mat <- pevent_mat[1:(l %% niter_out), ]
          write_now <- TRUE
        }
      }
      
      if (write_now) {
        first_write <- (ntrees_done + l == niter_out) || (ntrees_done + l < niter_out && l == ntrees)
        
        out_df <- data.frame(ngen_trees[ntrees_processed:l], pevent_mat)
        colnames(out_df) <- out_colnames
        write.table(out_df, file = pnevent_outpath, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
        
        out_df <- data.frame(ngen_trees[ntrees_processed:l], mutstrs)
        colnames(out_df) <- c("state", "history")
        write.table(out_df, file = hist_outpath, quote = FALSE, sep = "\t", row.names = FALSE, append = !first_write, col.names = first_write)
        
        ntrees_processed <- l + 1
        cat("processed gen ", ngen_trees[l], "\n")
        mutstrs <- character(niter_out)
        pevent_mat <- matrix(ncol = nstates * (nstates - 1L), nrow = niter_out)
        mode(pevent_mat) <- "integer"
        write_now <- FALSE
      }
    }
    
  }
}

args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]
trees_path <- args[2]
tnevent_outpath <- args[3]
pnevent_outpath <- args[4]
hist_outpath <- args[5]
pars_type <- "ACCTRAN"
if (length(args) >= 6 && toupper(args[6]) %in% c("ACCTRAN", "MPR")) {
  pars_type <- toupper(args[6])
}
ncores <- 1L
if (length(args) >= 7) {
  ncores <- as.integer(args[7])
}

get_pscore_pairwise(file_path = file_path, trees_path = trees_path, tnevent_outpath = tnevent_outpath, pnevent_outpath = pnevent_outpath, hist_outpath = hist_outpath, pars_type = pars_type, ncores = ncores)

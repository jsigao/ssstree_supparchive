# this script contains functions that generate posterior trees whose branch lengths are in unit of expected number substitutions based on the input time trees and inferred branch rates
# this posterior phylograms will then be used as input trees for computing the site likelihood
library(ape)
source('scripts/beasttree_readwrite/tree_proc_functions.R')

calcMu <- function(mean, stdev) {
  return(log(mean / ((1 + (stdev * stdev) / (mean * mean))^0.5)))
}
calcSigma <- function(mean, stdev) {
  return((log(1 + (stdev * stdev) / (mean * mean)))^0.5)
}

get_branchlength_tree <- function(trees_path, clock_path) {
  
  ntrees <- as.integer(system(paste0("grep -c 'tree STATE_' ", trees_path), intern = TRUE))
  sys_name <- Sys.info()["sysname"]
  if (grepl("Darwin", sys_name)) {
    comd_this <- paste0("grep -Eo 'tree STATE_\\d+' ", trees_path)
  } else if (grepl("Linux", sys_name)) {
    comd_this <- paste0("grep -o 'tree STATE_[0-9]\\+' ", trees_path)
  } else {
    message("Unsupported operating system.\n")
    return()
  }
  ngen_trees <- as.integer(gsub("tree STATE_", "", system(comd_this, intern = TRUE)))
  
  clock_row1 <- read.table(clock_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F, nrows = 1)
  col_classes <- rep("NULL", ncol(clock_row1))
  col_classes[1] <- "integer"
  if (all(c("ucld.mean", "ucld.stdev") %in% colnames(clock_row1))) {
    clock_type <- "relaxed"
    col_idx <- which(colnames(clock_row1) %in% c("ucld.mean", "ucld.stdev"))
  } else if (any(colnames(clock_row1) == "clock.rate")) {
    clock_type <- "strict"
    col_idx <- which(colnames(clock_row1) == "clock.rate")
  } else {
    col_idx <- grep("\\.rate$", colnames(clock_row1))
    if (length(col_idx) > 1) {
      clock_type <- "local"
    } else {
      stop("cannot figure out clock type")
    }
  }
  col_classes[col_idx] <- "numeric"
  clock_df <- read.table(clock_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F, colClasses = col_classes)
  clock_df <- clock_df[clock_df[, 1] %in% ngen_trees, ]
  
  if (clock_type == "relaxed") {
    trees <- readBeast(trees_path)
  } else {
    trees <- read.nexus.beast(trees_path)
  }
  cat("trees read\n")
  if (length(trees) != ntrees) stop("numbers of trees fetched from two sources do not match")
  
  bl_mins <- sapply(trees, function(x) min(x$edge.length))
  idx <- which(bl_mins == 0)
  if (length(idx) > 0) {
    bl_minpos <- sapply(trees, function(x) min(x$edge.length[x$edge.length > 0]))
    bl_tol <- max(min(bl_minpos) * 0.5, 1e-8)
    for (l in idx) {
      trees[[l]]$edge.length[trees[[l]]$edge.length == 0] <- bl_tol
    }
  }
  
  if (clock_type == "relaxed") {
    ntips <- length(trees[[1]]$tip.label)
    nedges <- ntips * 2 - 2L
    probs <- (1:nedges - 0.5) * (1 / nedges)
    
    for (j in 1:ntrees) {
      tree <- trees[[j]]
      rcats <- as.integer(gsub("rateCat=", "", tree$node.data)) + 1L
      m <- clock_df[j, 2]
      s <- clock_df[j, 3]
      rs <- qlnorm(probs[rcats], meanlog = calcMu(m, s), sdlog = calcSigma(m, s))
      tree$edge.length <- tree$edge.length * rs[tree$edge[, 2]]
      trees[[j]] <- tree
    }
  } else if (clock_type == "strict") {
    for (j in 1:ntrees) {
      tree <- trees[[j]]
      tree$edge.length <- tree$edge.length * clock_df[j, 2]
      trees[[j]] <- tree
    }
  } else {
    
    olcc <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    
    dir_path <- dirname(trees_path)
    xml_path <- list.files(dir_path, full.names = TRUE, pattern = "\\.xml$")
    if (length(xml_path) == 0) stop("cannot find xml file")
    xml_path <- xml_path[1]
    xml_text <- scan(xml_path, what = character(), sep = "\n", blank.lines.skip = F)
    
    taxon_lidx <- grep("<taxon id=\"", xml_text)
    taxons <- gsub(".*<taxon id=\"|\">", "", xml_text[taxon_lidx])
    
    taxa_lidx <- grep("<taxa id=\"", xml_text)
    if (length(taxa_lidx) <= 1) stop("cannot find monophyly constraint")
    taxa_lidx <- taxa_lidx[-1]
    constraint_names <- gsub(".*<taxa id=\"|\">", "", xml_text[taxa_lidx])
    nconstraints <- length(taxa_lidx)
    taxa_elidx <- grep("</taxa>", xml_text)[-1]
    taxa_elidx <- taxa_elidx[1:nconstraints]
    
    taxons_all <- vector("list", nconstraints)
    names(taxons_all) <- constraint_names
    for (k in 1:nconstraints) {
      slid <- taxa_lidx[k] + 1
      elid <- taxa_elidx[k] - 1
      taxon_lidx_this <- grep("<taxon idref=\"", xml_text[slid:elid]) + slid - 1L
      taxons_this <- gsub(".*<taxon idref=\"|\"/>", "", xml_text[taxon_lidx_this])
      taxons_all[[k]] <- sort(taxons_this)
    }
    
    rate_lidx <- grep(paste0("<parameter id.*=\"", ".*\\.rate\""), xml_text)
    localclock_slid <- grep("localClockModel id=\"", xml_text)
    localclock_elid <- grep("</localClockModel>", xml_text)
    rate_lidx <- rate_lidx[rate_lidx > localclock_slid & rate_lidx < localclock_elid]
    
    claderate_lidx <- rate_lidx + 1L
    ratewclade_lidx <- rate_lidx[grep("<taxa idref=\"", xml_text[claderate_lidx])]
    claderate_lidx <- claderate_lidx[grep("<taxa idref=\"", xml_text[claderate_lidx])]
    
    incstem_lidx <- rate_lidx - 1L
    incstem_lidx <- incstem_lidx[grep("<clade includeStem=\"", xml_text[incstem_lidx])]
    
    clade_names <- gsub(".*<taxa idref=\"|\"/>", "", xml_text[claderate_lidx])
    clade_incstems <- as.logical(gsub(".*<clade includeStem=\"|\">", "", xml_text[incstem_lidx]))
    ratewclade_names <- gsub(".*<parameter id.*=\"(.*)\\.rate\".*", "\\1.rate", xml_text[ratewclade_lidx])
    
    taxons_all_wrate <- taxons_all[clade_names]
    taxons_sorted <- sort(taxons)
    taxons_all_wrate <- lapply(taxons_all_wrate, function(x) match(x, taxons_sorted))
    
    claderate_df <- cbind(clade_names, ratewclade_names, clade_incstems, taxons_all_wrate)
    ncrs <- nrow(claderate_df)
    
    ratebg_lid <- rate_lidx[!rate_lidx %in% ratewclade_lidx]
    ratebg_name <- gsub(".*<parameter id.*=\"(.*)\\.rate\".*", "\\1.rate", xml_text[ratebg_lid])
    
    crs <- c(ratewclade_names, ratebg_name)
    rate_names <- colnames(clock_df)[-1]
    if (identical(sort(unique(crs)), sort(rate_names)) == FALSE) stop("rate names do not match")
    
    ntips <- length(trees[[1]]$tip.label)
    for (j in 1:ntrees) {
      tree <- trees[[j]]
      clades <- ape::prop.part(tree, check.labels = FALSE)
      clades <- c(clades[-1], as.list(1:ntips))
      nclades <- length(clades)
      ninters <- nclades - ntips
      
      cr_idx <- rep(ncrs + 1L, nclades)
      for (k in 1:nclades) {
        ntaxa_this <- length(clades[[k]])
        ntaxa_inter <- sapply(claderate_df[, 4], function(x) sum(clades[[k]] %in% x))
        tid <- which(ntaxa_inter > 0)
        if (length(tid) == 1) {
          ntaxa_clade <- length(claderate_df[[tid, 4]])
          if ((ntaxa_inter[tid] == ntaxa_clade && claderate_df[[tid, 3]]) || (ntaxa_inter[tid] == ntaxa_this && ntaxa_this < ntaxa_clade)) {
            cr_idx[k] <- tid
          }
        }
      }
      
      cr_names <- crs[cr_idx]
      rate_idx <- match(cr_names, rate_names)
      rates <- clock_df[j, rate_idx + 1]
      
      nidx <- c(sapply(clades[1:ninters], function(x) ape::getMRCA(tree, x)), 1:ntips)
      rs <- rates[match(tree$edge[, 2], nidx)]
      
      tree$edge.length <- tree$edge.length * rs
      trees[[j]] <- tree
    }
    
    Sys.setlocale("LC_COLLATE", olcc)
  }
  
  return(trees)
}

generate_branchlength_tree <- function(dir_path, rep_id, overwrite = FALSE) {
  ridstr <- formatC(rep_id, width = 2, format = "d", flag = 0)
  trees_path <- list.files(dir_path, full.name = T, pattern = paste0("^run", ridstr, "\\.trees$"))
  if (length(trees_path) != 1) stop("cannot find tree file")
  
  clock_path <- list.files(dir_path, full.name = T, pattern = paste0("^param_run", ridstr, "\\.log$"))
  if (length(clock_path) != 1) stop("cannot find param file")
  
  out_path <- paste0(dir_path, "/", gsub("^run", "bltree_run", basename(trees_path)))
  if (file.exists(out_path) == FALSE || overwrite) {
    cat("processing ", out_path, "\n")
    trees <- get_branchlength_tree(trees_path, clock_path)
    write.nexus.beast(trees, file = out_path)
    
    tl_path <- paste0(dir_path, "/", "treelength_expsub", "_run", ridstr, ".log")
    if (file.exists(tl_path) == FALSE || overwrite) {
      tls <- sapply(trees, function(tree) sum(tree$edge.length))
      ngen_trees <- as.integer(gsub("STATE_", "", names(trees)))
      tl_df <- data.frame(state = ngen_trees, TLexpsub = tls, stringsAsFactors = FALSE)
      write.table(tl_df, file = tl_path, quote = FALSE, sep = "\t", row.names = FALSE)
    }
  }
}

args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
rep_id <- as.integer(args[2])
generate_branchlength_tree(dir_path = dir_path, rep_id = rep_id, overwrite = TRUE)

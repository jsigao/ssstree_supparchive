# this script contains functions that calculate the pairwise RF tree distances

get_rf <- function(trees_path, out_path, ncores = 1L) {
  
  trees <- ape::read.nexus(trees_path)
  ntrees <- length(trees)
  ngen_idx <- gsub("STATE_", "", names(trees))
  line0 <- paste(c("state", ngen_idx), collapse = ",")
  line1 <- paste0(ngen_idx[1], paste(rep(",", ntrees), collapse = ""))
  
  ntrees_done <- 0L
  if (file.exists(out_path)) {
    nlines <- as.integer(system(paste0("wc -l < ", out_path), intern = TRUE))
    if (nlines == ntrees + 1L) return()
    
    row_last <- system(paste0("tail -n1 ", out_path), intern = T)
    row_last_nseps <- nchar(row_last) - nchar(gsub(",", "", row_last))
    row_last_lastchar <- system(paste0("tail -c 1 ", out_path), intern = TRUE)
    if (row_last_nseps < ntrees || (length(row_last_lastchar) == 1 && row_last_lastchar != "")) {
      out_path_tmp <- paste0(out_path, "tmp")
      system(paste0("head -n -1 ", out_path, " > ", out_path_tmp, " && mv ", out_path_tmp, " ", out_path))
      nlines <- nlines - 1L
    }
    
    ntrees_done <- nlines - 1L
    file_conn <- file(out_path, open = "at")
    if (ntrees_done == 0) {
      writeLines(line1, file_conn)
      flush(file_conn)
      ntrees_done <- 1L
    }
  } else {
    file_conn <- file(out_path, open = "wt")
    writeLines(line0, file_conn)
    writeLines(line1, file_conn)
    flush(file_conn)
    ntrees_done <- 1L
  }
  
  ntips <- length(trees[[1]]$tip.label)
  rfMax <- 2 * ntips - 4L
  
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl = cl, varlist = c("trees"), envir = environment())
    
    partStr <- parallel::parLapply(cl, 1:ntrees, function(i) {
      ps <- ape::prop.part(trees[[i]], check.labels = FALSE)
      vapply(ps[-1], function(x) paste0(x, collapse = ","), character(1L), USE.NAMES = FALSE)
    })
    
    parallel::stopCluster(cl)
  } else {
    partStr <- lapply(trees, function(tree) {
      ps <- ape::prop.part(tree, check.labels = FALSE)
      vapply(ps[-1], function(x) paste0(x, collapse = ","), character(1L), USE.NAMES = FALSE)
    })
  }
  
  rm(trees)
  gc()
  gc()
  
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl = cl, varlist = c("partStr"), envir = environment())
    sid <- ntrees_done + 1L
    ncores_multi <- 50L
    
    while (sid <= ntrees) {
      eid <- sid + ncores * ncores_multi - 1L
      if (eid > ntrees) eid <- ntrees
      nthis <- eid - sid + 1L
      
      pps <- parallel::parLapply(cl, sid:eid, function(i) {
        part_this <- partStr[[i]]
        vapply(partStr[1:(i - 1)], function(x) length(intersect(part_this, x)), integer(1L), USE.NAMES = FALSE)
      })
      
      for (i in 1:nthis) {
        ps <- rfMax - 2 * pps[[i]]
        pstrs <- character(ntrees)
        pstrs[1:(sid + i - 2)] <- as.character(ps)
        line_this <- paste0(ngen_idx[(sid + i - 1)], ",", paste(pstrs, collapse = ","))
        writeLines(line_this, file_conn)
      }
      
      flush(file_conn)
      sid <- sid + nthis
    }
    
    parallel::stopCluster(cl)
  } else {
    for (i in (ntrees_done + 1L):ntrees) {
      ps <- vapply(partStr[1:(i - 1)], function(x) length(intersect(partStr[[i]], x)), integer(1L), USE.NAMES = FALSE)
      ps <- rfMax - 2 * ps
      pstrs <- character(ntrees)
      pstrs[1:(i - 1)] <- as.character(ps)
      line_this <- paste0(ngen_idx[i], ",", paste(pstrs, collapse = ","))
      writeLines(line_this, file_conn)
      flush(file_conn)
    }
  }
  
  close(file_conn)
}


args <- commandArgs(trailingOnly = TRUE)
trees_path <- args[1]
out_path <- args[2]
ncores <- 1L
if (length(args) > 2) {
  ncores <- as.integer(args[3])
}

get_rf(trees_path = trees_path, out_path = out_path, ncores = ncores)

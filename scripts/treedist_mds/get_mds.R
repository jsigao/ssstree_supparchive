# this script contains functions that generate the MDS coordinates of tree space
source('scripts/treedist_mds/mds_customized.R')

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

process_dir_dmat <- function(dir_path, dist_str = "RF_combined", mds_method = "cmds", ngen_freq = 50000L, burnin = 0.5, ndim = 2, overwrite = FALSE) {
  
  if (burnin <= 0) {
    burnin_str <- ""
    burnin_str2 <- ""
  } else if (burnin < 1) {
    burnin_str <- paste0("_", burnin * 100, "burninPerc")
    burnin_str2 <- paste0("_", "burninPerc", burnin * 100)
  } else {
    burnin_str <- paste0("_", as.integer(burnin), "burninNgen")
    burnin_str2 <- paste0("_", "burninNgen", as.integer(burnin))
  }
  
  file_path <- list.files(dir_path, full.name = TRUE, pattern = paste0("pairwise", dist_str, "_resampfreq", as.integer(ngen_freq), burnin_str2, "\\.log$"))
  if (length(file_path) == 0) return()
  if (length(file_path) > 1) stop("multiple files found")
  cat("processing ", dir_path, "\n", sep = "")
  
  outdir_path <- gsub("/analyses", "/intermediate", paste0(getwd(), "/", gsub("^\\.", "", dir_path)))
  if (dir.exists(outdir_path) == FALSE) dir.create(outdir_path, recursive = TRUE, showWarnings = FALSE)
  
  out_path <- paste0(outdir_path, "/", gsub("\\_resampfreq.*", "_resampfreq", basename(file_path)), as.integer(ngen_freq), burnin_str, "_", mds_method, ndim, "d.rds")
  if (file.exists(out_path) == TRUE && overwrite == FALSE) return()
  
  dist_mat <- get_dist_mat(file_path)
  res <- get_mds_wdmat(dist_mat, mds_method = mds_method, ndim = ndim)
  saveRDS(res, file = out_path)
}

args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
mds_method <- args[2]
ngen_freq <- as.integer(args[3])
dist_str <- args[4]
burnin <- 0
if (length(args) >= 5) {
  burnin <- as.numeric(args[5])
}
ndim <- 2
if (length(args) >= 6) {
  ndim <- as.integer(args[6])
}

process_dir_dmat(dir_path, dist_str = dist_str, mds_method = mds_method, ngen_freq = ngen_freq, burnin = burnin, ndim = ndim)

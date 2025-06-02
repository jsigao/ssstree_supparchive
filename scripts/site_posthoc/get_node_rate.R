# this script contains functions that obtain the numeric branch rates based on the sampled relaxed-clock rate categories and the hyper-parameter values
calcMu <- function(mean, stdev) {
  return(log(mean / ((1 + (stdev * stdev) / (mean * mean))^0.5)))
}
calcSigma <- function(mean, stdev) {
  return((log(1 + (stdev * stdev) / (mean * mean)))^0.5)
}

get_node_rate <- function(rep_path, clock_path, probs) {
  rep_df <- read.table(rep_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  nrow_log <- nrow(rep_df)
  
  clock_row1 <- read.table(clock_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F, nrows = 1)
  col_classes <- rep("NULL", ncol(clock_row1))
  col_classes[1] <- "integer"
  col_classes[colnames(clock_row1) %in% c("ucld.mean", "ucld.stdev")] <- "numeric"
  
  clock_df <- read.table(clock_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F, colClasses = col_classes)
  clock_df <- clock_df[clock_df[, 1] %in% rep_df[, 1], ]
  
  rep_mat <- as.matrix(rep_df[, -1])
  out_mat <- matrix(NA, ncol = ncol(rep_mat), nrow = nrow_log)
  for (j in 1:nrow_log) {
    m <- clock_df[j, 2]
    s <- clock_df[j, 3]
    out_mat[j, ] <- qlnorm(probs[rep_mat[j, ] + 1L], meanlog = calcMu(m, s), sdlog = calcSigma(m, s))
  }
  
  out_df <- data.frame(out_mat, stringsAsFactors = FALSE)
  out_df <- cbind(rep_df[, 1], out_df)
  colnames(out_df) <- colnames(rep_df)
  return(out_df)
}


generate_node_rate <- function(dir_path, overwrite = FALSE) {
  rep_paths <- grep("/cladeshared_rateCat_", list.files(dir_path, full.name = T, pattern = "run\\d+\\.log$"), value = T)
  if (length(rep_paths) == 0) return()
  nreps <- length(rep_paths)
  clock_paths <- grep("/param_", list.files(dir_path, full.name = T, pattern = "run\\d+\\.log$"), value = T)
  if (length(clock_paths) == 0 || length(clock_paths) != nreps) return()
  
  xml_path <- list.files(dir_path, full.names = TRUE, pattern = "\\.xml$")[1]
  xml_text <- scan(xml_path, what = character(), sep = "\n", blank.lines.skip = F)
  ntips <- length(grep("<taxon id=\"", xml_text))
  nedges <- ntips * 2 - 2L
  probs <- (1:nedges - 0.5) * (1 / nedges)
  
  out_paths <- gsub("_rateCat_", "_rate_", rep_paths)
  for (i in 1:nreps) {
    out_path <- out_paths[i]
    if (file.exists(out_path) == FALSE || overwrite) {
      cat("processing ", dir_path, " rep", i, "\n", sep = "")
      out_df <- get_node_rate(rep_paths[i], clock_paths[i], probs)
      write.table(out_df, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE)
    }
  }
}


args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]
generate_node_rate(dir_path = dir_path)

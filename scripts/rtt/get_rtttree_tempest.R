# this script contains function that performs RTT regression test based on the ML tree and sequence sampling times
library(lubridate)
source("scripts/beasttree_readwrite/tree_proc_functions.R")

rtttree_tempest <- function(dir_path) {
  tree_path <- list.files(dir_path, full.names = TRUE, pattern = "ml\\.tree$")
  if (length(tree_path) == 0) next
  tree <- read.nexus.beast(tree_path)
  outtree_path <- paste0(dir_path, "/", "ml_rootbyrtt.tree")
  
  date_path <- list.files(dir_path, full.names = TRUE, pattern = "seqdate\\.tsv$")
  date_df <- read.table(date_path, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  if (identical(sort(tree$tip.label), sort(date_df[, 1])) == FALSE) stop("tip names do not match")
  
  ntips <- length(tree$tip.label)
  datestrs <- date_df[match(tree$tip.label, date_df[, 1]), 2]
  dates <- lubridate::Date(ntips)
  nparts <- lengths(strsplit(datestrs, "-"))
  if (any(nparts %in% c(1:3) == FALSE)) stop("unexpected date format")
  
  if (any(nparts == 3)) dates[nparts == 3] <- as.Date(datestrs[nparts == 3], format = "%Y-%m-%d")
  if (any(nparts == 2)) dates[nparts == 2] <- as.Date(paste0(datestrs[nparts == 2], "-01"), format = "%Y-%m-%d")
  if (any(nparts == 1)) dates[nparts == 1] <- as.Date(paste0(datestrs[nparts == 1], "-01-01"), format = "%Y-%m-%d")
  if (any(is.na(dates))) stop("unexpected date format")
  
  datenums <- lubridate::decimal_date(dates)
  if (any(nparts == 2)) datenums[nparts == 2] <- datenums[nparts == 2] + (1/12) / 2
  if (any(nparts == 1)) datenums[nparts == 1] <- datenums[nparts == 1] + 1 / 2
  
  dateuncs <- numeric(ntips)
  if (any(nparts == 2)) dateuncs[nparts == 2] <- 1 / 12
  if (any(nparts == 1)) dateuncs[nparts == 1] <- 1
  
  rtttree <- ape::rtt(tree, tip.dates = datenums, objective = "rms")
  write.nexus.beast(rtttree, file = outtree_path, digits = 10)
  rtttree <- read.nexus.beast(outtree_path)
  
  rttdivs <- ape::node.depth.edgelength(rtttree)[1:ntips]
  rttlm <- lm(rttdivs ~ datenums)
  rttlmsum <- summary(rttlm)
  
  res <- list(rttdivs = rttdivs, tipdates = datenums, tipdateuncs = dateuncs, slope = rttlmsum$coefficients[2, 1], xinter = -rttlmsum$coefficients[1, 1] / rttlmsum$coefficients[2, 1], 
              rsquared = rttlmsum$r.squared, rms = sum(rttlmsum$residuals^2) / (ntips - 2), residuals = rttlmsum$residuals, 
              pvalue = unname(pf(rttlmsum$fstatistic[1], rttlmsum$fstatistic[2], rttlmsum$fstatistic[3], lower.tail = FALSE)))
  out_path <- paste0(dir_path, "/", "rttvals.rds")
  saveRDS(res, file = out_path)
}

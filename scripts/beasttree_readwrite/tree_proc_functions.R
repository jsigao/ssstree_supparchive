# this script contains functions that read or write BEAST trees
library(ape)

#' This function writes trees in a file with the BEAST nexus tree file format.
#' It is modified from ape's write.nexus function; the added argument `digits` is passed as an argument of ape's write.tree function.
#' The modification mainly allows the taxon names and other formats to follow BEAST's default
write.nexus.beast <- function (..., file = "", translate = TRUE, digits = 5) 
{
  obj <- ape:::.getTreesFromDotdotdot(...)
  ntree <- length(obj)
  cat("#NEXUS\n\n", file = file)
  
  N <- length(obj[[1]]$tip.label)
  cat("Begin taxa;\n", file = file, append = TRUE)
  cat(paste("\tDimensions ntax=", N, ";\n", sep = ""), file = file, append = TRUE)
  cat("\tTaxlabels\n", file = file, append = TRUE)
  tiplabel <- cleanTaxonName(obj[[1]]$tip.label)
  cat(paste("\t\t", tiplabel, sep = ""), sep = "\n", file = file, append = TRUE)
  cat("\t\t;\n", file = file, append = TRUE)
  cat("End;\n\n", file = file, append = TRUE)
  cat("Begin trees;\n", file = file, append = TRUE)
  
  if (translate) {
    cat("\tTranslate\n", file = file, append = TRUE)
    obj <- ape::.compressTipLabel(obj)
    
    tiplabel <- cleanTaxonName(attr(obj, "TipLabel"))
    X <- paste("\t\t", 1:N, " ", tiplabel, ",", sep = "")
    X[length(X)] <- gsub(",", "", X[length(X)])
    cat(X, file = file, append = TRUE, sep = "\n")
    cat("\t\t;\n", file = file, append = TRUE)
    class(obj) <- NULL
    for (i in 1:ntree) obj[[i]]$tip.label <- as.character(1:N)
  } else {
    if (is.null(attr(obj, "TipLabel"))) {
      for (i in 1:ntree) obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
    } else {
      attr(obj, "TipLabel") <- checkLabel(attr(obj, "TipLabel"))
      obj <- ape::.uncompressTipLabel(obj)
    }
  }
  
  title <- names(obj)
  if (is.null(title)) 
    title <- rep("UNTITLED", ntree)
  else {
    if (any(s <- title == "")) 
      title[s] <- "UNTITLED"
  }
  for (i in 1:ntree) {
    if (!inherits(obj[[i]], "phylo")) next
    root.tag <- ifelse(is.rooted(obj[[i]]), "= [&R] ", "= [&U] ")
    cat("tree", title[i], root.tag, file = file, append = TRUE)
    cat(write.tree(obj[[i]], file = "", digits = digits), "\n", sep = "", file = file, append = TRUE)
  }
  cat("End;\n", file = file, append = TRUE)
}

#' This is a helper function that format taxon names to deal with special characters in them
#' @param tiplabel taxon names
cleanTaxonName <- function(tiplabel) {
  special_characters_regex <- ".*[\\s\\.;,\"\'].*"
  spidx <- which(grepl(special_characters_regex, tiplabel, perl = TRUE))
  nsps <- length(spidx)
  if (nsps > 0) {
    for (i in 1:nsps) {
      if (grepl("\'", tiplabel[spidx[i]])) {
        if (grepl("\"", tiplabel[spidx[i]])) {
          stop("Illegal taxon name - contains both single and double quotes")
        }
        tiplabel[spidx[i]] <- paste0("\"", tiplabel[spidx[i]], "\"")
      } else {
        tiplabel[spidx[i]] <- paste0("\'", tiplabel[spidx[i]], "\'")
      }
    }
  }
  
  return(tiplabel)
}

#' This function read trees (the topology and branch lengths) from a BEAST nexus tree file.
#' It is modified from ape's read.nexus function to allow quotes in the taxon names to be correctly parsed
read.nexus.beast <- function (file, tree.names = NULL, force.multi = FALSE) 
{
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  LEFT <- grep("\\[", X)
  RIGHT <- grep("\\]", X)
  if (length(LEFT)) {
    w <- LEFT == RIGHT
    if (any(w)) {
      s <- LEFT[w]
      X[s] <- gsub("\\[[^]]*\\]", "", X[s])
    }
    w <- !w
    if (any(w)) {
      s <- LEFT[w]
      X[s] <- gsub("\\[.*", "", X[s])
      sb <- RIGHT[w]
      X[sb] <- gsub(".*\\]", "", X[sb])
      if (any(s < sb - 1)) X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
    }
  }
  endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
  semico <- grep(";", X)
  i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
  i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
  translation <- ifelse(length(i2) == 1 && i2 > i1, TRUE, FALSE) 
  if (translation) {
    end <- semico[semico > i2][1]
    x <- X[(i2 + 1):end]
    x <- gsub("^\\s+", "", x)
    x <- gsub("[,;]", "", x)
    x <- unlist(regmatches(x, regexpr("\\s+", x), invert = TRUE))
    x <- x[nzchar(x)]
    TRANS <- matrix(x, ncol = 2, byrow = TRUE)
    # TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2]) # this would incorrectly remove quote/prime inside a name
    TRANS[, 2] <- trimws(TRANS[, 2], which = "both", whitespace = "['\"]")
    n <- dim(TRANS)[1]
  }
  start <- ifelse(translation, semico[semico > i2][1] + 1, i1 + 1)
  end <- endblock[endblock > i1][1] - 1
  tree <- X[start:end]
  rm(X)
  tree <- tree[tree != ""]
  semico <- grep(";", tree)
  Ntree <- length(semico)
  if (Ntree == 1 && length(tree) > 1) 
    STRING <- paste(tree, collapse = "")
  else {
    if (any(diff(semico) != 1)) {
      STRING <- character(Ntree)
      s <- c(1, semico[-Ntree] + 1)
      j <- mapply(":", s, semico)
      if (is.list(j)) {
        for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], collapse = "")
      } else {
        for (i in 1:Ntree) STRING[i] <- paste(tree[j[, i]], collapse = "")
      }
    }
    else STRING <- tree
  }
  rm(tree)
  STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
  Ntree <- length(STRING)
  nms.trees <- sub(" *= *.*", "", STRING)
  nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", nms.trees, ignore.case = TRUE)
  STRING <- sub("^.*= *", "", STRING)
  STRING <- gsub(" ", "", STRING)
  colon <- grep(":", STRING)
  if (!length(colon)) {
    trees <- lapply(STRING, ape:::.cladoBuild)
  } else if (length(colon) == Ntree) {
    trees <- if (translation) 
      lapply(STRING, ape:::.treeBuildWithTokens)
    else lapply(STRING, ape:::.treeBuild)
  } else {
    trees <- vector("list", Ntree)
    trees[colon] <- lapply(STRING[colon], ape:::.treeBuild)
    nocolon <- (1:Ntree)[!1:Ntree %in% colon]
    trees[nocolon] <- lapply(STRING[nocolon], ape:::.cladoBuild)
    if (translation) {
      for (i in 1:Ntree) {
        tr <- trees[[i]]
        for (j in 1:n) {
          ind <- which(tr$tip.label[j] == TRANS[, 1])
          tr$tip.label[j] <- TRANS[ind, 2]
        }
        if (!is.null(tr$node.label)) {
          for (j in 1:length(tr$node.label)) {
            ind <- which(tr$node.label[j] == TRANS[, 1])
            tr$node.label[j] <- TRANS[ind, 2]
          }
        }
        trees[[i]] <- tr
      }
      translation <- FALSE
    }
  }
  for (i in 1:Ntree) {
    tr <- trees[[i]]
    if (!translation) n <- length(tr$tip.label)
  }
  if (Ntree == 1 && !force.multi) {
    trees <- trees[[1]]
    if (translation) {
      trees$tip.label <- if (length(colon)) 
        TRANS[, 2]
      else TRANS[, 2][as.numeric(trees$tip.label)]
    }
  } else {
    if (!is.null(tree.names)) 
      names(trees) <- tree.names
    if (translation) {
      if (length(colon) == Ntree) 
        attr(trees, "TipLabel") <- TRANS[, 2]
      else {
        for (i in 1:Ntree) trees[[i]]$tip.label <- TRANS[, 2][as.numeric(trees[[i]]$tip.label)]
        trees <- ape::.compressTipLabel(trees)
      }
    }
    class(trees) <- "multiPhylo"
    if (!all(nms.trees == "")) 
      names(trees) <- nms.trees
  }
  trees
}



#' Read in BEAST tree(s) log files and then parse each extended newick string (additional information other than the branch lengths) into a phylo object
#' @param file_path path to the BEAST tree(s) log file
#' @param tree_text text content of the BEAST tree(s) log file
#' @return a (list of) phylo object(s)
readBeast <- function(file_path = NULL, tree_text = NULL) {
  
  # read in the tree file content as vector of strings
  if (!is.null(file_path) && !is.null(tree_text)) {
    message("either file_path or tree_text needs to be provided, but not both; we will only use file_path then.\n")
    tree_text <- NULL
  }
  
  if (!is.null(file_path) && length(file_path) == 1 && file.exists(file_path)) {
    x_raw <- scan(file_path, quiet = T, sep = "\n", what = character(), quote = "'\"", blank.lines.skip = FALSE)
  } else if (!is.null(tree_text)) {
    x_raw <- tree_text
  } else {
    stop("don't know what to do.\n")
  }
  
  # only retain the relevant tree text
  start_linenum <- grep("^begin trees;$", tolower(x_raw))
  if (length(start_linenum) != 1) {
    stop("cannot find the tree starting line")
  }
  end_linenum <- grep("^end;$", tolower(x_raw))
  end_linenum <- min(end_linenum[end_linenum > start_linenum])
  if (length(end_linenum) != 1) {
    stop("cannot find the tree ending line")
  }
  x_raw <- x_raw[start_linenum:end_linenum]
  
  # new just take the tree lines and parse them
  treeline_num <- grep('\\[&R\\].*?;', x_raw)
  x <- x_raw[treeline_num]
  x <- gsub("^.*\\[&R\\] |^.*\\[&R\\]", "", x)
  my.data <- lapply(x, parseBeastNewick) # parsing
  tree_names <- gsub("tree (.*) = \\[&R\\].*", "\\1", x_raw[treeline_num])
  names(my.data) <- tree_names
  
  # decide whether tip label translation is needed
  tl_int <- suppressWarnings(as.integer(my.data[[1]]$tip.label))
  if (all(is.na(tl_int))) {
    translate <- FALSE
  } else if (all(!is.na(tl_int))) {
    translate <- TRUE
  } else {
    stop("the tip labels are mixed of pure numbers and characters")
  }
  
  # order all the tips (and of course node.data) the same way 1:tip.num
  if (translate) {
    ntips <- length(my.data[[1]]$tip.label)
    for (i in 1:length(my.data)) {
      tip.num <- as.integer(my.data[[i]]$tip.label)
      my.data[[i]]$edge[, 2][my.data[[i]]$edge[, 2] <= ntips] <- tip.num[my.data[[i]]$edge[, 2][my.data[[i]]$edge[, 2] <= ntips]]
      
      my.data[[i]]$node.data[tip.num] <- my.data[[i]]$node.data[1:ntips]
      my.data[[i]]$tip.label <- 1:ntips
    }
  }
  
  if (!translate) { # Trees have names in them
    for (i in 1:length(my.data)) {
      my.data[[i]]$tip.label <- trimws(my.data[[i]]$tip.label, which = "both", whitespace = "['\"]")
    }
    
    # Trees use numbers not names, we will use taxon name block of format "number name" to 
    # link the numbers used in the tree file to the names
  } else {
    
    x_raw <- x_raw[1:(treeline_num[1] - 1)]
    first <- grep('translate', tolower(x_raw)) + 1L
    last <- grep(';', x_raw)
    last <- min(last[last > first]) - 1L
    x <- trimws(x_raw[first:last], which = "both", whitespace = "[\\s,]")
    tip_ids <- as.integer(gsub("^(\\d+).*", "\\1", x))
    tip_names <- trimws(gsub("^\\d+ ", "", x), which = "both", whitespace = "['\"]")
    
    for (i in 1:length(my.data)) {
      my.data[[i]]$tip.label <- tip_names[match(my.data[[i]]$tip.label, tip_ids)]
    }
  }
  
  for (i in 1:length(my.data)) {
    class(my.data[[i]]) <- "phylo"
  }
  
  if (length(my.data) > 1) {
    return(my.data)
  } else {
    return(my.data[[1]])
  }
}


#' Read in an extended newick string, parse it to convert it into a phylo object
#' @param string the extended newick string
#' @return a phylo object
parseBeastNewick <- function(string) {
  
  # parse around front brackets and before back brackets while preserve them
  string <- gsub("\\(", "++(++", string)
  string <- gsub("\\)", "++)", string)
  string <- gsub("++++", "++", string, fixed = T)
  tree.vec <- strsplit(string, "++", fixed = T)[[1]][-1]
  
  n.tips <- length(grep("\\(", tree.vec)) + 1L # assuming at least two tips
  n.nodes <- as.integer(2 * n.tips - 1) # assuming fully bifurcating rooted tree
  n.edges <- as.integer(2 * n.tips - 2)
  edges <- matrix(0L, ncol = 2, nrow = n.edges, byrow = T)
  
  max.node <- n.tips
  at.node <- n.tips
  tip.node <- 0L # assuming fully bifurcating
  edges_rownum <- 1L
  
  tip.names <- character(n.tips)
  node_bls <- numeric(n.nodes)
  node_annotations <- character(n.nodes)
  
  # strip out extended newick part
  newicks <- gsub("\\[(.*?)\\]", "", tree.vec)
  newickexts <- gsub("]:[&", ",", tree.vec, fixed = T)
  extnewick_exists <- grepl("\\[&", tree.vec)
  
  is_firstvisits <- newicks == "("
  is_finalvisits <- grepl("\\)", newicks)
  has_subtendingbranchs <- grepl(":", newicks)
  
  for (i in 1:length(tree.vec)) {
    
    newick <- newicks[i]
    newickext <- newickexts[i]
    
    if (is_firstvisits[i]) { # adding a new node we've never seen before, guaranteed to be internal
      
      if (at.node != n.tips) {
        edges[edges_rownum, ] <- c(at.node, max.node + 1L)
        edges_rownum <- edges_rownum + 1L
      }
      max.node <- max.node + 1L
      at.node <- max.node
      
    } else if (is_finalvisits[i]) { # we're going back through a previously visited internal node
      
      old.node <- at.node
      
      if (has_subtendingbranchs[i]) { # with a subtending branch (i.e., non-root node)
        
        at.node <- edges[edges[, 2] == old.node, 1] # traverse back one node to their parent
        
        coloncommasplit <- strsplit(newick, ":|,|;")[[1]][-1]
        node_bls[old.node] <- as.numeric(coloncommasplit[1])
        
        if (length(coloncommasplit) == 3) { # there is a tip that's sister to this internal node, so this internal node is not root
          
          tip.node <- tip.node + 1L
          tip.names[tip.node] <- coloncommasplit[2]
          node_bls[tip.node] <- as.numeric(coloncommasplit[3])
          
          edges[edges_rownum, ] <- c(at.node, tip.node)
          edges_rownum <- edges_rownum + 1L
          
          if (extnewick_exists[i]) { # parse node annotation, the extended newick part, if there exists
            # when there are multiple extended newick parts
            newickext <- strsplit(newickext, ",(?![^[]*])", perl = T)[[1]]
            
            # first node in the string is the node we just passed through
            # second node in the string is the newly found tip
            for (j in 1:2) {
              if (grepl("\\[&", newickext[j])) {
                node_annotations[c(old.node, tip.node)[j]] <- gsub(".*\\[&|\\].*", "", newickext[j])
              }
            }
          }
          
        } else if (length(coloncommasplit) == 1) { # sister to another internal node or is the root
          
          if (extnewick_exists[i]) {
            node_annotations[old.node] <- gsub(".*\\[&|\\].*", "", newickext)
          }
        } else {
          stop("the format of this internal node string is not what we expect")
        }
      } else { # root
        node_bls[old.node] <- 0
        if (extnewick_exists[i]) {
          node_annotations[old.node] <- gsub(".*\\[&|\\].*", "", newickext)
        }
      }
      
    } else if (has_subtendingbranchs[i]) { # there must be one or two tips descending from the current node
      
      coloncommasplit <- strsplit(newick, ":|,")[[1]]
      if (!length(coloncommasplit) %in% c(2, 4)) {
        stop("the format of this tip node string is not what we expect")
      }
      
      tip.node_old <- tip.node
      for (j in 1:(length(coloncommasplit) / 2)) {
        tip.node <- tip.node + 1L
        tip.names[tip.node] <- coloncommasplit[j * 2 - 1]
        node_bls[tip.node] <- as.numeric(coloncommasplit[j * 2])
        
        edges[edges_rownum, ] <- c(at.node, tip.node)
        edges_rownum <- edges_rownum + 1L
      }
      
      if (extnewick_exists[i]) { # associated node info exists
        newickext <- strsplit(newickext, ",(?![^[]*])", perl = T)[[1]]
        
        for (j in 1:length(newickext)) {
          if (grepl("\\[&", newickext[j])) {
            node_annotations[tip.node_old + j] <- gsub(".*\\[&|\\].*", "", newickext[j])
          }
        }
      }
      
    } else {
      stop("cannot understand some components of this newick string")
    }
  }
  
  if (any(is.na(node_bls))) {
    stop("some branch lengths cannot be found")
  }
  
  # format the phylo object
  edge.lengths <- numeric(n.edges)
  node_indices <- (1:n.nodes)[-(n.tips + 1)]
  edge.lengths[match(node_indices, edges[, 2])] <- node_bls[node_indices]
  
  my.tree <- list(edge = edges, tip.label = tip.names, edge.length = edge.lengths, Nnode = n.tips - 1, node.data = node_annotations)
  if (node_bls[n.tips + 1] > 0) {
    my.tree$root.edge <- node_bls[n.tips + 1]
  }
  
  return(my.tree)
}

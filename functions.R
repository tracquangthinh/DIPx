fn.gsum <- function(gene, xdat) {
  ##  sum(expression) of a set of genes,
  pick <- gene %in% rownames(xdat)
  gene <- gene[pick]
  if (length(gene) > 1) gsum <- apply(xdat[gene, drop = FALSE, ], 2, sum)
  if (length(gene) == 1) gsum <- xdat[gene, ]
  if (length(gene) == 0) gsum <- rep(0, ncol(xdat))
  return(gsum)
}

generate_pas <- function(gex_path, synergy_path, target_updown_path, geneset_path, driver_gene_path) {

  load(target_updown_path)
  load(synergy_path)
  load(gex_path)
  load(geneset_path)
  load(driver_gene_path)

  list_driver_gene_names <- names(list_driver_gene)

  ### get drug target networks
  down <- Map(setdiff, tgt_down_all, tgt_up_all)
  down <- Map(setdiff, down, target)
  onup <- Map(union, target, tgt_up_all)

  # get up/down sets from drug target networks and pathways
  gs_onup = foreach(i=1:length(onup)) %dopar% {
    res = sapply(geneset,intersect, onup[[i]]) 
    res
  }
  names(gs_onup) = names(onup)

  gs_down = foreach(i=1:length(down)) %dopar% {
    res = sapply(geneset,intersect, down[[i]])
    res
  }
  names(gs_down) = names(down)

  driver_gene <- unique(unlist(list_driver_gene))

  all_drugs <- names(onup)

  res = foreach(i = 1:length(all_drugs)) %dopar% {
    tryCatch({
      drug <- all_drugs[i]
      cls <- synergy_df$CELL_LINE[which(synergy_df$COMBINATION_ID == drug)]

      # upstream pas
      ingene <- gs_onup[[drug]]
      gex_cl <- gex[, cls]
      insum.gex <- t(sapply(ingene, fn.gsum, gex_cl))
      insum.gex[is.na(insum.gex)] <- 0
      onup.gex.t <- insum.gex
      if (nrow(onup.gex.t) == 1) {
        onup.gex.t[1, ] <- as.numeric(onup.gex.t[1, ])
        onup.gex.t <- t(onup.gex.t)
      } else {
        onup.gex.t <- apply(insum.gex, 2, as.numeric)
      }

      # down pas
      ingene <- gs_down[[drug]]
      gex_cl <- gex[, cls]
      insum.gex <- t(sapply(ingene, fn.gsum, gex_cl))
      insum.gex[is.na(insum.gex)] <- 0
      down.gex.t <- insum.gex
      if (nrow(down.gex.t) == 1) {
        down.gex.t[1, ] <- as.numeric(down.gex.t[1, ])
        down.gex.t <- t(down.gex.t)
      } else {
        down.gex.t <- apply(insum.gex, 2, as.numeric)
      }

      # get driver genes of the corresponding cell lines
      drivers <- list_driver_gene[cls]
      gex_drivers <- gex_cl[
                      rownames(gex_cl) %in% unlist(drivers), drop = FALSE, ]

      for (j in 1:ncol(gex_drivers)) {
        pick <- rownames(gex_drivers) %in% drivers[[j]]
        # set zero for the expression of genes which are 
        # not drivers of the cell line
        gex_drivers[!pick, j] <- 0
      }

      # pas3: intersection between driver genes and pathways
      ingene <- geneset
      insum.gex <- t(sapply(ingene, fn.gsum, gex_drivers))
      insum.gex[is.na(insum.gex)] <- 0
      pas3.gex.t <- insum.gex
      if (nrow(pas3.gex.t) == 1) {
        pas3.gex.t[1, ] <- as.numeric(pas3.gex.t[1, ])
        pas3.gex.t <- t(pas3.gex.t)
      } else {
        pas3.gex.t <- apply(pas3.gex.t, 2, as.numeric)
      }


      colnames(onup.gex.t) <- colnames(down.gex.t) <-
      colnames(pas3.gex.t) <- paste0(cls, ".", drug)

      list(onup.gex.t, down.gex.t, pas3.gex.t)
    }, error = function(cond) {
      print(i)
      return(NULL)
    })
  }

  pas_onup <- lapply(res, function(r) r[[1]])
  pas_onup <- do.call(cbind, pas_onup)
  rownames(pas_onup) <- names(geneset)

  pas_down <- lapply(res, function(r) r[[2]])
  pas_down <- do.call(cbind, pas_down)
  rownames(pas_down) <- names(geneset)

  pas_driver <- lapply(res, function(r) r[[3]])
  pas_driver <- do.call(cbind, pas_driver)
  rownames(pas_driver) <- names(geneset)

  path = "data/created_pas.RData"
  save(pas_onup, pas_down, pas_driver,
       file = path)

  path
}


get_pas_matrix <- function(pas_path, synergy_path, geneset_path) {
  load(pas_path)
  load(synergy_path)
  load(geneset_path)

  pathways = names(geneset)

  df <- t(pas_onup[pathways, ])
  colnames(df) <- paste0("ONUP_", colnames(df))
  x <- as.data.frame(df)

  df <- t(pas_down[pathways, ])
  colnames(df) <- paste0("DOWN_", colnames(df))
  x <- cbind(x, df)

  df <- t(pas_driver[pathways, ])
  colnames(df) <- paste0("PAS3_", colnames(df))
  x <- cbind(x, df)

  x <- log2(x + 1)

  pick <- match(rownames(x), synergy_df$experiment)
  synergy_df <- synergy_df[pick, ]
  x$SYNERGY_SCORE <- synergy_df$SYNERGY_SCORE

  return(x)
}


#' Identify desired neighbors and calculate weighted expression levels
#'
#' This function calculates the weighted expression level of a gene to spots lying at a given distance.
#'
#' @param neighbor_dist A numeric input that indicates the distance between desired neighboring pairs
#' @param round.dist The rounded distance matrix
#' @param ST.mat The ST matrix
#' @return A matrix
#' @export
getNeighbor <- function(neighbor_dist,
                        round.dist,
                        ST.mat){
  round.dist <- t(apply(round.dist == neighbor_dist, 1, probalized))
  res <- round.dist %*% ST.mat
  return(res)
}


#' Calculte a part of SPER: spatial weighted expression level
#'
#' This function calculates the weighted expression level of a gene to spots lying at a given distance.
#'
#' @param round.dist The rounded distance matrix
#' @param ST.mat The ST matrix
#' @param dist.list A vector containing the list of desired distance values
#' @return A list of matrix
#' @export
getSPERmat <- function(round.dist,
                       ST.mat,
                       dist.list){
  res <- lapply(dist.list,
                getNeighbor,
                round.dist = round.dist,
                ST.mat = ST.mat)
  return(res)
}



#' Calculate the SPER score list
#'
#' This function calculates the SPER score list with the input of the distance matrix, the ST data and
#' the cell-type spatial CoDa.
#'
#' @param dist.mat  Distance matrix after rounded
#' @param dist.list Distance list which includes the value of interest in the distance matrix
#' @param ST.mat    The ST matrix
#' @param CoDa.data Cell-type spatial compositional data
#' @param gene.list Gene ID list; to specify SPER running on a given list of genes instead of all genes in the spatial transciptomics matrix.
#' @return A list of matrix, each representing the SPER subscores at a given distance
#' @export
SPER <- function(dist.mat,
                 dist.list,
                 ST.mat,
                 CoDa.data,
                 gene.list = NULL){
  if(!is.matrix(CoDa.data)){
    CoDa.data <- as.matrix(CoDa.data)
  }
  CoDa.data <- apply(CoDa.data, 2, function(x){x / sum(x)})
  if(!is.null(gene.list)){
    ST.mat <- as.matrix(ST.mat[,gene.list])
  }
  pair_cor_list <- getSPERmat(dist.mat,
                              ST.mat,
                              dist.list)
  pair.list <- list()
  for(i in 1:length(pair_cor_list)){
    pair.list[[i]] <- t(CoDa.data) %*% pair_cor_list[[i]]
  }

  reshape.pair.list <- list()
  for(i in 1:ncol(CoDa.data)){
    ##  Initialization
    reshape.pair.list[[i]] <- matrix(0, length(pair.list), ncol(pair.list[[1]]))
    rownames(reshape.pair.list[[i]]) <- dist.list
    colnames(reshape.pair.list[[i]]) <- colnames(pair.list[[1]])

    ##  Reshape the data: extract rows from original data list, forming new dfs based on cell type
    for(j in 1:length(pair.list)){
      reshape.pair.list[[i]][j,] <- pair.list[[j]][i,]
    }
    ##  Scale data
    reshape.pair.list[[i]] <- t(t(reshape.pair.list[[i]]) / colMeans(ST.mat))
  }
  names(reshape.pair.list) <- colnames(CoDa.data)
  return(reshape.pair.list)
}

#' Weight SPER list to one score matrix
#'
#' Apply a weight to the SPER list and integrate the subscores at different distance into one matrix
#'
#' @param reshape.pair.list SPER list of subscore matrices
#' @param weight            A weight to apply. Must have the same length as the distance list.
#' @return A matrix containing the SPER scores of all pairs
#' @export
weightSPER <- function(reshape.pair.list,
                       weight){
  res <- sapply(reshape.pair.list, function(X){weight %*% X})
  rownames(res) <- colnames(reshape.pair.list[[1]])
  return(res)
}

#' Find the putative paracrine signals based on SPER scores
#'
#' This function helps find the potential signals based on the given threshold in both SPER scores and
#' expression prevalence. One can provide the marker (cell-type-specific) gene list and ligand-receptor pair data to
#' limit the search range and get a more precise result.
#'
#' @param target.cell     Target ell Type Name
#' @param score.mat       Score matrix (SPER, or correlations)
#' @param gene.set        The gene list of interest: could be extracellular or ligand gene sets
#' @param expr.frac.mat   Expression prevalence matrix: matrix values from 0 to 1, showing the percentage of cells in that cell type expressing the gene.
#' @param LRP.data        Ligand-receptor pair data; used to find the receptors of found paracrine signals; should contain at least two columns: 'gene1' for the ligands and 'gene2' for receptors
#' @param marker.gene     The marker gene matrix with a column 'Type' indicating the cell type information and 'Gene' indicating the gene ID
#' @param score.threshold Threshold on the score: pairs whose scores larger than the threshold will be kept
#' @param expr.fraction   Threshold on expression prevalence: pairs whose prevalence smaller than the threshold will be kept
# #' @param method.name     Name of the method; default is SPER
#' @return A gene list containing the spatial dependent genes; if LRP data provided, return a data frame containing both ligand and receptor information
#' @export
findSPERsignals <- function(target.cell,
                            score.mat,
                            expr.frac.mat,
                            marker.gene = NULL,
                            gene.set = NULL,
                            LRP.data = NULL,
                            score.threshold = 1.5,
                            expr.fraction = 0.15
                            # only.extra = T,
                            # plot.feature = T,
                            # return.numbers = F,
                            # calc.receptor.frac = F,
                            # except.marker = T,
                            # method.name = "SPER"
                            ){

  feature.gene.list <- rownames(score.mat)[which(score.mat[,target.cell] > score.threshold)]

  ##  First, remove the marker genes from our candidates
  if(!is.null(marker.gene)){
    marker_id <- intersect(marker.gene$Gene[which(marker.gene$Type == target.cell)], rownames(score.mat))
    feature.gene.list <- setdiff(feature.gene.list, marker_id)
  }

  ##  Return if no SDG is found
  if(length(feature.gene.list) == 0){
    print("No feature gene found.")
    return()
  }

  ##  If scRNA-seq expression matrix provided, find SDG only below given expression fraction threshold
  if(!is.null(expr.frac.mat)){
    expr_frac_filtered_gene <- rownames(expr.frac.mat)[which(expr.frac.mat[,target.cell] < expr.fraction)]
    feature.gene.list <- intersect(feature.gene.list, expr_frac_filtered_gene)
    ##  Return if no SDG is found
    if(length(feature.gene.list) == 0){
      print("No feature gene found.")
      return()
    }
  }

  ##  If interested gene set provided, find SDG only within the set
  if(!is.null(gene.set)){
    feature.gene.list <- intersect(feature.gene.list, gene.set)
    ##  Return if no SDG is found
    if(length(feature.gene.list) == 0){
      print("No feature gene found.")
      return()
    }
  }

  ##  If LRP data provided, find SDG only in ligand list
  if(!is.null(LRP.data)){
    feature.gene.list <- intersect(feature.gene.list, LRP.data$gene1)

    ##  Return if no SDG is found
    if(length(feature.gene.list) == 0){
      print("No putative paracrine gene found.")
      return()
    }
    ##  Furthermore, show the expression of receptors in target cells: maybe in another function called 'findSDGreceptor'
    receptor_list <- LRP.data[which(LRP.data$gene1 %in% feature.gene.list), ]
    tmp_expr_frac_mat <- expr.frac.mat[,target.cell]
    names(tmp_expr_frac_mat) <- rownames(expr.frac.mat)
    filtered_LRP <- matrix(0, 1, 3)
    for(i in 1:nrow(receptor_list)){
      receptor_ID <- receptor_list[i, "gene2"]

      if(!(receptor_ID %in% rownames(expr.frac.mat))){next}

      if(tmp_expr_frac_mat[receptor_ID] > 0.05){
        filtered_LRP <- rbind(filtered_LRP,
                              c(receptor_list[i, "gene1"],
                                receptor_list[i, "gene2"],
                                round(tmp_expr_frac_mat[receptor_ID], 4)))
      }
    }
    if(nrow(filtered_LRP) == 1){
      return("No putative paracrine ligand found.")
    }
    filtered_LRP <- as.data.frame(filtered_LRP)[-1,]
    colnames(filtered_LRP) <- c("SPER_ligand", "SPER_receptor", "Receptor_frac")
    filtered_LRP$SPER_score <- round(score.mat[filtered_LRP$SPER_ligand, target.cell], 4)
    filtered_LRP <- filtered_LRP[order(filtered_LRP[,"SPER_score"], filtered_LRP[,"Receptor_frac"], decreasing = T),]

    SPER_sorted_name <- rownames(score.mat)[order(score.mat[, target.cell], decreasing = T)]
    SPER_rank <- which(SPER_sorted_name %in% filtered_LRP$SPER_ligand)
    names(SPER_rank) <- SPER_sorted_name[SPER_rank]
    for(i in 1:nrow(filtered_LRP)){
      filtered_LRP$Score_rank[i] <- SPER_rank[filtered_LRP$SPER_ligand[i]]
    }

    rownames(filtered_LRP) <- 1:nrow(filtered_LRP)
    filtered_LRP$Cell_type <- target.cell
    return(filtered_LRP)
  }

  ##  Return the SDG gene list
  return(feature.gene.list)
}

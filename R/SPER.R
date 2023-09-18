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

####  Plot should be in another function, or just another command after getting this SDG list
#' Find the spatial dependent gene (SDG) based on SPER scores
#'
#' This function helps find the SDG based on the given threshold in both SPER scores and
#' expression prevalence.
#'
#' @import dplyr
#' @param cell.type       Cell Type Name
#' @param score.mat       Score matrix (SPER, or correlations)
#' @param extra.gene.list The gene list of interest: could be extracellular or ligand gene sets
# #' @param expr.mat        Expression prevalence matrix
#' @param LRP.data        Ligand-receptor pair data; used to find the receptors of found paracrine signals
#' @param marker.gene     The marker gene matrix with a column 'Type' indicating the cell type information and 'Gene' indicating the gene ID
#' @param ratio.threshold Threshold on the score: pairs whose scores larger than the threshold will be kept
#' @param expr.fraction   Threshold on expression prevalence: pairs whose prevalence smaller than the threshold will be kept
#' @param only.extra      Boolean value: only consider the candidate genes
# #' @param plot.feature    Boolean value: plot the SDG using 'plotFeatureGene'
# #' @param return.numbers  Boolean value: return the statistical results of hypergeometric test
# #' @param calc.receptor.frac Boolean value: calculate the fraction of receptors that expressed in the target cell types. A column of these receptors must be provided in the extra.gene.list as 'gene2'.
# #' @param except.marker   Boolean value: ignore the marker of the cell types
#' @param method.name     Name of the method; default is SPER
#' @return A gene list containing the spatial dependent genes
#' @export
findSDG <- function(cell.type,
                    score.mat,
                    extra.gene.list,
                    # expr.mat,
                    LRP.data,
                    marker.gene,
                    ratio.threshold = 1.5,
                    expr.fraction = 0.15,
                    only.extra = T,
                    # plot.feature = T,
                    # return.numbers = F,
                    # calc.receptor.frac = F,
                    # except.marker = T,
                    method.name = "SPER"
                    ){
  marker_id <- intersect(marker.gene$Gene[which(marker.gene$Type == cell.type)], rownames(score.mat))
  feature.gene.list <- setdiff(rownames(score.mat)[which(score.mat[,cell.type] > ratio.threshold)],
                               marker_id)
  # feature.gene.list <- intersect(feature.gene.list,
  #                                rownames(expr_frac_matrix)[which(expr_frac_matrix[,cell.type] < expr.fraction)])
  if(length(feature.gene.list) == 0){
    print("No feature gene found.")
    #if(return.numbers){return(c(1,0,0))}
    #if(calc.receptor.frac){return(c(1,0,0))}
    return(NULL)
  }
  if(only.extra){
    feature.gene.list <- intersect(feature.gene.list, extra.gene.list)
  }

  ##  Plot the feature gene's spatial plot
  # if(plot.feature){
  #   plot.dir <- paste0("./Plots/Feature Gene/", method.name, "/", cell.type, "/")
  #   dir.create(plot.dir, recursive = T, showWarnings = F)
  #   tmp <- sapply(feature.gene.list,
  #                 plotSDG,
  #                 cell.type = cell.type,
  #                 plot.dir = plot.dir)
  # }


  # res <- stats::phyper(length(extra.feature.gene.list) - 1,
  #                      length(extra.gene.list),
  #                      nrow(score.mat) - length(extra.gene.list),
  #                      length(feature.gene.list),
  #                      lower.tail = FALSE)

  ##  Return both p-value and receptor fraction
  # if(calc.receptor.frac){
  #   if(is.null(LRP.data$gene2)){
  #     print("Please provide receptor gene as 'gene2'.")
  #     return(res)
  #   }
  #   LRP.data <- LRP.data %>%
  #     filter(gene2 %in% rownames(score.mat)) %>%
  #     filter(gene1 %in% rownames(score.mat))
  #
  #   tmp <- LRP.data %>%
  #     filter(gene1 %in% extra.feature.gene.list)
  #   potential.receptor <- unique(tmp$gene2)
  #   potential.receptor <- potential.receptor[which(expr.mat[potential.receptor,cell.type] > 0.5)]
  #   tmp2 <- tmp %>%
  #     filter(gene2 %in% potential.receptor)
  #
  #   potential.receptor <- unique(LRP.data$gene2)
  #   potential.receptor <- potential.receptor[which(expr.mat[potential.receptor,cell.type] > 0.5)]
  #   tmp2_remain <- LRP.data %>%
  #     filter(gene2 %in% potential.receptor)
  #
  #   if(nrow(tmp2) == 0){return(c(1,0,0))}
  #   res <- stats::phyper(nrow(tmp2) - 1,
  #                        nrow(tmp2_remain),
  #                        nrow(LRP.data) - nrow(tmp2_remain),
  #                        nrow(tmp),
  #                        lower.tail = FALSE)
  #   return(c(res, nrow(tmp2), nrow(tmp)))
  #
  # }

  ##  Return both p-value and number of catch
  # if(return.numbers){
  #   return(c(res, length(extra.feature.gene.list), length(feature.gene.list)))
  # }


  ##  Return just the p-value
  return(feature.gene.list)
}

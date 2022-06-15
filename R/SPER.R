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
#' @return A list of matrix, each representing the SPER subscores at a given distance
#' @export
SPER <- function(dist.mat,
                 dist.list,
                 ST.mat,
                 CoDa.data){
  if(!is.matrix(CoDa.data)){
    CoDa.data <- as.matrix(CoDa.data)
  }
  CoDa.data <- apply(CoDa.data, 2, function(x){x / sum(x)})
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
  return(res)
}



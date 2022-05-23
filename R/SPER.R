#' Identify desired neighbors and calculate weighted expression levels
#'
#' This function calculates the weighted expression level of a gene to spots lying at a given distance.
#'
#' @param neighbor_dist A numeric input that indicates the distance between desired neighboring pairs
#' @param round.dist The rounded distance matrix
#' @param expr.mat The expression matrix
#' @return A matrix
#' @export
getNeighbor <- function(neighbor_dist,
                        round.dist,
                        expr.mat){
  round.dist <- t(apply(round.dist == neighbor_dist, 1, probalized))
  res <- round.dist %*% expr.mat
  return(res)
}


#' Calculate the SPER for a distance
#'
#' This function calculates the weighted expression level of a gene to spots lying at a given distance.
#'
#' @param round.dist The rounded distance matrix
#' @param expr.mat The expression matrix
#' @param dist.list A vector containing the list of desired distance values
#' @return A list of matrix
#' @export
getSPER <- function(round.dist,
                    expr.mat,
                    dist.list){
  res <- lapply(dist.list, 
                getNeighbor, 
                round.dist = round.dist,
                expr.mat = expr.mat)
  return(res)
}
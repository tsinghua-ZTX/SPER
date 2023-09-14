#' Normalize a vector: L1 norm
#'
#' This function normalizes a given vector so that the sum of every elements is 1
#'
#' @param x An Input vector
#' @return Normalized vector
#' @export
probalized <- function(x){
  if(sum(x) == 0){
    return(x)
  }else{
    return(x / sum(x))
  }
}


#' Normalize a vector: L2 norm
#'
#' This function normalizes a given vector so that the length (L2) of a vector is 1
#'
#' @param x An Input vector
#' @return Normalized vector
#' @export
normalizd <- function(x){
  if(sum(x^2) == 0){
    return(x)
  }else{
    return(x / sqrt(sum(x^2)))
  }
}


#' Scale the range of a list of numbers from 0 to 1
#'
#' Return the scaled vector whose range is from 0 to 1
#'
#' @param x An Input vector
#' @return Scaled vector
#' @export
range_1 <- function(x){(x-min(x))/(max(x)-min(x))}




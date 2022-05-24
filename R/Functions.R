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



#' Find the spatial dependent gene (SDG) based on SPER scores
#'
#' This function helps find the SDG based on the given threshold in both SPER scores and 
#' expression prevalence. 
#' 
#' @import dplyr
#' @param cell.type       Cell Type Name
#' @param score.mat       Score matrix (SPER, or correlations)
#' @param extra.gene.list The gene list of interest: could be extracellular or ligand gene sets
#' @param expr.mat        Expression prevalence matrix
#' @param ratio.threshold Threshold on the score: pairs whose scores larger than the threshold will be kept
#' @param expr.fraction   Threshold on expression prevalence: pairs whose prevalence smaller 
#' than the threshold will be kept
#' @param only.extra      Boolean value: only consider the candidate genes
#' @param plot.feature    Boolean value: plot the SDG using 'plotFeatureGene'
#' @param return.numbers  Boolean value: return the statistical results of hypergeometric test
#' @param calc.receptor.frac Boolean value: calculate the fraction of receptors that expressed in the target
#' cell types. A column of these receptors must be provided in the extra.gene.list as 'gene2'.
#' @param except.marker   Boolean value: ignore the marker of the cell types
#' @param marker.gene     The marker gene matrix with a column 'Type' indicating the cell type information and
#' 'Gene' indicating the gene ID
#' @param LRP.data        Ligand-receptor pair data
#' @param method.name     Name of the method
#' @return Scaled vector
#' @export
findSDG <- function(cell.type,
                    score.mat,
                    extra.gene.list = gene_id,
                    expr.mat = expr_frac_matrix,
                    ratio.threshold = 1.5,
                    expr.fraction = 0.15,
                    method.name = "Multiply",
                    only.extra = T,
                    plot.feature = T,
                    return.numbers = F,
                    calc.receptor.frac = F,
                    LRP.data = LRP_data,
                    except.marker = T,
                    marker.gene = marker_gene){
  marker_id <- intersect(marker.gene$Gene[which(marker.gene$Type == cell.type)], rownames(score.mat))
  feature.gene.list <- setdiff(rownames(score.mat)[which(score.mat[,cell.type] > ratio.threshold)], 
                               marker_id)
  feature.gene.list <- intersect(feature.gene.list, 
                                 rownames(expr_frac_matrix)[which(expr_frac_matrix[,cell.type] < expr.fraction)])
  if(length(feature.gene.list) == 0){
    print("No feature gene found.")
    if(return.numbers){return(c(1,0,0))}
    if(calc.receptor.frac){return(c(1,0,0))}
    return(1)
  }
  if(only.extra){
    extra.feature.gene.list <- intersect(feature.gene.list, extra.gene.list)
  }
  res <- stats::phyper(length(extra.feature.gene.list) - 1, 
                       length(extra.gene.list), 
                       nrow(score.mat) - length(extra.gene.list),
                       length(feature.gene.list),
                       lower.tail = FALSE)
  ## Return the odds ratio if wanted
  # res2 <- (length(extra.feature.gene.list) / length(feature.gene.list)) / 
  #   (length(extra.gene.list) / nrow(score.mat))

  ##  Plot the feature gene's spatial plot
  if(plot.feature){
    plot.dir <- paste0("Plots/Feature Gene/", method.name, "/", cell.type, "/")
    dir.create(plot.dir, recursive = T, showWarnings = F)
    if(only.extra){
      tmp <- sapply(extra.feature.gene.list, 
                    plotFeatureGene,
                    cell.type = cell.type,
                    plot.dir = plot.dir)
    }else{
      tmp <- sapply(feature.gene.list, 
                    plotFeatureGene,
                    cell.type = cell.type,
                    plot.dir = plot.dir)
    }
  }
  
  ##  Return both p-value and receptor fraction
  if(calc.receptor.frac){
    if(is.null(LRP.data$gene2)){
      print("Please provide receptor gene as 'gene2'.")
      return(res)
    }
    LRP.data <- LRP.data %>%
      filter(gene2 %in% rownames(score.mat)) %>%
      filter(gene1 %in% rownames(score.mat))
    
    tmp <- LRP.data %>%
      filter(gene1 %in% extra.feature.gene.list)
    potential.receptor <- unique(tmp$gene2)
    # tmp_id <- intersect(potential.receptor, rownames(expr.mat))
    potential.receptor <- potential.receptor[which(expr.mat[potential.receptor,cell.type] > 0.5)]
    tmp2 <- tmp %>%
      filter(gene2 %in% potential.receptor)
    
    
    potential.receptor <- unique(LRP.data$gene2)
    potential.receptor <- potential.receptor[which(expr.mat[potential.receptor,cell.type] > 0.5)]
    tmp2_remain <- LRP.data %>%
      filter(gene2 %in% potential.receptor)
    
    if(nrow(tmp2) == 0){return(c(1,0,0))}
    res <- stats::phyper(nrow(tmp2) - 1, 
                         nrow(tmp2_remain), 
                         nrow(LRP.data) - nrow(tmp2_remain),
                         nrow(tmp),
                         lower.tail = FALSE)
    return(c(res, nrow(tmp2), nrow(tmp)))
    
    # max_rep_frac <- rep(0, length(extra.feature.gene.list))
    # for(K in 1:length(extra.feature.gene.list)){
    #   tmp <- LRP.data %>%
    #     filter(gene1 == extra.feature.gene.list[K])
    #   possible_rep <- intersect(rownames(expr.mat), tmp$gene2)
    #   if(length(possible_rep) == 0){
    #     max_rep_frac[K] <- 0
    #     next
    #   }
    #   rep_expr <- expr.mat[possible_rep, cell.type]
    #   max_rep_frac[K] <- max(rep_expr)
    # }
    # return(c(res, mean(max_rep_frac)))
  }
  
  ##  Return both p-value and number of catch
  if(return.numbers){
    return(c(res, length(extra.feature.gene.list), length(feature.gene.list)))
  }
  
  
  ##  Return just the p-value
  return(res)
}
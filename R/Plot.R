#' Plot the spatial figures for genes / cell types
#'
#' With spatial data, coordinates and given specific feature (gene or cell type), this function plots the  
#' map of this feature's spatial distribution.
#'
#' @param feature Specific spatial feature to plot
#' @param spatial_df A data.frame with rows as spatial spots and columns as each features
#' @param coord Spatial coordinates that match the spatial data
#' @return a ggplot2 object
#' @export
gg_gene <- function(feature, spatial_df, coord, type = "Estimated"){
  p <- ggplot() + 
    geom_point(aes(x = coord$imagecol, 
                   y = coord$imagerow, 
                   col = spatial_df[,feature])) + 
    scale_color_gradient2(low = "darkblue", 
                          high = "yellow", 
                          mfeature = "purple", 
                          mfeaturepoint = max(spatial_df[,feature], na.rm = T) / 2) +
    theme_bw() + 
    labs(x = "X / μm", 
         y = "Y / μm", 
         title = paste0(feature, " ", type,  " Spatial Distribution"), 
         col = feature) +
    theme(aspect.ratio = 1)
  return(p)
}

#' Visualize SPER scores of a gene/cell-type pair
#'
#' Draw the plot of: 1) the SPER curve of the indicated pair; 2) the expression level of the gene in scRNA-seq
#' reference data; 3) spatial map using gg_gene
#'
#' @param gene.id   Gene ID
#' @param cell.type Name of cell type
#' @param plot.dir  Directory to save the plots
#' @param dist.list Distance list when plotting the SPER curve
#' @param reference scRNA-seq reference; should be a Seurat object
#' @return Null
#' @export
visualizePairCor <- function(gene.id,
                             cell.type,
                             plot.dir,
                             dist.list = dist_list,
                             reference = allen_reference){
  dir.create(plot.dir, showWarnings = F, recursive = T)
  jpeg(paste0(plot.dir, "/", gene.id, "_", cell.type, "_lineplot.jpg"), quality = 100)
  plot(x = dist.list, 
       y = rf_pair_cor_list[[cell.type]][,gene.id], 
       type = "b",
       main = paste0(gene.id, " : ", cell.type),
       ylab = "Average Expression Densitty",
       xlab = "Distance / μm")
  dev.off()
  p1 <- VlnPlot(reference,
                features = gene.id,
                pt.size = 0.2,
                ncol = 1)
  Sys.sleep(0.5)
  ggsave(paste0(plot.dir, "/", gene.id, "_expression.jpg"), p1, width = 8, height = 5)
  p2 <- gg_gene(gene.id, ST_expression, coordinates)
  Sys.sleep(0.5)
  ggsave(paste0(plot.dir, "/", gene.id, "_ST.jpg"), p2, width = 5, height = 5)
  return(NULL)
}

#' Plot the spatial dependent gene (SDG) based on SPER scores
#'
#' Generate a series of plots of SDGs discovered by 'findSDG': the SPER curve is plotted, and a spatial map
#' showing the spatial of both the gene and the cell type is generated. 
#'
#' @param cell.type       Cell Type Name
#' @param gene.id         Gene ID
#' @param score           SPER score shown in the plot title
#' @param thershold       Threshold for the proportion of given cell type; spots whose value larger than this
#' threshold will be marked in the spatial map
#' @param coord           Spatial coordinates for 'gg_gene'
#' @param ST.data         Spatial transcriptomics data, used for the spatial map of the gene's expression
#' @param CoDa.data       Cell-type spatial compositional data, used for the spatial map of cell-type 
#' distribution
#' @param reference.data  scRNA-seq reference; should be a Seurat object
#' @param SPER.data       SPER list to plot the SPER curve
#' @param plot.dir        Directory to save the plots
#' @return Null
#' @export
plotSDG <- function(cell.type, 
                    gene.id, 
                    score = NULL, 
                    thershold = 0.3,
                    coord = coordinates_test,
                    ST.data = ST_expression,
                    CoDa.data = KP_local$Z,
                    reference.data = allen_reference,
                    SPER.data = rf_pair_cor_list,
                    plot.dir = "Plots/Feature Gene/"){
  # gg_gene(gene.id, ST.data, coordinates)
  # gg_gene(cell.type, CoDa.data, coordinates)
  pdf(paste0(plot.dir, gene.id, " curve.pdf"), width = 5, height = 4)
  plot(x = dist_list, 
       y = SPER.data[[cell.type]][,gene.id], type = "b", 
       xlab = "Distance / um",
       ylab = "Expression Ratio",
       main = paste0(cell.type, ": ", gene.id))
  dev.off()
  
  VlnPlot(reference.data,
          features = gene.id,
          pt.size = 0.2,
          ncol = 1)
  ggsave(paste0(plot.dir, gene.id, " expr.pdf"), width = 8, height = 5)
  int_signal <- (ST.data[,gene.id] > 0) + as.numeric(CoDa.data[,cell.type] > thershold) * 2
  p1 <- ggplot() + 
    geom_point(aes(x = coord$imagecol, 
                   y = coord$imagerow, 
                   col = as.factor(int_signal))) + 
    theme_bw() + 
    scale_color_manual(name = "Spot Class", 
                       breaks = 0:3,
                       values = c("gray", "darkred", "darkblue", "purple"), 
                       labels = c("Other", paste0(gene.id), 
                                  paste0(cell.type, " > 30%"), "Overlap")) + 
    labs(x = "X / um", 
         y = "Y / um", 
         title = paste0(gene.id, "/", cell.type,": ", score)) +
    theme(aspect.ratio = 1)
  ggsave(paste0(plot.dir, gene.id, " spatial.pdf"), p1, width = 6, height = 5)
  return()
}
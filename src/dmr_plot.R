dmr_plot <- function(dmr_gr, meth_data, anno_gr, meth_groups, important_cpgs = NULL, flanks = NULL, title = NULL, legend = T, highlight = NULL, rug = T, gridlines = F){
  require(reshape2)
  
  if(class(dmr_gr) != "GRanges") stop("dmr_gr must be a GRanges object")
  if(is.null(flanks)) flanks <- width(dmr_gr)/2
  if(!is.null(highlight)){
    if(length(highlight) == 1){
      print("Highlight drawn around a single point")
      highlight <- c(highlight, highlight)
    } else if(length(highlight) > 2){
      stop("Highlight can be one or two values")
    }
    highlight <- sort(highlight)
    highlight <- data.frame(x0 = highlight[1], x1 = highlight[2])
  }
  
  plotrange <- range(dmr_gr + flanks)
  
  cg_ids <- queryHits(findOverlaps(query = anno_gr, subject = dmr_gr))
  cg_ids <- c(min(cg_ids)-1, cg_ids, max(cg_ids)+1)
  anno_sub <- anno_gr[cg_ids,]
  
  meth_df <- data.frame(pos = start(anno_sub), 
                        meth_data[names(anno_sub),],
                        label = names(anno_sub))
  
  meth_df_melt <- melt(meth_df, id = c("pos", "label"), variable.name = "sample", value.name = "methylation")
  meth_df_melt$Group <- rep(meth_groups, each = length(cg_ids))
  
  if(!is.null(important_cpgs)){
    #Dilation factor of the highlighted region
    dil_factor <- ceiling(width(plotrange)/1000)
    if(all(!important_cpgs %in% names(anno_gr))) stop("important_cpgs cannot be found in anno_gr")
    if(all(!important_cpgs %in% meth_df_melt$label)) stop("important_cpgs cannot be found in dmr_gr")
    important_cpgs_coords <- data.frame(imp_cpgs_start = start(anno_gr[important_cpgs,])-dil_factor,
                                        imp_cpgs_end = start(anno_gr[important_cpgs,])+dil_factor)
  } 
  
  if(is.null(title)) title <- "Methylation"
  
  #Plot object
  gplot_obj <- ggplot(meth_df_melt, aes(x = pos, y = methylation))
  
  if(!is.null(highlight)){
    #Dilation factor of the highlighted region
    dil_factor <- ceiling(width(plotrange)/150)
    
    gplot_obj <- gplot_obj + 
      geom_rect(data = highlight,
                inherit.aes = FALSE,
                aes(xmin = x0 - dil_factor, xmax = x1 + dil_factor, ymin = -Inf, ymax = Inf), alpha = 0.15)
  }
  
  gplot_obj <- gplot_obj + 
    geom_point(aes(col = Group, group = Group)) +
    stat_summary(aes(col = Group, group = Group), fun.y = mean, geom = "smooth") +
    coord_cartesian(xlim = c(start(plotrange), end(plotrange))) +
    ylim(0,1) +
    ylab("Beta") +
    labs(title = title,
         subtitle = as.character(dmr_gr)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
  

  #Label CpGs of interest
  if(!is.null(important_cpgs)){
    gplot_obj <- gplot_obj + 
      geom_rect(data = important_cpgs_coords,
                inherit.aes = F,
                mapping = aes(xmin = imp_cpgs_start, 
                              xmax = imp_cpgs_end, 
                              ymin = 0, 
                              ymax = 1),
                alpha = 0.5)
  } 
  #Add rugplot at the bottom
  if(rug) gplot_obj <- gplot_obj + geom_rug(sides = "b")
  #Add gridline
  if(!gridlines) gplot_obj <- gplot_obj + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #Remove legend
  if(!legend) gplot_obj <- gplot_obj + theme(legend.position = "none")
  
  return(gplot_obj)
}

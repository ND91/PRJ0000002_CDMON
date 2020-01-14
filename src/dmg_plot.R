dmg_plot <- function(gene_of_interest, tophits_gr, ylim = NULL, significance = T, smooth = T){
  require(ggbio)
  require(Homo.sapiens)
  
  cpgs <- tophits_gr[grep(paste0("(^|;)", as.character(gene_of_interest), "($|;)"), tophits_gr$UCSC_RefGene_Name),] 
  
  #methylation difference
  mdiff_track <- ggplot(cpgs, aes(x = start, y = Beta)) +
    #geom_alignment() + 
    geom_hline(yintercept = 0) +
    #geom_line(alpha = 0.5) +
    
    theme_bw() +
    theme(legend.pos = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 12))
  
  if(!is.null(ylim)) mdiff_track <- mdiff_track + ylim(ylim)
  
  if(smooth){
    mdiff_track <- mdiff_track + geom_smooth(method = "loess", alpha = 0.5)
  } else{
    mdiff_track <- mdiff_track + geom_line(alpha = 0.5)
  }
  
  if(significance){
    cpgs$significance <- cpgs$adj.P.Val<0.05
    mdiff_track <- mdiff_track + geom_point(aes(alpha = -log10(P.Value), fill = significance), shape = 21)
  } else{
    mdiff_track <- mdiff_track + geom_point(aes(alpha = -log10(P.Value))) 
  }
  
  gene_track <- ggplot() + 
    geom_alignment(Homo.sapiens, which = range(cpgs, ignore.strand = T), columns = "GENEID") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 12))

  plotobj <- tracks("Gene" = gene_track, 
                    "Methylation difference" = mdiff_track, 
                    heights = c(3, 4))
}

dmg_plot <- function(gene_of_interest, tophits_gr, ylim = NULL, significance = T, smooth = T, stat = "reduce"){
  require(ggbio)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  data(genesymbol, package = "biovizBase")
  
  cpgs <- tophits_gr[grep(paste0("(^|;)", as.character(gene_of_interest), "($|;)"), tophits_gr$UCSC_RefGene_Name),] 
  
  #methylation difference
  mdiff_track <- ggplot(cpgs, aes(x = start, y = Beta)) +
    #geom_alignment() + 
    geom_hline(yintercept = 0) +
    #geom_line(alpha = 0.5) +
    theme_bw() +
    theme(legend.pos = "top",
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 12),
          panel.border = element_blank())
  
  if(!is.null(ylim)) mdiff_track <- mdiff_track + ylim(ylim)
  
  if(smooth){
    mdiff_track <- mdiff_track + geom_smooth(method = "loess")
  } else{
    mdiff_track <- mdiff_track + geom_line()
  }
  
  if(significance){
    cpgs$significance <- cpgs$adj.P.Val<0.05
    mdiff_track <- mdiff_track + geom_point(aes(fill = significance), shape = 21, alpha = 0.2)
  } else{
    mdiff_track <- mdiff_track + geom_point(alpha = 0.2) 
  }
  
  gene_track <- ggplot() + 
    #geom_alignment(Homo.sapiens, which = range(cpgs, ignore.strand = T), columns = "GENEID", stat = "reduce") + 
    geom_alignment(TxDb.Hsapiens.UCSC.hg19.knownGene, which = genesymbol[gene_of_interest], stat = stat) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 12),
          panel.border = element_blank())

  plotobj <- tracks("Gene" = gene_track, 
                    "Methylation difference" = mdiff_track, 
                    heights = c(1, 4))
}

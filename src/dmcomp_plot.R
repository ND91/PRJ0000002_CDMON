dmcomp_plot <- function(gene_of_interest, anno_gr, tophits1, tophits2, comparison1_name, comparison2_name, meth_data, meth_groups, ylim = NULL, smooth_diff = T, smooth_m = T, stat = "reduce"){
  require(ggbio)
  require(reshape2)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  data(genesymbol, package = "biovizBase")
  
  cpgs <- anno_gr[grep(paste0("(^|;)", as.character(gene_of_interest), "($|;)"), anno_gr$UCSC_RefGene_Name),] 
  
  #gene track
  gene_track <- ggplot() + 
    geom_alignment(TxDb.Hsapiens.UCSC.hg19.knownGene, which = genesymbol[gene_of_interest], stat = stat) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 12),
          panel.border = element_blank())
  
  #methylation difference track
  tophits_merged <- rbind(data.frame(tophits1[names(cpgs), c("Beta", "P.Value", "Name")], Comparison = comparison1_name),
                          data.frame(tophits2[names(cpgs), c("Beta", "P.Value", "Name")], Comparison = comparison2_name))
  tophits_merged$Pos <- start(anno_gr[tophits_merged$Name])
  
  mdiff_track <- ggplot(tophits_merged, aes(x = Pos, y = Beta, group = Comparison)) +
    #geom_alignment() + 
    geom_hline(yintercept = 0) +
    #geom_line(alpha = 0.5) +
    ylab("% Methylation difference") +
    scale_color_manual(values=c("#008000", "#ff00ff")) +
    theme_bw() +
    theme(legend.pos = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size = 12))
  
  if(!is.null(ylim)) mdiff_track <- mdiff_track + ylim(ylim)
  
  if(smooth_diff){
    mdiff_track <- mdiff_track + geom_smooth(method = "loess", alpha = 0.5, aes(col = Comparison))
  } else{
    mdiff_track <- mdiff_track + geom_line(alpha = 0.5, aes(col = Comparison))
  }
  
  #mdiff_track <- mdiff_track + geom_point(aes(alpha = -log10(P.Value), col = Comparison))
  mdiff_track <- mdiff_track + geom_point(aes(col = Comparison), alpha = 0.1) 
  
  #methylation track
  meth_df <- reshape2::melt(data.frame(Group = meth_groups, t(meth_data[names(cpgs),])), variable.name = "CpG", value.name = "Methylation")
  meth_df$Pos <- start(anno_gr[meth_df$CpG])
  
  m_track <- ggplot(meth_df, aes(x = Pos, y = Methylation)) +
    geom_point(aes(col = Group, group = Group), alpha = 0.1) +
    #coord_cartesian(xlim = c(start(plotrange), end(plotrange))) +
    #ylim(0,1) +
    scale_color_manual(values=c("#00dae0", "#0000ff", "#ff0000")) +
    ylab("% Methylation") +
    theme_bw() +
    theme(legend.pos = "bottom",
          plot.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if(smooth_m){
    m_track <- m_track + geom_smooth(method = "auto", alpha = 0.5, aes(col = Group))
  } else{
    m_track <- m_track + stat_summary(aes(col = Group, group = Group), fun.y = mean, geom = "smooth")
  }
  
  # Summary
  plotobj <- tracks("Gene" = gene_track, 
                    "% Methylation difference" = mdiff_track, 
                    "% Methylation" = m_track,
                    heights = c(1, 4, 4))
}

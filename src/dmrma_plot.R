dmrma_plot <- function(roi_gr, ownmeth_data, ownmeth_anno_gr, ownmeth_groups, pubmeth_data, pubmeth_anno_gr, pubmeth_groups){
  require(ggbio)
  require(Homo.sapiens)
  
  ownmeth_ol <- subjectHits(findOverlaps(roi_gr, ownmeth_anno_gr))
  
  if(length(ownmeth_ol) == 0) stop("No overlap found in own data")
  
  pubmeth_ol <- subjectHits(findOverlaps(roi_gr, pubmeth_anno_gr))
  
  if(length(pubmeth_ol) == 0) stop("No overlap found in public data")
  
  #gene track
  gene_track <- ggplot() + 
    geom_alignment(Homo.sapiens, which = range(roi_gr, ignore.strand = T), columns = "GENEID") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 12))
  
  #own meth track
  ownmeth_df <- reshape2::melt(data.frame(Group = ownmeth_groups, t(ownmeth_data[ownmeth_ol,])), variable.name = "CpG", value.name = "Methylation")
  ownmeth_df$Pos <- start(ownmeth_anno_gr[ownmeth_df$CpG])
  
  ownmeth_track <- ggplot(ownmeth_df, aes(x = Pos, y = Methylation)) +
    geom_point(aes(col = Group, group = Group), alpha = 0.1) +
    #stat_summary(aes(col = Group, group = Group), fun.y = mean, geom = "smooth") +
    geom_smooth(method = "loess", aes(col = Group), se = F) +
    #ylim(0,1) +
    scale_color_manual(values=c("#00dae0", "#0000ff", "#ff0000")) +
    ylab("% Methylation") +
    theme_bw() +
    theme(legend.pos = "bottom",
          plot.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  #public meth track
  pubmeth_df <- reshape2::melt(data.frame(Group = pubmeth_groups, t(pubmeth_data[pubmeth_ol,])), variable.name = "CpG", value.name = "Methylation")
  pubmeth_df$Pos <- start(pubmeth_anno_gr[as.numeric(gsub("X", "", pubmeth_df$CpG))])+3
  
  pubmeth_track <- ggplot(pubmeth_df, aes(x = Pos, y = Methylation)) +
    geom_point(aes(col = Group, group = Group), alpha = 0.1) +
    geom_smooth(method = "loess", aes(col = Group), se = F) +
    #stat_summary(aes(col = Group, group = Group), fun.y = mean, geom = "smooth") +
    #ylim(0,1) +
    ylab("NML") +
    theme_bw() +
    theme(legend.pos = "bottom",
          plot.title = element_text(face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  #Compiled
  plotobj <- tracks("Gene" = gene_track, 
                    "Own" = ownmeth_track, 
                    "GSE73788" = pubmeth_track,
                    heights = c(3, 5, 5))
}

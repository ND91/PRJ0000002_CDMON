#This script is used to analyze the Crohn 2013 data and find potential DMPs and DMRs
#The package used for data analysis is called minfi and requires all the .idat files to appear in pairs in respective subfolders. 
#This can be achieved using the shell script: idat_finder.sh
require(IlluminaHumanMethylation450kmanifest)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(limma)
require(minfi)
require(missMethyl)
require(shinyMethyl)
require(MethylAid)
require(FlowSorted.Blood.450k)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(BSgenome.Hsapiens.UCSC.hg19)
require(FDb.InfiniumMethylation.hg19)
require(doParallel)
require(pkgmaker)
require(miscTools)
require(rngtools)
require(GenomicRanges)
require(Gviz)
require(rafalib)

sourceDir = '/home/ayliyim/Dropbox/Epimac/Data/Crohn/450k/Monocyte Faltose/Monocytes/data'
set.seed(1)

options(stringsAsFactors = F)
registerDoParallel(detectCores()-4)

# Import the PhenoData which read.450k.exp uses to import the targets. 
# Important: The script will look specifically for a .csv that starts with "Pheno"
targets <- read.450k.sheet(base = sourceDir, pattern = "Pheno")
# For some reason females are considered 2 and males are considered 1. 
# For the blood distribution this needs to be changed to F and M respectively.
targets$gender[targets$gender == 1] <- "F" #Females
targets$gender[targets$gender == 2] <- "M" #Males
targets$Slide <- as.numeric(targets$Slide)

RGset <- read.450k.exp(targets = targets)
phenoDataFrame <- pData(object = RGset)
Mset.raw <- preprocessRaw(rgSet = RGset)

#QC-check: MethylAid
targets.correct <- methylAid_nonvisual(RGset = RGset)
if(nrow(targets.correct) != 0){
  RGset.correct <- read.450k.exp(targets = targets.correct)
  Mset.raw <- preprocessRaw(rgSet = RGset.correct)
} else{
  RGset.correct <- RGset
}

if(methylaid_QC){
  targets.summary <- summarize(targets)
  visualize(targets.summary)
}

#################
# Normalization #
#################
GMset <- preprocessFunnorm(rgSet = RGset.correct)
#GMset <- preprocessSWAN(rgSet = RGset.correct)
#GMset <- mapToGenome(ratioConvert(GMset, what = "both", keepCN = T))

if(visualize_normalization){
  #I think that this way of checking the normalization gives a false sense of security; based on such plots Quantile Normalization would yield
  #the best results, but that need not be true as was described in Irizarry et al. 2014, i.e. a lot of data might be thrown out. 
  
  cols <- c("red", "blue")
  groups.uniq <- unique(pData(Mset.raw)$Cohort)
  m <- match(pData(Mset.raw)$Cohort, groups.uniq)
  
  #Raw
  boxplot(x = getBeta(Mset.raw), col = cols[m], main = "Beta (Raw)", xlab = "Sample Name", ylab = "Beta", xaxt = "n", pch = 20)
  #Remove orignal ticks on the x-axis and replace it with the sample labels
  axis(side = 1, at = seq(1, length(pData(Mset.raw)$Sample_Name), by = 1), labels = FALSE)
  text(x = seq(1, length(pData(Mset.raw)$Sample_Name), by = 1), y = -0.1, labels = pData(Mset.raw)$Sample_Name, srt = 45, pos = 2, xpd = TRUE, cex = 0.6)
  
  #Normalized
  boxplot(x = getBeta(GMset), col = cols[m], main = "Beta (Normalized)", xlab = "Sample Name", ylab = "Beta", xaxt = "n", pch = 20)
  #Remove orignal ticks on the x-axis and replace it with the sample labels
  axis(side = 1, at = seq(1, length(pData(GMset)$Sample_Name), by = 1), labels = FALSE)
  text(x = seq(1, length(pData(GMset)$Sample_Name), by = 1), y = -0.1, labels = pData(GMset)$Sample_Name, srt = 45, pos = 2, xpd = TRUE, cex = 0.6)
}

##########################
# Find annotation probes #
##########################

# Find the 450k genes annotations associated to the CpGs
annotation.total <- getAnnotation(GMset)
annotation.total.gr <- makeGRangesFromDataFrame(annotation.culled, keep.extra.columns = T, start.field = "pos", end.field = "pos")
# Find the transcripts of the genes
gene.transcripts <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

#########################
# Culling of the probes #
#########################

# Drop CpG loci on sex chromosomes
annotation.noSex <- annotation.total[-which(annotation.total$chr %in% c("chrY", "chrX")),]
GMset.noSex <- GMset[rownames(annotation.noSex), ]
# Drop CpG loci with known SNPs # NOT ALL SNPS ARE REMOVED, CHECK WITH UCSC GENOME BROWSER
GMset.noSNPs <- dropLociWithSnps(GMset.noSex, snps = c("SBE","CpG"), maf = 0)
# Drop CpG loci associated to CpG probes that were found to be promiscuous: UNDER CONSTRUCTION
promProbes <- read.csv("~/Dropbox/Epimac/Data/HumanMethylation450k/Non-specific-probes-Illumina450k.csv")[,1]
GMset.nopromProbes <- GMset.noSNPs[featureNames(GMset.noSNPs) %in% rownames(annotation.total[!rownames(annotation.total) %in% promProbes, ]), ]

GMset.culled <- GMset.nopromProbes
annotation.culled <- getAnnotation(GMset.culled)
annotation.culled.gr <- makeGRangesFromDataFrame(annotation.culled, keep.extra.columns = T, start.field = "pos", end.field = "pos")

M.culled <- getM(GMset.culled)
# I have no idea how some samples have negative infinity. It means that the unmethylated channel was 0 (log2(0) == -Inf).
# I will throw these CpGs out.
M.culled <- M.culled[-which(M.culled == -Inf, arr.ind = T, useNames = F)[,1],] 

phenoData.culled <- pData(GMset.culled)

#Annotate factor of interest
fac_int <- gsub("-.*", "", phenoData.culled$inclusionnr)
fac_int <- gsub("[0-9]", "", fac_int)

#######################################################
# Dora Explora! (Exploratory plots for batch effects) #
#######################################################

#Correlation matrix
library(colorRamps)
library(RColorBrewer)
mypar()
M.cor <- cor(M.culled)
cols <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
image(1:ncol(M.culled), 1:ncol(M.culled), M.cor, col = cols, main = "Correlation Matrix", xaxt = "n", xlab = "", yaxt = "n", ylab = "")
axis(2, 1:ncol(M.culled), as.numeric(as.factor(pData(GMset)$Cohort)), las = 2)
axis(1, 1:ncol(M.culled), as.numeric(as.factor(pData(GMset)$Cohort)), las = 2)

#SVD/PCA
M.culled.centered <- M.culled-rowMeans(M.culled)
M.svd <- svd(M.culled.centered)
#Plot shows how much of the variance is explained by the principal components. 
#What we expect to see if no batch effects were present
y0 <- matrix(rnorm(nrow(M.culled)*ncol(M.culled)), nrow(M.culled), ncol(M.culled))
d0 <- svd(y0)$d

LIM <- range(c(d0^2/sum(d0^2), M.svd$d^2/sum(M.svd$d^2)))
plot(d0^2/sum(d0^2), ylab = "Variance Explained", xlab = "Principal Component", pch = 16, ylim = LIM)
#What we actually see
points(M.svd$d^2/sum(M.svd$d^2), ylab = "Variance Explained", xlab = "Principal Component", pch = 16, col = "red")
legend("topright", c("Monocytes2015", "Randomized_Monocytes2015"), col = c("red", "black"), pch = 16)

#Cohort correlation
if(cohort){
  cols <- as.numeric(as.factor(fac_int))
  pchs <- 16
  mypar(1,1)
  plot(M.svd$v[,1], M.svd$v[,2], col = cols, pch = 16, xlab = "PC1", ylab = "PC2")
  legend("bottomleft", legend = unique(fac_int), pch = 16, col = unique(cols))
  
  #Outlier on PC1
  outliers_PC1 <- pData(GMset.culled)[which(M.svd$v[,1] < -.2),]
  outliers_PC2 <- pData(GMset.culled)[which(M.svd$v[,2] < -.35),]
  
  #Boxplot of PCs
  variable <- as.numeric(as.factor(fac_int))
  mypar(5,5)
  for(i in 1:nrow(M.svd$v)){
    boxplot(split(M.svd$v[,i], variable), las = 2, range = 0, main = paste0("PC",i))
    stripchart(split(M.svd$v[,i], variable), add = T, vertical = T, pch = 1, cex = 0.5, col = 1)
    abline(h = 0, cex = 0.5)
  }
  #Correlation of PCs with cohort
  corr <- sapply(1:ncol(M.svd$v), function(i){
    fit <- lm(M.svd$v[,i]~as.factor(fac_int))
    return(summary(fit)$adj.r.squared)
  })
  mypar()
  plot(corr, pch = 21, bg = "red", main = "Correlation Cohort and PCs", ylab = "Correlation", xlab = "PCs")
}

#Age correlation
if(age){
  cols <- as.numeric(round(pData(GMset.culled)$AGE/10))
  pchs <- 16
  mypar(1,1)
  plot(M.svd$v[,1], M.svd$v[,2], col = cols, pch = 16, xlab = "PC1", ylab = "PC2")
  legend("bottomleft", legend = sort(unique(cols*10)), pch = 16, col = sort(unique(cols)))
  #Boxplot of PCs
  variable <- cols
  mypar(5,5)
  for(i in 1:nrow(M.svd$v)){
    boxplot(split(M.svd$v[,i], variable), las = 2, range = 0, main = paste0("PC",i))
    stripchart(split(M.svd$v[,i], variable), add = T, vertical = T, pch = 1, cex = 0.5, col = 1)
    abline(h = 0, cex = 0.5)
  }
  #Correlation of PCs with cohort
  corr <- sapply(1:ncol(M.svd$v), function(i){
    fit <- lm(M.svd$v[,i]~as.factor(round(pData(GMset.culled)$AGE/10)))
    return(summary(fit)$adj.r.squared)
  })
  mypar()
  plot(corr, pch = 21, bg = "red", main = "Correlation Age and PCs", ylab = "Correlation", xlab = "PCs")
}

#Estimate Blood Cell Distribution

if(Bloodcelldistr){
  cellCounts <- estimateCellCounts(rgSet = RGset.correct, compositeCellType = "Blood", cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"), verbose = T, returnAll = F, meanPlot = T)
  CC.df <- data.frame(names = rownames(cellCounts), cohort = fac_int, cellCounts)
  
  CC.df.melted <- melt(CC.df, id.vars = c("names", "cohort"))
  CC.means <- ddply(CC.df.melted, c("cohort", "variable"), summarise,mean=mean(value))
  CC.means.barplot <- qplot(x=cohort, y=mean, fill=variable, data=CC.means, geom="bar", stat="identity", position="dodge")
}


################
# Linear model #
################

if(remove_outlier){
  fac_int <- fac_int[-which(colnames(M.culled) %in% rownames(outliers_PC1))]
  M.culled <- M.culled[,-which(colnames(M.culled) %in% rownames(outliers_PC1))]
  GMset.culled <- GMset.culled[,-which(colnames(GMset.culled) %in% rownames(outliers_PC1))]
} else if(mark_outlier){
  fac_int1 <- fac_int #Backup
  fac_int[which(colnames(M.culled) %in% rownames(outliers_PC1))] <- paste0(fac_int[which(colnames(M.culled) %in% rownames(outliers_PC1))], "_OUT_PC1")
}

#Method 1: Use limma's builtin contrast function
if(method_1){
  #Fit initial model to find negative control features
  mv.group <- as.factor(fac_int)
  mv.design <- model.matrix(~0 + mv.group + pData(GMset.culled)$AGE)
  
  #RUVfit: UNDER CONSTRUCTION NOT EVEN SURE IF POSSIBLE
  lmfit1 <- lmFit(M.culled, mv.design)
  design.cont <- makeContrasts(CDACT-CDREM, CDREM-HC, CDACT-HC, levels = mv.design)
  lmfit1.cont <- contrasts.fit(lmfit1, design.cont)
  lmfit1.ebayes <- eBayes(lmfit1.cont) 
  
  #Choosing the negative controls: I am at a crossroad here, either I can decide to choose negative controls based on whether they are non-significant (P.adj > 0.5) in all groups, or just each group separately.
  lmfit1.top <- topTable(lmfit1.ebayes, coef = 1:3, number = Inf, adjust.method = "BH")
  negcons <- rownames(M.culled) %in% rownames(lmfit1.top[which(lmfit1.top$adj.P.Val > 0.5),])
  
  ruvfit1 <- RUVfit(data = M.culled, design = mv.design, coef = 1, ctl = negcons)
  ruvfit1.adj <- RUVadj(ruvfit)
  ruvfit1.top <- topRUV(ruvfit1.adj)
  
#   lmfit1.topCDACTCDREM <- topTable(lmfit1.ebayes, coef = 1, adjust.method = "BH")
#   lmfit1.topCDREMHC <- topTable(lmfit1.ebayes, coef = 2, adjust.method = "BH")
#   lmfit1.topCDACTHC <- topTable(lmfit1.ebayes, coef = 3, adjust.method = "BH")
#   negcons <- rownames(lmfit1.topCDACTCDREM$adj.P.Value > 0.5 & lmfit1.topCDREMHC$adj.P.Value > 0.5 & lmfit1.topCDACTHC$adj.P.Val > 0.5)
  
  #SVA
  library(sva)
  mv.svafit <- sva(M.culled, mv.design)
  mv.design.sva <- cbind(mv.design, mv.svafit$sv)
  colnames(mv.design.sva) <- c("CDACT", "CDREM", "HC", "Age", "B1", "B2", "B3")
  
  mv.lmfit <- lmFit(M.culled, mv.design.sva)
  
  #mv.lmfit.residuals <- residuals.MArrayLM(object = mv.lmfit, y = M.culled)
  
  #Contrast the different groups:
  # CDACT vs CDREM
  # CDREM vs HC
  # CDACT vs HC
  mv.design.cont <- makeContrasts(CDACT-CDREM, CDREM-HC, CDACT-HC, levels = mv.design.sva)
  mv.lmfit.cont <- contrasts.fit(mv.lmfit, mv.design.cont)
  mv.lmfit.ebayes <- eBayes(mv.lmfit.cont)
  
  mv.lmfit.ebayes.CDACT_CDREM <- topTable(mv.lmfit.ebayes, coef = 1, num = Inf, adjust.method = "BH")
  mv.lmfit.ebayes.CDACT_CDREM <- cbind(mv.lmfit.ebayes.CDACT_CDREM, annotation.total[rownames(mv.lmfit.ebayes.CDACT_CDREM),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  #dmps.CDACT_CDREM <- mv.lmfit.ebayes.CDACT_CDREM[which(mv.lmfit.ebayes.CDACT_CDREM$p.ebayes.BH < 0.05),]
  #dmps.CDACT_CDREM <- cbind(dmps.CDACT_CDREM, annotation.total[rownames(dmps.CDACT_CDREM),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  mv.lmfit.ebayes.CDREM_HC <- topTable(mv.lmfit.ebayes, coef = 2, num = Inf, adjust.method = "BH")
  mv.lmfit.ebayes.CDREM_HC <- cbind(mv.lmfit.ebayes.CDREM_HC, annotation.total[rownames(mv.lmfit.ebayes.CDREM_HC),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  #dmps.CDREM_HC <- mv.lmfit.ebayes.CDREM_HC[which(mv.lmfit.ebayes.CDREM_HC$p.ebayes.BH < 0.05),]
  #dmps.CDREM_HC <- cbind(dmps.CDREM_HC, annotation.total[rownames(dmps.CDREM_HC),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  mv.lmfit.ebayes.CDACT_HC <- topTable(mv.lmfit.ebayes, coef = 3, num = Inf, adjust.method = "BH")
  mv.lmfit.ebayes.CDACT_HC <- cbind(mv.lmfit.ebayes.CDACT_HC, annotation.total[rownames(mv.lmfit.ebayes.CDACT_HC),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  #dmps.CDACT_HC <- mv.lmfit.ebayes.CDACT_HC[which(mv.lmfit.ebayes.CDACT_HC$p.ebayes.BH < 0.05),]
  #dmps.CDACT_HC <- cbind(dmps.CDACT_HC, annotation.total[rownames(dmps.CDACT_HC),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  #CDACT and CDREM vs HC
  mv.group2 <- as.factor(phenoData.culled$disease)
  mv.design2 <- model.matrix(~mv.group2)
  mv.lmfit2 <- lmFit(M.culled, mv.design2)
  mv.lmfit2.ebayes <- eBayes(mv.lmfit2)
  mv.lmfit2.ebayes.top <- topTable(mv.lmfit2.ebayes, coef = 2, num = Inf, sort.by = "P")
}

#Method 2: split M.culled into sub matrices and calculate the linear model for each submatrix. NOT RECOMMENDED
if(method_2){
  
  #CDACT vs CDREM
  GMset.CDACT_CDREM <- GMset.culled[,which(fac_int == "CDACT" | fac_int == "CDREM")] 
  M.CDACT_CDREM <- M.culled[,which(fac_int == "CDACT" | fac_int == "CDREM")] 
  cont.CDACT_CDREM <- fac_int[which(fac_int == "CDACT" | fac_int == "CDREM")]
  design.CDACTCDREM <- model.matrix(~as.factor(cont.CDACT_CDREM))
  lmfit.CDACTCDREM <- lmFit(M.CDACT_CDREM, design.CDACTCDREM)
  lmfit.CDACTCDREM.ebayes <- eBayes(lmfit.CDACTCDREM)
  lmfit.CDACTCDREM.ebayes.top <- topTable(lmfit.CDACTCDREM.ebayes, coef = 2, num = Inf, sort.by = "P")
  lmfit.CDACTCDREM.ebayes.annot <- cbind(lmfit.CDACTCDREM.ebayes.top, annotation.total[rownames(lmfit.CDACTCDREM.ebayes.top),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  dmps.CDACTCDREM <- lmfit.CDACTCDREM.ebayes.top[which(lmfit.CDACTCDREM.ebayes.top$p.ebayes.BH < 0.05),]
  #SVA
  if(SVA){
    svafit.CDACTCDREM <- sva(M.CDACT_CDREM, design.CDACTCDREM)
    
    sva.CDACT_CDREM.design <- model.matrix(~as.factor(cont.CDACT_CDREM)+svafit.CDACTCDREM$sv)
    sva.CDACT_CDREM.lmfit <- lmFit(M.CDACT_CDREM, sva.CDACT_CDREM.design)
    sva.CDACT_CDREM.tt <- sva.CDACT_CDREM.lmfit$coef[,2]*sqrt(sva.CDACT_CDREM.lmfit$df.residual)/(2*sva.CDACT_CDREM.lmfit$sigma)
    sva.CDACT_CDREM.res <- data.frame(dm = -sva.CDACT_CDREM.lmfit$coef[,2], p.value = 2*(1-pt(abs(sva.CDACT_CDREM.tt), sva.CDACT_CDREM.lmfit$df.residual[1])))
    
    mypar(1,2)
    hist(sva.CDACT_CDREM.res$p.value, breaks = 100, main = "SVA", xlab = "P-values")
    plot(sva.CDACT_CDREM.res$dm, -log10(sva.CDACT_CDREM.res$p.value))
    
    sva.CDACT_CDREM.batch <- sva.CDACT_CDREM.lmfit$coef[,-c(1,2)]%*%t(sva.CDACT_CDREM.design[,-c(1,2)])
    sva.CDACT_CDREM.signal <- sva.CDACT_CDREM.lmfit$coef[,c(1,2)]%*%t(sva.CDACT_CDREM.design[,c(1,2)])
    sva.CDACT_CDREM.error <- M.CDACT_CDREM-sva.CDACT_CDREM.signal-sva.CDACT_CDREM.batch
    
    sva.CDACT_CDREM.ebayes <- eBayes(sva.CDACT_CDREM.lmfit)
    sva.CDACT_CDREM.top <- topTable(fit = sva.CDACT_CDREM.ebayes, coef = 2, number = nrow(sva.CDACT_CDREM.ebayes), adjust.method = "BH")
    sva.CDACT_CDREM.top <- cbind(sva.CDACT_CDREM.top, annotation.total[rownames(sva.CDACT_CDREM.top),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
    
    #Plotting the decomposition
    #Demean for plot
    M.CDACT_CDREM.demean <- M.CDACT_CDREM-rowMeans(M.CDACT_CDREM)
    sva.CDACT_CDREM.signal.demean <- sva.CDACT_CDREM.signal - rowMeans(sva.CDACT_CDREM.signal)
    
    mypar(1,4, mar = c(2.75, 4.5, 2.6, 1.1))
    image(t(M.CDACT_CDREM.demean[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Total")
    image(t(sva.CDACT_CDREM.signal.demean[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Case vs Control")
    image(t(sva.CDACT_CDREM.batch[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Batches")
    image(t(sva.CDACT_CDREM.error[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Unaccounted")
    
    #DMR
    sva.CDACT_CDREM.dmr <- bumphunter(object = GMset.CDACT_CDREM, design = sva.CDACT_CDREM.design, coef = 2, cutoff = 0.08, B = 100, nullMethod = "bootstrap", smooth = T, smoothFunction = loessByCluster)
    sva.CDACT_CDREM.dmr.culled <- sva.CDACT_CDREM.dmr$table[sva.CDACT_CDREM.dmr$table$L >=2,]
    sva.CDACT_CDREM.dmr.culled.gr <- makeGRangesFromDataFrame(sva.CDACT_CDREM.dmr.culled, keep.extra.columns = TRUE)
    #Find nearest genes and append information
    sva.CDACT_CDREM.nearestgenes.gr <- gene.transcripts[nearest(x = sva.CDACT_CDREM.dmr.culled.gr, subject = gene.transcripts),]
    sva.CDACT_CDREM.dmr.culled <- cbind(sva.CDACT_CDREM.dmr.culled, Gene_start = as.data.frame(ranges(sva.CDACT_CDREM.nearestgenes.gr))$start, Gene_end = as.data.frame(ranges(sva.CDACT_CDREM.nearestgenes.gr))$end, Gene_Entrez = sva.CDACT_CDREM.nearestgenes.gr$Entrez, Gene_Name = sva.CDACT_CDREM.nearestgenes.gr$Gene)
  }
  
  #CDREM vs HC
  GMset.CDREM_HC <- GMset.culled[,which(fac_int == "CDREM" | fac_int == "HC")] 
  M.CDREM_HC <- M.culled[,which(fac_int == "CDREM" | fac_int == "HC")]
  cont.CDREM_HC <- fac_int[which(fac_int == "CDREM" | fac_int == "HC")]
  design.CDREMHC <- model.matrix(~as.factor(cont.CDREM_HC))
  lmfit.CDREMHC <- lmFit(M.CDREM_HC, design.CDREMHC)
  lmfit.CDREMHC.ebayes <- eBayes(lmfit.CDREMHC)
  lmfit.CDREMHC.ebayes.top <- topTable(lmfit.CDREMHC.ebayes, coef = 2, num = Inf, sort.by = "P")
  lmfit.CDREMHC.ebayes.annot <- cbind(lmfit.CDREMHC.ebayes.top, annotation.total[rownames(lmfit.CDREMHC.ebayes.top),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  dmps.CDREMHC <- lmfit.CDREMHC.ebayes.top[which(lmfit.CDREMHC.ebayes.top$p.ebayes.BH < 0.05),]
  
  #SVA
  if(SVA){
    svafit.CDREM_HC <- sva(M.CDREM_HC, design.CDREMHC)
    
    sva.CDREM_HC.design <- model.matrix(~as.factor(cont.CDREM_HC)+svafit.CDREM_HC$sv)
    sva.CDREM_HC.lmfit <- lmFit(M.CDREM_HC, sva.CDREM_HC.design)
    sva.CDREM_HC.tt <- sva.CDREM_HC.lmfit$coef[,2]*sqrt(sva.CDREM_HC.lmfit$df.residual)/(2*sva.CDREM_HC.lmfit$sigma)
    sva.CDREM_HC.res <- data.frame(dm = -sva.CDREM_HC.lmfit$coef[,2], p.value = 2*(1-pt(abs(sva.CDREM_HC.tt), sva.CDREM_HC.lmfit$df.residual[1])))
    
    mypar(1,2)
    hist(sva.CDREM_HC.res$p.value, breaks = 100, main = "P-values", xlab = "P-values")
    plot(sva.CDREM_HC.res$dm, -log10(sva.CDREM_HC.res$p.value), main = "Volcano")
    
    sva.CDREM_HC.batch <- sva.CDREM_HC.lmfit$coef[,-c(1,2)]%*%t(sva.CDREM_HC.design[,-c(1,2)])
    sva.CDREM_HC.signal <- sva.CDREM_HC.lmfit$coef[,c(1,2)]%*%t(sva.CDREM_HC.design[,c(1,2)])
    sva.CDREM_HC.error <- M.CDREM_HC-sva.CDREM_HC.signal-sva.CDREM_HC.batch
    
    sva.CDREM_HC.ebayes <- eBayes(sva.CDREM_HC.lmfit)
    sva.CDREM_HC.top <- topTable(fit = sva.CDREM_HC.ebayes, coef = 2, number = nrow(sva.CDREM_HC.ebayes), adjust.method = "BH")
    sva.CDREM_HC.top <- cbind(sva.CDREM_HC.top, annotation.total[rownames(sva.CDREM_HC.top),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
    
    #Plotting the decomposition
    #Demean for plot
    M.CDREM_HC.demean <- M.CDREM_HC-rowMeans(M.CDREM_HC)
    sva.CDREM_HC.signal.demean <- sva.CDREM_HC.signal - rowMeans(sva.CDREM_HC.signal)
    
    mypar(1,4, mar = c(2.75, 4.5, 2.6, 1.1))
    image(t(M.CDREM_HC.demean[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Total")
    image(t(sva.CDREM_HC.signal.demean[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Case vs Control")
    image(t(sva.CDREM_HC.batch[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Batches")
    image(t(sva.CDREM_HC.error[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Unaccounted")
    
    #DMR
    sva.CDREM_HC.dmr <- bumphunter(object = GMset.CDREM_HC, design = sva.CDREM_HC.design, coef = 2, cutoff = 0.08, B = 100, nullMethod = "bootstrap", smooth = T, smoothFunction = loessByCluster)
    sva.CDREM_HC.dmr.culled <- sva.CDREM_HC.dmr$table[sva.CDREM_HC.dmr$table$L >=2,]
    sva.CDREM_HC.dmr.culled.gr <- makeGRangesFromDataFrame(sva.CDREM_HC.dmr.culled, keep.extra.columns = TRUE)
    #Find nearest genes and append information
    sva.CDREM_HC.nearestgenes.gr <- gene.transcripts[nearest(x = sva.CDREM_HC.dmr.culled.gr, subject = gene.transcripts),]
    sva.CDREM_HC.dmr.culled <- cbind(sva.CDREM_HC.dmr.culled, Gene_start = as.data.frame(ranges(sva.CDREM_HC.nearestgenes.gr))$start, Gene_end = as.data.frame(ranges(sva.CDREM_HC.nearestgenes.gr))$end, Gene_Entrez = sva.CDREM_HC.nearestgenes.gr$Entrez, Gene_Name = sva.CDREM_HC.nearestgenes.gr$Gene)
    
  }
  
  #CDACT vs HC
  M.CDACT_HC <- M.culled[,which(fac_int == "CDACT" | fac_int == "HC")]
  cont.CDACT_HC <- fac_int[which(fac_int == "CDACT" | fac_int == "HC")]
  design.CDACTHC <- model.matrix(~as.factor(cont.CDACT_HC))
  lmfit.CDACTHC <- lmFit(M.CDACT_HC, design.CDACTHC)
  lmfit.CDACTHC.ebayes <- eBayes(lmfit.CDACTHC)
  lmfit.CDACTHC.ebayes.top <- topTable(lmfit.CDACTHC.ebayes, coef = 2, num = Inf, sort.by = "P")
  lmfit.CDACTHC.ebayes.annot <- cbind(lmfit.CDACTHC.ebayes.top, annotation.total[rownames(lmfit.CDACTHC.ebayes.top),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  dmps.CDACTHC <- lmfit.CDACTHC.ebayes.top[which(lmfit.CDACTHC.ebayes.top$p.ebayes.BH < 0.05),]
  #SVA
  if(SVA){
    svafit.CDACT_HC <- sva(M.CDACT_HC, design.CDACTHC)
    
    sva.CDACT_HC.design <- model.matrix(~as.factor(cont.CDACT_HC)+svafit.CDACT_HC$sv)
    sva.CDACT_HC.lmfit <- lmFit(M.CDACT_HC, sva.CDACT_HC.design)
    sva.CDACT_HC.tt <- sva.CDACT_HC.lmfit$coef[,2]*sqrt(sva.CDACT_HC.lmfit$df.residual)/(2*sva.CDACT_HC.lmfit$sigma)
    sva.CDACT_HC.res <- data.frame(dm = -sva.CDACT_HC.lmfit$coef[,2], p.value = 2*(1-pt(abs(sva.CDACT_HC.tt), sva.CDACT_HC.lmfit$df.residual[1])))
    
    mypar(1,2)
    hist(sva.CDACT_HC.res$p.value, breaks = 100, main = "P-values", xlab = "P-values")
    plot(sva.CDACT_HC.res$dm, -log10(sva.CDACT_HC.res$p.value), main = "Volcano")
    
    sva.CDACT_HC.batch <- sva.CDACT_HC.lmfit$coef[,-c(1,2)]%*%t(sva.CDACT_HC.design[,-c(1,2)])
    sva.CDACT_HC.signal <- sva.CDACT_HC.lmfit$coef[,c(1,2)]%*%t(sva.CDACT_HC.design[,c(1,2)])
    sva.CDACT_HC.error <- M.CDACT_HC-sva.CDACT_HC.signal-sva.CDACT_HC.batch
    
    sva.CDACT_HC.ebayes <- eBayes(sva.CDACT_HC.lmfit)
    sva.CDACT_HC.top <- topTable(fit = sva.CDACT_HC.ebayes, coef = 2, number = nrow(sva.CDACT_HC.ebayes), adjust.method = "BH")
    sva.CDACT_HC.top <- cbind(sva.CDACT_HC.top, annotation.total[rownames(sva.CDACT_HC.top),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
    
    #Plotting the decomposition
    #Demean for plot
    M.CDACT_HC.demean <- M.CDACT_HC-rowMeans(M.CDACT_HC)
    sva.CDACT_HC.signal.demean <- sva.CDACT_HC.signal - rowMeans(sva.CDACT_HC.signal)
    
    mypar(1,4, mar = c(2.75, 4.5, 2.6, 1.1))
    image(t(M.CDACT_HC.demean[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Total")
    image(t(sva.CDACT_HC.signal.demean[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Case vs Control")
    image(t(sva.CDACT_HC.batch[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Batches")
    image(t(sva.CDACT_HC.error[1:1000,]), col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100), zlim = c(-5,5), xaxt = "n", yaxt = "n", main = "Unaccounted")
  }
}

#Method 3: split M.culled according to cases and controls (no three levels anymore; only Crohn's or Control)
if(method_3){
  
  pheno.3 <- gsub("(CDACT|CDREM)", "Crohn", fac_int)
  bin.design <- model.matrix(~relevel(factor(pheno.3), "HC"))
  
  #No correction
  lmfit.3 <- lmFit(object = M.culled, design = bin.design)
  lmfit.3.ebayes <- eBayes(lmfit.3)
  top.3 <- topTable(fit = lmfit.3.ebayes, coef = 2, number = Inf, adjust.method = "BH")
  top.3 <- cbind(top.3, annotation.total[rownames(top.3),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  #RUVfit
  #negcons.3 <- !(rownames(top.3) %in% rownames(top.3[nrow(top.3)-1000:nrow(top.3), ]))
  negcons.3 <- rownames(M.culled) %in% rownames(top.3[which(top.3$adj.P.Val > 0.5),])
  ruvfit.3 <- RUVfit(data = M.culled, design = bin.design, coef = 2, ctl = negcons.3)
  ruvfit.3.adj <- RUVadj(ruvfit.3)
  ruvfit.3.top <- topRUV(ruvfit.3.adj, number = Inf)
  ruvfit.3.top <- cbind(ruvfit.3.top, annotation.total[rownames(ruvfit.3.top),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  #SVA
  svafit.3 <- sva(M.culled, bin.design)
  design.sva.3 <- cbind(bin.design, svafit.3$sv)
  colnames(design.sva.3) <- c("(Intercept)", "Crohn", "1", "2", "3", "4")
  
  lmfit.sva.3 <- lmFit(M.culled, design.sva.3)
  lmfit.sva.ebayes.3 <- eBayes(lmfit.sva.3)
  lmfit.sva.top.3 <- topTable(fit = lmfit.sva.ebayes.3, coef = 2, number = Inf, adjust.method = "BH")
  lmfit.sva.top.3 <- cbind(lmfit.sva.top.3, annotation.total[rownames(lmfit.sva.top.3),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  #GO
  require(BiasedUrn)
  require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  go.sva <- gometh(sig.cpg = rownames(lmfit.sva.top.3[which(lmfit.sva.top.3$adj.P.Val < 0.05),]), all.cpg = rownames(annotation.culled), plot.bias = T, prior.prob = T)
  go.sva.sig <- go.sva[go.sva$FDR < 0.05,]
  
  #Comparison
  hist(top.3$P.Value, breaks = 1000, main = "None")
  qqnorm(top.3$t, pch = 16, main = "None", ylim = c(-10,10))
  abline(0,1)
  hist(lmfit.sva.top.3$P.Value, breaks = 1000, main = "SVA")
  qqnorm(lmfit.sva.top.3$t, pch = 16, main = "SVA", ylim = c(-10,10))
  abline(0,1)
  hist(ruvfit.3.top$p, breaks = 1000, main = "RUV")
  qqnorm(ruvfit.3.top$t, pch = 16, main = "RUV", ylim = c(-10,10))
  abline(0,1)
  
}

#####################
# Hypothesis driven #
#####################

GWAS_genes <- read.csv('/home/ayliyim/Dropbox/Epimac/Data/Crohn Whole Blood 2013/Target_Genes/Candidate Genes V3_GWAS.csv')
EWAS_genes <- read.csv('/home/ayliyim/Dropbox/Epimac/Data/Crohn Whole Blood 2013/Target_Genes/Candidate Genes V3_EWAS.csv')

ewasdata <- unlist(unique(EWAS_genes$GENE))
ewasdata <- ewasdata[-which(ewasdata == "")]
gwasdata <- unlist(unique(GWAS_genes$GENE))
gwasdata <- gwasdata[-which(gwasdata == "")]

#Hypothesis-driven: Genes obtained from EWAS and GWAS
unique_genes <- c(as.vector(ewasdata), as.vector(gwasdata))

#Using the annotated gene data
cpg <- sapply(unique_genes, function(x) names(annotation.culled.gr[grep(paste0("(^|;)*", x, "($|;)*"), annotation.culled.gr$UCSC_RefGene_Name), ]))
cpg.vector <- unique(as.vector(unlist(cpg)))

require(dplyr)

#None
top3.hypdriv <- top.3[cpg.vector,]
top3.hypdriv$p.BH.hypdriv <- p.adjust(top3.hypdriv$P.Value, method = "BH")
top3.hypdriv <- top3.hypdriv[order(top3.hypdriv$p.BH.hypdriv),]
#SVA
svafit.3.hypdriv <- lmfit.sva.top.3[cpg.vector,]
svafit.3.hypdriv$p.BH.hypdriv <- p.adjust(svafit.3.hypdriv$P.Value, method = "BH")
svafit.3.hypdriv <- svafit.3.hypdriv[order(svafit.3.hypdriv$p.BH.hypdriv),]
#RUV
ruvfit.3.hypdriv <- ruvfit.3.top[cpg.vector,]
ruvfit.3.hypdriv$p.BH.hypdriv <- p.adjust(ruvfit.3.hypdriv$p, method = "BH")
ruvfit.3.hypdriv <- ruvfit.3.hypdriv[order(ruvfit.3.hypdriv$p.BH.hypdriv),]

##############################################
# RUV: See missMethyl package (multivariate) #
##############################################

#Negative controls 
neg_CDACTCDREM <- rownames(mv.lmfit.ebayes.CDACT_CDREM[mv.lmfit.ebayes.CDACT_CDREM$adj.P.Val > 0.5,])
neg_CDREMHC <- rownames(mv.lmfit.ebayes.CDREM_HC[mv.lmfit.ebayes.CDREM_HC$adj.P.Val > 0.5,])
neg_CDACTHC <- rownames(mv.lmfit.ebayes.CDACT_HC[mv.lmfit.ebayes.CDACT_HC$adj.P.Val > 0.5,])

#CDACT vs CDREM
#Perform RUV adjustment and fit
if(CDACT_v_CDREM){
  negcons <- rownames(M.culled) %in% neg_CDACTCDREM
  
  ruvfit <- RUVfit(data = M.culled, design = mv.design, coef = 1, ctl = negcons)
  cat("Found", dim(ruvfit$W)[1], "unwanted variables")
  ruvfit.ebayes <- RUVadj(ruvfit)
  
  LIM <- range(-log10(mv.lmfit.ebayes$p.value[,1]), -log10(ruvfit.ebayes$p.ebayes))
  #Uncorrected
  #Volcano
  plot(mv.lmfit.ebayes$coef[,1], -log10(mv.lmfit.ebayes$p.value[,1]), cex = 0.8, pch = 21, bg = "red", main = "Uncorrected Volcano Plot", xlab = "Coefficient", ylab = "-log10(Pval)", xlim = c(-3, 3), ylim = LIM)
  abline(h = -log10(0.05))
  #Check P-value histogram prior to correction
  hist(mv.lmfit.ebayes$p.val[,2], breaks = 100, main = "P-vals")
  #Negative controls
  hist(mv.lmfit.ebayes$p.val[negcons,2], breaks = 100, main = "P-val negative controls")
  #Corrected
  #Volcano (P-values)
  plot(ruvfit.ebayes$coef, -log10(ruvfit.ebayes$p.ebayes), cex = 0.8, pch = 21, bg = "red", main = "RUV Volcano Plot", xlab = "Coefficient", ylab = "-log10(Pval)", ylim = LIM)
  abline(h = -log10(0.05), lty = 1)
  #Volcano (P-adjusted)
  LIM <- range(-log10(ruvfit.ebayes$p.ebayes.BH))
  plot(ruvfit.ebayes$coef, -log10(ruvfit.ebayes$p.ebayes.BH), cex = 0.8, pch = 21, bg = "red", main = "RUV Volcano Plot", xlab = "Coefficient", ylab = "-log10(P.adj)", xlim = c(-3, 3), ylim = LIM)
  abline(h = -log10(0.05), lty = 1)
  
  #Check P-value histogram after correction
  hist(ruvfit.ebayes$p.ebayes, breaks = 100, main = "P-vals")
  #Negative controls
  hist(ruvfit.ebayes$p.ebayes[negcons], breaks = 100, main = "P-vals negative controls")
  
  #Find DMPs
  ruv.dmps <- topRUV(ruvfit.ebayes, number = Inf)
  ruv.dmps.significant <- ruv.dmps[which(ruv.dmps$p.ebayes.BH < 0.05),]
  ruv.dmps.significant <- cbind(ruv.dmps.significant, annotation.total[rownames(ruv.dmps.significant),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  #T-statistic
  hist(ruvfit.ebayes$t, breaks = 1000, main = "RUV Correction T-distribution")
  qqnorm(ruvfit.ebayes$t)
  qqline(ruvfit.ebayes$t)
  #Lambda
  ruv.chisq <- qchisq(p = 1-as.vector(ruvfit.ebayes$p), df = 1)
  ruv.lambda <- median(ruv.chisq)/qchisq(0.5,1)
  #Kolmogorov-Smirnov
  ks.test(ruvfit.ebayes$p, "punif", 0, 1)
  
  #Find DMRs
  ruv.design <- cbind(mv.design, t(ruvfit$W))
  
  #Smoothing = F
  #ruv.dmrs <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, B = 100, nullMethod = "bootstrap", what = "M")
  ruv.dmrs <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, nullMethod = "bootstrap", what = "M")
  if(!is.null(ruv.dmrs$table)){
    ruv.dmrs.culled <- ruv.dmrs$tab[ruv.dmrs$table$L >=2,]
    #Find nearest genes
    ruv.dmrs.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.culled, keep.extra.columns = TRUE)
    ruv.nearestgenes.granges <- gene.transcripts[nearest(x = ruv.dmrs.culled.gr, subject = gene.transcripts),]
    #Append ranges, Entrez and Gene to GRanges object
    ruv.dmrs.culled <- cbind(ruv.dmrs.culled, Gene_start = as.data.frame(ranges(ruv.nearestgenes.granges))$start, Gene_end = as.data.frame(ranges(ruv.nearestgenes.granges))$end, Gene_Entrez = ruv.nearestgenes.granges$Entrez, Gene_Name = ruv.nearestgenes.granges$Gene)
    ruv.dmrs.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.culled, keep.extra.columns = TRUE)
  }
  
  #Smoothing = T
  #ruv.dmrs.smooth <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, B = 100, nullMethod = "bootstrap", smooth = T, smoothFunction = loessByCluster)
  ruv.dmrs.smooth <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, nullMethod = "bootstrap", smooth = T, smoothFunction = loessByCluster)
  if(!is.null(ruv.dmrs.smooth$table)){
    ruv.dmrs.smooth.culled <- ruv.dmrs.smooth$tab[ruv.dmrs.smooth$table$L >=2,]
    #Find nearest genes
    ruv.dmrs.smooth.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.smooth.culled, keep.extra.columns = TRUE)
    ruv.nearestgenes.granges <- gene.transcripts[nearest(x = ruv.dmrs.smooth.culled.gr, subject = gene.transcripts),]
    #Append ranges, Entrez and Gene to GRanges object
    ruv.dmrs.smooth.culled <- cbind(ruv.dmrs.smooth.culled, Gene_start = as.data.frame(ranges(ruv.nearestgenes.granges))$start, Gene_end = as.data.frame(ranges(ruv.nearestgenes.granges))$end, Gene_Entrez = ruv.nearestgenes.granges$Entrez, Gene_Name = ruv.nearestgenes.granges$Gene)
    ruv.dmrs.smooth.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.smooth.culled, keep.extra.columns = TRUE)
  }
}
if(CDACT_v_HC){
  negcons <- rownames(M.culled) %in% neg_CDACTHC
  
  ruvfit <- RUVfit(data = M.culled, design = mv.design, coef = 3, ctl = negcons)
  cat("Found", dim(ruvfit$W)[1], "unwanted variables")
  ruvfit.ebayes <- RUVadj(ruvfit)
  
  LIM <- range(-log10(mv.lmfit.ebayes$p.value[,1]), -log10(ruvfit.ebayes$p.ebayes))
  #Uncorrected
  #Volcano
  plot(mv.lmfit.ebayes$coef[,1], -log10(mv.lmfit.ebayes$p.value[,1]), cex = 0.8, pch = 21, bg = "red", main = "Uncorrected Volcano Plot", xlab = "Coefficient", ylab = "-log10(Pval)", xlim = c(-3, 3), ylim = LIM)
  abline(h = -log10(0.05))
  #Check P-value histogram prior to correction
  hist(mv.lmfit.ebayes$p.val[,2], breaks = 100, main = "P-vals")
  #Negative controls
  hist(mv.lmfit.ebayes$p.val[negcons,2], breaks = 100, main = "P-val negative controls")
  #Corrected
  #Volcano (P-values)
  plot(ruvfit.ebayes$coef, -log10(ruvfit.ebayes$p.ebayes), cex = 0.8, pch = 21, bg = "red", main = "RUV Volcano Plot", xlab = "Coefficient", ylab = "-log10(Pval)", ylim = LIM)
  abline(h = -log10(0.05), lty = 1)
  #Volcano (P-adjusted)
  LIM <- range(-log10(ruvfit.ebayes$p.ebayes.BH))
  plot(ruvfit.ebayes$coef, -log10(ruvfit.ebayes$p.ebayes.BH), cex = 0.8, pch = 21, bg = "red", main = "RUV Volcano Plot", xlab = "Coefficient", ylab = "-log10(P.adj)", xlim = c(-3, 3), ylim = LIM)
  abline(h = -log10(0.05), lty = 1)
  
  #Check P-value histogram after correction
  hist(ruvfit.ebayes$p.ebayes, breaks = 100, main = "P-vals")
  #Negative controls
  hist(ruvfit.ebayes$p.ebayes[negcons], breaks = 100, main = "P-vals negative controls")
  
  #Find DMPs
  ruv.dmps <- topRUV(ruvfit.ebayes, number = Inf)
  ruv.dmps.significant <- ruv.dmps[which(ruv.dmps$p.ebayes.BH < 0.05),]
  ruv.dmps.significant <- cbind(ruv.dmps.significant, annotation.total[rownames(ruv.dmps.significant),c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
  
  #T-statistic
  hist(ruvfit.ebayes$t, breaks = 1000, main = "RUV Correction T-distribution")
  qqnorm(ruvfit.ebayes$t)
  qqline(ruvfit.ebayes$t)
  #Lambda
  ruv.chisq <- qchisq(p = 1-as.vector(ruvfit.ebayes$p), df = 1)
  ruv.lambda <- median(ruv.chisq)/qchisq(0.5,1)
  #Kolmogorov-Smirnov
  ks.test(ruvfit.ebayes$p, "punif", 0, 1)
  
  #Find DMRs
  ruv.design <- cbind(mv.design, t(ruvfit$W))
  
  #Smoothing = F
  #ruv.dmrs <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, B = 100, nullMethod = "bootstrap", what = "M")
  ruv.dmrs <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, nullMethod = "bootstrap", what = "M")
  if(!is.null(ruv.dmrs$table)){
    ruv.dmrs.culled <- ruv.dmrs$tab[ruv.dmrs$table$L >=2,]
    #Find nearest genes
    ruv.dmrs.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.culled, keep.extra.columns = TRUE)
    ruv.nearestgenes.granges <- gene.transcripts[nearest(x = ruv.dmrs.culled.gr, subject = gene.transcripts),]
    #Append ranges, Entrez and Gene to GRanges object
    ruv.dmrs.culled <- cbind(ruv.dmrs.culled, Gene_start = as.data.frame(ranges(ruv.nearestgenes.granges))$start, Gene_end = as.data.frame(ranges(ruv.nearestgenes.granges))$end, Gene_Entrez = ruv.nearestgenes.granges$Entrez, Gene_Name = ruv.nearestgenes.granges$Gene)
    ruv.dmrs.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.culled, keep.extra.columns = TRUE)
  }
  
  #Smoothing = T
  #ruv.dmrs.smooth <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, B = 100, nullMethod = "bootstrap", smooth = T, smoothFunction = loessByCluster)
  ruv.dmrs.smooth <- bumphunter(object = GMset.culled, design = ruv.design, coef = 2, cutoff = 0.08, nullMethod = "bootstrap", smooth = T, smoothFunction = loessByCluster)
  if(!is.null(ruv.dmrs.smooth$table)){
    ruv.dmrs.smooth.culled <- ruv.dmrs.smooth$tab[ruv.dmrs.smooth$table$L >=2,]
    #Find nearest genes
    ruv.dmrs.smooth.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.smooth.culled, keep.extra.columns = TRUE)
    ruv.nearestgenes.granges <- gene.transcripts[nearest(x = ruv.dmrs.smooth.culled.gr, subject = gene.transcripts),]
    #Append ranges, Entrez and Gene to GRanges object
    ruv.dmrs.smooth.culled <- cbind(ruv.dmrs.smooth.culled, Gene_start = as.data.frame(ranges(ruv.nearestgenes.granges))$start, Gene_end = as.data.frame(ranges(ruv.nearestgenes.granges))$end, Gene_Entrez = ruv.nearestgenes.granges$Entrez, Gene_Name = ruv.nearestgenes.granges$Gene)
    ruv.dmrs.smooth.culled.gr <- makeGRangesFromDataFrame(ruv.dmrs.smooth.culled, keep.extra.columns = TRUE)
  }
}

#######################
# Visualization DMRS #
#######################

#DMR
write.csv(as.data.frame(annotation.culled.gr[which(start(annotation.culled.gr) >= start_dmr & start(annotation.culled.gr) <= end_dmr & seqnames(annotation.culled.gr) == chr),]), "LINC0062_DMR_DMPs.csv")

#Where are the DMRs are located?
for(i in 1:10){
  ind_dmps <- dmr_dmps(ruv.dmrs.culled.gr[i,], annotation.culled.gr)
}

length(ruv.dmrs.culled.gr)
#Plot RUV DMR
dir.create("DMR_data")
for(i in 1:dim(ruv.dmrs.culled)[1]){
  ruv_dmrs <- dmr_plotter(name = paste(ruv.dmrs.culled[i,]$Gene_Name), 
                          chr = paste(ruv.dmrs.culled[i,]$chr), 
                          gmset = GMset.culled, 
                          annotation.gr = annotation.culled.gr,
                          genome_version = "hg19",
                          start_dmr = ruv.dmrs.culled[i,]$start,
                          end_dmr = ruv.dmrs.culled[i,]$end,
                          flanks = 1000)
  ruv_dmrs_df <- as.data.frame(ruv_dmrs)
  ruv_dmrs_df <- cbind(ruv_dmrs_df, ruv.dmps[names(ruv_dmrs),])
  write.csv(ruv_dmrs_df, paste0(getwd(), "/DMR_data/", ruv.dmrs.culled[i,]$Gene_Name, ".csv"))
}

##########################################
# Hypothesis-driven approach: GWAS check #
##########################################

#GWAS genes obtained from the articles by Franke et al. 2012 and Jostins et al. 2012
GWAS_genes <- read.csv('/home/ayliyim/Dropbox/Epimac/Data/Crohn/Target_Genes/Candidate_genes_v2.csv')
unique_genes <- unique(as.vector(t(GWAS_genes)))

#Add single names to the annotation.culled.gr
annotation.culled.gr$single_UCSC_Name <- gsub("(.*?);.*", "\\1", annotation.culled.gr$UCSC_RefGene_Name)

#Using the annotated gene data
#cpg <- sapply(unique_genes, function(x) names(annotation.culled.gr[grep(paste0(x, ), annotation.culled.gr$UCSC_RefGene_Name), ]))
cpg <- sapply(unique_genes, function(x) names(annotation.culled.gr[grep(paste0("^", x, "$"), annotation.culled.gr$single_UCSC_Name), ]))

cpg.vector <- unique(as.vector(unlist(cpg)))
selected_DMPs <- data.frame(seqnames = seqnames(annotation.culled.gr[cpg.vector, ]),
                            pos = start(annotation.culled.gr[cpg.vector, ]),
                            strand = strand(annotation.culled.gr[cpg.vector, ]),
                            UCSC_gene = annotation.culled.gr[cpg.vector, ]$single_UCSC_Name, 
                            #UCSC_Accession = annotation.culled.gr[cpg.vector, ]$UCSC_RefGene_Accession, 
                            UCSC_gene_all = annotation.culled.gr[cpg.vector, ]$UCSC_RefGene_Name, 
                            UCSC_group = annotation.culled.gr[cpg.vector, ]$UCSC_RefGene_Group, 
                            Regulatory_Feature = annotation.culled.gr[cpg.vector, ]$Regulatory_Feature_Group,
                            coefficient = as.vector(ruvfit.ebayes$coef[cpg.vector,]),
                            p.ebayes = ruv.dmps[cpg.vector, "p.ebayes"], 
                            p.ebayes.BH_hypfree = ruv.dmps[cpg.vector, "p.ebayes.BH"],
                            p.ebayes.BH_hypdriv = p.adjust(ruv.dmps[cpg.vector, "p.ebayes"], method = "BH"))
rownames(selected_DMPs) <- cpg.vector
#nofilter
selected_DMPs <- selected_DMPs[order(selected_DMPs$p.ebayes.BH_hypdriv, decreasing = F),]

#Missing genes in list
missing_genes <- unique_genes[which(!unique_genes %in% unique(selected_DMPs$UCSC_gene))]
annotation.culled.gr[which(annotation.culled.gr$UCSC_RefGene_Name == missing_genes),]

#Volcano plot with colored chosen CpGs
indices_chosen <- rownames(ruvfit.ebayes$coef) %in% cpg.vector
plot(ruvfit.ebayes$coef[!indices_chosen], -log10(ruvfit.ebayes$p.ebayes[!indices_chosen,]), cex = 0.8, pch = 21, bg = "red", main = "RUV Volcano Plot", xlab = "Effect size", ylab = "-log10(Pval)", xlim = c(-3,3))
points(ruvfit.ebayes$coef[indices_chosen,], -log10(ruvfit.ebayes$p.ebayes[indices_chosen,]), cex = 0.8, pch = 21, bg = "yellow")
abline(h = -log10(0.05), lty = 1)
legend("bottomright", "Selected CpGs", col = "yellow", pch = 16)

#filter
selected_DMPs.sig <- selected_DMPs[which(selected_DMPs$p.ebayes.BH_hypdriv <= 0.05),]
selected_DMPs.sig <- selected_DMPs.sig[order(selected_DMPs.sig$p.ebayes.BH_hypdriv, decreasing = F),]
selected_DMPs.sig.gr <- makeGRangesFromDataFrame(selected_DMPs.sig, keep.extra.columns = T, start.field = "pos", end.field = "pos")

dmr_plotter(name = "SP140", chr = "chr2", gmset = GMset.culled, annotation.gr = annotation.culled.gr, genome_version = "hg19", start_dmr = 231090640, end_dmr = 231090640, flanks = 1000)
dmr_plotter(name = "SPDEF", chr = "chr6", gmset = GMset.culled, annotation.gr = annotation.culled.gr, genome_version = "hg19", start_dmr = 34524597, end_dmr = 34524597, flanks = 1000)
dmr_plotter(name = "HDAC4", chr = "chr2", gmset = GMset.culled, annotation.gr = annotation.culled.gr, genome_version = "hg19", start_dmr = 239984030, end_dmr = 239984030, flanks = 1000)

# #Defining the gene regions myself
# gene.transcripts <- annotateTranscripts(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)
# GWAS_genes.granges <- gene.transcripts[which(gene.transcripts$Gene == "NOD2"),]
# GWAS_gene <- GWAS_genes.granges[which(width(GWAS_genes.granges) == max(width(GWAS_genes.granges))),]
# #Add 500 bps up- and downstream of the gene
# GWAS_gene_flanks <- GWAS_gene + 500
# which(annotation.culled.gr)
# start(GWAS_gene_flanks)
# end(GWAS_gene_flanks)

##############################################################################################################################################################
# Functions

############################
# Find the DMPs in the DMR #
############################

dmr_dmps <- function(dmr, annotation.gr){
  
  start_dmr <- start(dmr)
  end_dmr <- end(dmr)
  chr <- as.character(seqnames(dmr))

  dmr.hg19.annotation <- annotation.gr[which(seqnames(annotation.gr) == chr),]
  dmr.hg19.annotation <- dmr.hg19.annotation[which(start(dmr.hg19.annotation) >= start_dmr & start(dmr.hg19.annotation) <= end_dmr),]
  print(dmr.hg19.annotation)
}

################################
# DMR Plotter: Visualizes DMRs #
################################

dmr_plotter <- function(name, chr, gmset, annotation.gr, pheno, start_dmr, end_dmr, flanks = 1000, exp_data = NULL, genome_version = "hg19", pdf = F){
  #Check DMR
  dmr.hg19.annotation <- annotation.gr[which(seqnames(annotation.gr) == chr),]
  dmr.hg19.annotation <- dmr.hg19.annotation[which(start(dmr.hg19.annotation) >= start_dmr & start(dmr.hg19.annotation) <= end_dmr),]
  dmr.hg19.annotation.clean <- dmr.hg19.annotation
  print(dmr.hg19.annotation.clean)
  
  #DMR track
  dmrtrack <- AnnotationTrack(start = start_dmr, end = end_dmr, name = "DMR", chromosome = chr, genome = genome_version)
  
  #DMP track
  dmp.hg19.annotation <- dmr.hg19.annotation
  
  #Add the G and direction of each CpG
  start(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "-"]) <- start(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "-"])-1
  end(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "+"]) <- end(dmp.hg19.annotation[strand(dmp.hg19.annotation) == "+"])+1
  dmptrack <- AnnotationTrack(dmp.hg19.annotation, name = "DMP", col = dmp.hg19.annotation$red_end)
  
  #Add flanks
  start_dmr <- start_dmr-flanks
  end_dmr <- end_dmr+flanks
  dmr.hg19.annotation <- annotation.gr[which(seqnames(annotation.gr) == chr),]
  dmr.hg19.annotation <- dmr.hg19.annotation[which(start(dmr.hg19.annotation) >= start_dmr & start(dmr.hg19.annotation) <= end_dmr),]
  
  #Betas
  GMset.beta <- getBeta(gmset)
  dmr.450k <- data.frame(seqnames = as.character(seqnames(dmr.hg19.annotation)), pos = start(dmr.hg19.annotation), GMset.beta[names(dmr.hg19.annotation),])
  colnames(dmr.450k) <- c("seqnames", "pos", pData(gmset)$Sample_Name)
  dmr.450k.gr <- makeGRangesFromDataFrame(df = dmr.450k, keep.extra.columns = T, start.field = "pos", end.field = "pos")
  
  #Cohort
  cohort.450k <- pheno
  
  #Plot
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome_version, chromosome = chr)
#   grtrack <- UcscTrack(genome = genome_version, 
#                        chromosome = chr, 
#                        track = "knownGene", 
#                        from = start_dmr, 
#                        to = end_dmr, 
#                        trackType = "GeneRegionTrack", 
#                        rstarts = "exonStarts", 
#                        rends = "exonEnds", 
#                        gene = "name", 
#                        symbol = "name2", 
#                        transcript = "name", 
#                        strand = "strand", 
#                        fill = "#8282d2", 
#                        name = "UCSC Genes",
#                        showId = F)
  snptrack <- UcscTrack(genome = genome_version, 
                        chromosome = chr, 
                        track = "Common SNPs(141)", 
                        from = start_dmr, 
                        to = end_dmr, 
                        trackType = "AnnotationTrack", 
                        start = "chromStart", 
                        end = "chromEnd", 
                        gene = "name", 
                        feature = "func", 
                        strand = "strand", 
                        shape = "box",
                        stacking = "dense",
                        fill = "black", 
                        name = "SNPs")
  dtrack.450k <- DataTrack(range = dmr.450k.gr, name = "450K", ylim = c(0,1), groups = as.factor(cohort.450k), type = c("a", "p", "g"), legend = TRUE)
  
  if(pdf){
    pdf(file = paste0(getwd(), "/", name, ".pdf"))
  }
  plotTracks(list(itrack, 
                  gtrack, 
                  # grtrack, 
                  dtrack.450k, 
                  dmrtrack, 
                  dmptrack, 
                  snptrack), 
             from = start_dmr, 
             to = end_dmr, 
             main = name,
             cex.title = 0.9,
             cex.axis = 1,
             add53 = T,
             add35 = T,
             littleTicks = T)
  if(pdf){
    dev.off()  
  }
  
  return(dmr_annotation = dmr.hg19.annotation.clean)
}

###################################################################
# Quality Control: Taken from MethylAid without the visualization #
###################################################################

methylAid_nonvisual <- function(RGset, MU_threshold = 10.50, NP_threshold = 11.75, BS_threshold = 12.75, HC_threshold = 13.25, DP_threshold = 0.95){
  #Calculate detection p-value and frequency of probe passing threshold
  DP <- detectionP(RGset)
  DPfreq <- colSums(DP < 0.01, na.rm=TRUE)/nrow(DP)
  
  #Medians of the overall methylated and unmethylated intensities
  MU.full <- matrix(0.0, nrow = 2, ncol = ncol(RGset))
  M.full <- getMeth(Mset.raw)
  U.full <- getUnmeth(Mset.raw)
  MU.full[1,] <- colMedians(M.full, na.rm = TRUE)
  MU.full[2,] <- colMedians(U.full, na.rm = TRUE)
  colnames(MU.full) <- colnames(RGset)
  rownames(MU.full) <- c("Methylated", "Unmethylated")
  
  #Find the red and green intensities of the control probes
  data(hm450.controls, package="FDb.InfiniumMethylation.hg19", envir=environment())
  Red.full <- getRed(RGset)
  Green.full <- getGreen(RGset)
  id <- intersect(hm450.controls$Address, rownames(Red.full))
  Red.controls <- Red.full[rownames(Red.full) %in% id,]
  Green.controls <- Green.full[rownames(Green.full) %in% id,]
  hm450.controls <- hm450.controls[hm450.controls$Address %in% id,]
  hm450.controls <- hm450.controls[order(hm450.controls$Address), ]
  
  Red.controls.log <- log2(Red.controls)
  Green.controls.log <- log2(Green.controls)
  
  hm450.controls <- hm450.controls[!(hm450.controls$Type %in% c("NORM_A", "NORM_G", "NORM_C", "NORM_T")), ]
  control.data <- data.frame(Address=rep(rownames(Red.controls.log), ncol(Red.controls.log)), Samples=rep(colnames(Red.controls.log), each=nrow(Red.controls.log)), IntRed=as.vector(Red.controls.log),IntGrn=as.vector(Green.controls.log))
  hm450.control.data <- merge(hm450.controls, control.data)
  
  #Methylated/Unmethylated
  MU.full.log2 <- t(as.data.frame(log2(MU.full)))
  MU.outliers <- names(which(MU.full.log2[,1] <= MU_threshold))
  
  #Non-Polymorphic
  NP.control <- hm450.control.data[grepl("^NON-POLYMORPHIC$", hm450.control.data$Type),]
  NP.green <- NP.control[NP.control$Name %in% c("NP (C)", "NP (G)"), c(1:5,7)]
  NP.avg.green <- tapply(NP.green$IntGrn, NP.green$Samples, mean)
  NP.outliers <- names(which(NP.avg.green <= NP_threshold))
  
  #Bisulfite Conversion I
  BS.control <- hm450.control.data[grepl("^BISULFITE CONVERSION I$", hm450.control.data$Type),]
  BS.green <- BS.control[grepl("C1|C2|C3", BS.control$Name), c(1:5,7)]
  BS.avg.green <- tapply(BS.green$IntGrn, BS.green$Samples, mean)
  BS.outliers <- names(which(BS.avg.green <= BS_threshold))
  
  #Hybridization
  HC.control <- hm450.control.data[grepl("^HYBRIDIZATION$", hm450.control.data$Type),]
  HC.control <- HC.control[order(HC.control$Samples),]
  HC.green <- 0.5*(HC.control$IntGrn[grepl("High", HC.control$Name)] + HC.control$IntGrn[grepl("Low", HC.control$Name)])
  names(HC.green) <- HC.control$Samples[grepl("High", HC.control$Name)]
  HC.outliers <- names(which(HC.green <= HC_threshold))
  
  #Detection P-values
  DP.outliers <- names(which(DPfreq <= DP_threshold))
  
  #Total outliers
  total.outliers <- unique(c(MU.outliers, NP.outliers, HC.outliers, DP.outliers, BS.outliers))
  
  cat("Methylated outliers are:", MU.outliers, 
      "\nNon-polymorphic outliers are:", NP.outliers, 
      "\nBisulfite Conversion outliers are:", BS.outliers, 
      "\nHybridization outliers are:", HC.outliers, 
      "\nP-value outliers are:", DP.outliers, 
      "\n\nTotal outliers are:", total.outliers,
      "\nSample names:", phenoDataFrame[total.outliers,"Sample_Name"])
  
  #Find the outliers, remove them from the samplesheet and reread them
  indices.outliers <- which(targets$Sample_Name %in% phenoDataFrame[total.outliers,"Sample_Name"])
  targets.correct <- targets[-indices.outliers,]
  
  #Return the corrected target sheet
  return(targets.correct = targets.correct)
}

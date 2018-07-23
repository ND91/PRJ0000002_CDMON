#Find the SXS codes

pheno <- read.csv('/home/ayliyim/Dropbox/Epimac/Data/Crohn/450k/Monocyte Faltose/Monocytes/monocyte_450k_phenotype_V2.csv')[,-c(1,2)]
h450k <- read.csv('/media/ayliyim/ND_Backup/Data_Backup/102645_1D235-b5/Sample information/SampleSheet_Infinium_HDMethylation_102645-Batch001_102350-Batch005_Henneman.csv', skip = 10)

#Monocytes
pheno$inclusionnr <- gsub("CDACT", "CDACT-", pheno$inclusionnr)
pheno$inclusionnr <- gsub("CDREM", "CDREM-", pheno$inclusionnr)
monocytes <- h450k[h450k$Customer_Sample_Name %in% pheno$inclusionnr,]

monocytes <- dplyr::arrange(monocytes, Customer_Sample_Name)
pheno <- dplyr::arrange(pheno, inclusionnr)

Monocytes_Pheno <- data.frame(Sentrix_ID = monocytes$Sentrix_ID, Sentrix_Position = monocytes$Sentrix_Position, pheno)
write.csv(Monocytes_Pheno, "Phenodata_Monocytes.csv")

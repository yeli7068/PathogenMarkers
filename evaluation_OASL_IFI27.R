# License: GPLv3
# Author: Yang LI 
# Email: yeli7068@outlook.com
# install.packages("GEOquery", "dplyr", "BiocManager", "ggpubr")
# BiocManager::install("pROC", "illuminaHumanv3.db", "hgu133a.db", "clusterProfiler")
library(GEOquery)
library(dplyr)
library(ggpubr)
library(pROC)
library(illuminaHumanv3.db)  #GPL6947 GSE42046
library(hgu133a.db) # GPL96 GSE6269
library(clusterProfiler)

annotProbes <- function(expr, symbol){
    if (typeof(symbol) != "list") {
        ids <- toTable(symbol)
    } else (
        ids <- symbol
    )
    print(nrow(expr))
    expr <- expr[rownames(expr) %in% ids$probe_id, ]
    ids <- ids[match(rownames(expr), ids$probe_id),]
    tmp <- by(expr,
              ids$symbol,
              function(x) rownames(x)[which.max(rowMeans(x))])
    probes <- as.character(tmp)
    expr <- expr[rownames(expr) %in% probes, ] 
    rownames(expr) <- ids[match(rownames(expr), ids$probe_id), 2]
    print(nrow(expr))
    return(expr)
}


log2expr <- function(expr){
    # modified from GEO2R http://ncbi.nlm.nih.gov/geo/geo2r/
    ex <- expr
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    expr <- log2(ex) }
    return(expr)
}

plotGeneExpression <- function(expr_annoted, traits, gene, GSE) {
    library(ggpubr)
    traits$expression <- expr_annoted[gene, ]
    p <- ggboxplot(traits,
                   x = "pathogen",
                   y = "expression",
                   col = "pathogen",
                   palette = "jco",
                   add = "jitter",
                   title = paste("Expression of", gene, "in", GSE)
    )
    return(p)
}


plotROCgenes <- function(expr_annoted, traits, target = "flu", gene1 = "IFI27", gene2 ="OASL", GSE) {
    library(pROC)
    
    traits$roc_tag <- ifelse(grepl(target, traits$pathogen),"1", "0")
    traits$expression <- expr_annoted[gene1, ]
    
    roc1 <- roc(traits$roc_tag, traits$expression, ci =T, direction = "<")
    traits$expression <- expr_annoted[gene2, ]
    roc2 <- roc(traits$roc_tag, traits$expression, ci =T, direction = "<")
    #
    p <- ggroc(list(roc1,roc2), aes = "colour",legacy.axes = T) + 
        scale_color_discrete(name = "Gene", labels=c(paste0(gene1," AUC:", roc1$auc %>% round(digits = 3)), paste0(gene2," AUC:", roc2$auc %>% round(digits = 3)))) + 
        theme(legend.justification=c(1,0), legend.position=c(1,0)) +
        ggtitle(paste("Influenza versus non-Influenza on \n", GSE)) +
        theme(plot.title = element_text(hjust = 0.5))
    return(p)
}

if(T){
    # GSE6269 FigS
    GSE6269 <- getGEO("GSE6269", AnnotGPL = F, getGPL = F, destdir = "./")
    GSE6269expr <- exprs(GSE6269[[3]])
    GSE6269trait <- pData(GSE6269[[3]])
    library(hgu133a.db)
    GSE6269trait$pathogen <- GSE6269trait$`Pathogen:ch1`
    GSE6269trait$pathogen <- gsub("\\w\\..*", "Bacterial Infection",GSE6269trait$pathogen, perl = T)
    GSE6269trait$pathogen <- gsub(".*None.*", "Control", GSE6269trait$pathogen)
    GSE6269expr <- annotProbes(GSE6269expr, hgu133aSYMBOL)
    GSE6269expr <- log2expr(GSE6269expr)
    
    ROC_6269_P <- plotROCgenes(GSE6269expr, GSE6269trait, GSE = "GSE6269")
    
    GENE_IFI27_6269_P <- plotGeneExpression(expr_annoted = GSE6269expr, 
                                            traits = GSE6269trait, 
                                            gene = "OASL", GSE = "GSE6269")
    
    GENE_OASL_6269_P <- plotGeneExpression(expr_annoted = GSE6269expr, 
                                           traits = GSE6269trait, 
                                           gene = "IFI27", GSE = "GSE6269")
    
    multoplot(ROC_6269_P, GENE_OASL_6269_P, GENE_IFI27_6269_P, GSE = "GSE6269")
    
}

if(T){
    GSE42026 <- getGEO("GSE42026", AnnotGPL = F, getGPL = F, destdir = "./")
    GSE42026expr <- GSE42026[[1]] %>% exprs()
    GSE42026trait <- GSE42026[[1]] %>% pData()
    
    GSE42026trait$pathogen <- GSE42026trait$`infecting pathogen:ch1`
    GSE42026trait <- GSE42026trait[grep("bacterial", GSE42026trait$pathogen, invert = T),]
    GSE42026trait$pathogen <- gsub(".*Influenza.*", "Influenza A", GSE42026trait$pathogen)
    GSE42026trait$pathogen <- gsub(".*none.*", "Control", GSE42026trait$pathogen)
    
    GSE42026expr <- GSE42026expr[,rownames(GSE42026trait)]
    library(illuminaHumanv3.db)
    GSE42026expr <- annotProbes(GSE42026expr, illuminaHumanv3SYMBOL)
    GSE42026expr <- log2expr(GSE42026expr)
    
    GENE_IFI27_42026_P <- plotGeneExpression(expr_annoted = GSE42026expr, 
                                             traits = GSE42026trait, 
                                             gene = "IFI27", GSE = "GSE42026") 
    
    GENE_OASL_42026_P <- plotGeneExpression(expr_annoted = GSE42026expr, 
                                            traits = GSE42026trait, 
                                            gene = "OASL", GSE = "GSE42026")
    
    ROC_42026_p <- plotROCgenes(GSE42026expr, GSE42026trait, gene2 = "OASL", GSE = "GSE42026")
    
    multoplot(ROC_42026_p, GENE_OASL_42026_P, GENE_IFI27_42026_P, GSE="42026")
    
}

if(T){
    GSE38900 <- getGEO("GSE38900", AnnotGPL = F, getGPL = F, destdir = "./")
    GPL6884 <- getGEO("GPL6884", destdir="./")
    GPL6884SYMBOL <- Table(GPL6884)[,c("ID", "Entrez_Gene_ID", "Symbol")]
    GPL6884SYMBOL <- na.omit(GPL6884SYMBOL)
    GPL6884SYMBOL <- GPL6884SYMBOL[,c("ID", "Symbol")]
    colnames(GPL6884SYMBOL) <- c("probe_id", "symbol")
    
    GSE38900expr <- GSE38900[[2]] %>% exprs()
    GSE38900trait <- GSE38900[[2]] %>% pData()
    
    GSE38900trait$pathogen <- GSE38900trait$source_name_ch1
    GSE38900trait <- GSE38900trait[grep("follow up", GSE38900trait$pathogen, invert = T), ]
    GSE38900trait <- GSE38900trait[grep("Finnish", GSE38900trait$pathogen, invert = T), ]
    
    GSE38900trait$pathogen <- gsub(".*Influenza.*", "Influenza A", GSE38900trait$pathogen)
    GSE38900trait$pathogen <- gsub(".*HRV*", "HRV", GSE38900trait$pathogen)
    GSE38900trait$pathogen <- gsub(".*RSV*", "RSV", GSE38900trait$pathogen)
    GSE38900trait$pathogen <- gsub(".*healthy.*", "Control", GSE38900trait$pathogen)
    
    GSE38900expr <- GSE38900expr[,rownames(GSE38900trait)]
    
    GSE38900expr <- annotProbes(GSE38900expr, GPL6884SYMBOL)
    GSE38900expr <- log2expr(GSE38900expr)
    
    
    
    GENE_IFI27_38900_P <- plotGeneExpression(expr_annoted = GSE38900expr, 
                                             traits = GSE38900trait, 
                                             gene = "IFI27", GSE = "GSE38900") 
    
    GENE_OASL_38900_P <- plotGeneExpression(expr_annoted = GSE38900expr, 
                                            traits = GSE38900trait, 
                                            gene = "OASL", GSE = "GSE38900")
    
    ROC_38900_P <- plotROCgenes(GSE38900expr, GSE38900trait, GSE = "GSE38900")
    multoplot(ROC_38900_P, GENE_OASL_38900_P,GENE_IFI27_38900_P, GSE = "GSE38900")
}

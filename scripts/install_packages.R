b <- "https://bioconductor.org/packages/3.22"
options(repos=c(BioCsoft=paste0(b,"/bioc"), BioCann=paste0(b,"/data/annotation"), CRAN="https://cloud.r-project.org"))
pk <- c("data.table","ggplot2","ggrepel","pheatmap","survival","survminer","rmarkdown","limma","edgeR","clusterProfiler","enrichplot","org.Hs.eg.db")
install.packages(pk[!sapply(pk, requireNamespace, quietly=TRUE)])

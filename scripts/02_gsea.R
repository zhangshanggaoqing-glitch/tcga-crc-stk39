# GSEA (GO BP + KEGG) on the limma moderated t-statistic, STK39-high vs low
suppressPackageStartupMessages({library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2)})
set.seed(42)
setwd(normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE))),"..")))
tt <- fread("results/DE_STK39_high_vs_low_limma_voom.csv")
tt[, entrez := mapIds(org.Hs.eg.db, ensembl, "ENTREZID", "ENSEMBL")]
tt <- tt[!is.na(entrez)][order(-abs(t))][!duplicated(entrez)]
rk <- sort(setNames(tt$t, tt$entrez), decreasing=TRUE)
g1 <- gseGO(rk, OrgDb=org.Hs.eg.db, ont="BP", minGSSize=15, maxGSSize=500, pvalueCutoff=0.05, seed=TRUE, eps=0)
g1 <- setReadable(g1, org.Hs.eg.db, "ENTREZID")
g2 <- tryCatch(setReadable(gseKEGG(rk, organism="hsa", minGSSize=15, pvalueCutoff=0.05, seed=TRUE, eps=0), org.Hs.eg.db, "ENTREZID"), error=function(e){message("KEGG failed: ",conditionMessage(e)); NULL})
out <- list(GSEA_GO_BP=g1, GSEA_KEGG=g2)
res <- list()
for (n in names(out)) if (!is.null(out[[n]])) {
  d <- as.data.frame(out[[n]]); fwrite(d, sprintf("results/%s.csv", n))
  res[[n]] <- list(n_sig=nrow(d), n_pos=sum(d$NES>0), n_neg=sum(d$NES<0),
                   top_pos=head(d[d$NES>0, c("Description","NES","p.adjust")][order(-d$NES[d$NES>0]),],8),
                   top_neg=head(d[d$NES<0, c("Description","NES","p.adjust")][order(d$NES[d$NES<0]),],8))
  # dot plots are drawn in scripts/03_enrich_plots.R
}
print(res); saveRDS(res, "results/gsea_summary.rds"); cat("DONE\n")

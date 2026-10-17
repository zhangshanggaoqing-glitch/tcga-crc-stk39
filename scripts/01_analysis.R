# TCGA COAD+READ: STK39 expression, high-vs-low DE, enrichment, survival
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(ggplot2); library(ggrepel)
  library(pheatmap); library(survival); library(survminer)
  library(clusterProfiler); library(org.Hs.eg.db)
})
set.seed(42)
setwd(here <- normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE))),"..")))
dir.create("results",FALSE); dir.create("figures",FALSE)
STK39 <- "ENSG00000198648"

## ---- 1. Load data (UCSC Xena GDC hub; values are log2(count+1)) ----
rd <- function(p){ x <- fread(sprintf("data/TCGA-%s.star_counts.tsv.gz",p)); m <- as.matrix(x[,-1]); rownames(m) <- x$Ensembl_ID; m }
lc <- cbind(rd("COAD"), rd("READ"))
counts <- round(2^lc - 1); rm(lc); gc()
rownames(counts) <- sub("\\..*","",rownames(counts))
pm <- fread("data/gencode.v36.annotation.gtf.gene.probemap")
pm[, ens := sub("\\..*","",id)]
sym <- setNames(pm$gene, pm$ens)

clin <- rbindlist(lapply(c("COAD","READ"), function(p){
  d <- fread(sprintf("data/TCGA-%s.clinical.tsv.gz",p), select=c("sample","gender.demographic","age_at_index.demographic","ajcc_pathologic_stage.diagnoses"))
  d[, project := p]}))
setnames(clin, c("sample","gender","age","stage","project"))
surv <- rbindlist(lapply(c("COAD","READ"), function(p) fread(sprintf("data/TCGA-%s.survival.tsv.gz",p))))

meta <- data.table(sample=colnames(counts))
meta[, code := substr(sample,14,15)][, patient := substr(sample,1,12)]
meta <- meta[code %in% c("01","11")]
meta[, type := ifelse(code=="01","Tumor","Normal")]
# one sample per patient per type (first vial alphabetically)
setorder(meta, sample); meta <- meta[!duplicated(meta[,.(patient,type)])]
counts <- counts[, meta$sample]
cat("Samples:", table(meta$type), "\n")

## ---- 2. Normalise ----
dge <- DGEList(counts); dge <- calcNormFactors(dge)
lcpm <- cpm(dge, log=TRUE, prior.count=1)
meta[, STK39 := lcpm[STK39, sample]]

## ---- 3. Tumor vs Normal ----
w_all <- wilcox.test(STK39 ~ type, data=meta)
pairs <- dcast(meta[patient %in% meta[type=="Normal",patient]], patient ~ type, value.var="STK39")[complete.cases(Tumor,Normal)]
w_pair <- wilcox.test(pairs$Tumor, pairs$Normal, paired=TRUE)
tn <- list(n_tumor=sum(meta$type=="Tumor"), n_normal=sum(meta$type=="Normal"),
  median_tumor=median(meta[type=="Tumor",STK39]), median_normal=median(meta[type=="Normal",STK39]),
  p_unpaired=w_all$p.value, n_pairs=nrow(pairs),
  median_pair_diff=median(pairs$Tumor-pairs$Normal), p_paired=w_pair$p.value,
  frac_pairs_lower_in_tumor=mean(pairs$Tumor<pairs$Normal))
print(tn)
p1 <- ggplot(meta, aes(type, STK39, fill=type)) + geom_violin(alpha=.3, colour=NA) +
  geom_boxplot(width=.25, outlier.size=.5) + scale_fill_manual(values=c(Normal="#4DBBD5",Tumor="#E64B35")) +
  labs(x=NULL, y="STK39 log2 CPM (TMM)", title="STK39 in TCGA COAD+READ",
       subtitle=sprintf("Tumor n=%d vs Normal n=%d; Wilcoxon p=%.2g", tn$n_tumor, tn$n_normal, tn$p_unpaired)) +
  theme_classic(14) + theme(legend.position="none")
ggsave("figures/fig1_boxplot_tumor_normal.pdf", p1, width=4.5, height=5); ggsave("figures/fig1_boxplot_tumor_normal.png", p1, width=4.5, height=5, dpi=300)
pl <- melt(pairs, id.vars="patient", variable.name="type", value.name="STK39")
p1b <- ggplot(pl, aes(type, STK39, group=patient)) + geom_line(alpha=.3) + geom_point(size=1) +
  labs(x=NULL, y="STK39 log2 CPM", title="Matched pairs", subtitle=sprintf("n=%d; paired Wilcoxon p=%.2g", tn$n_pairs, tn$p_paired)) + theme_classic(14)
ggsave("figures/fig1b_paired.pdf", p1b, width=4, height=5); ggsave("figures/fig1b_paired.png", p1b, width=4, height=5, dpi=300)

## ---- 4. High vs low STK39 within tumors (limma-voom) ----
tm <- meta[type=="Tumor"]
tm[, group := factor(ifelse(STK39 > median(STK39), "High", "Low"), levels=c("Low","High"))]
tm[, project := clin$project[match(sample, clin$sample)]]
tm[is.na(project), project := ifelse(sample %in% fread("data/TCGA-READ.survival.tsv.gz")$sample,"READ","COAD")]
# restrict to protein-coding genes (org.Hs.eg.db GENETYPE) to avoid pseudogene/multimapping artefacts
gt <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, rownames(dge), c("GENETYPE"), "ENSEMBL"))
pc <- unique(gt$ENSEMBL[gt$GENETYPE %in% "protein-coding"])
d2 <- dge[intersect(rownames(dge), pc), tm$sample]
design <- model.matrix(~ project + group, data=tm)
keep <- filterByExpr(d2, design); d2 <- d2[keep,, keep.lib.sizes=FALSE]; d2 <- calcNormFactors(d2)
v <- voom(d2, design); fit <- eBayes(lmFit(v, design))
tt <- topTable(fit, coef="groupHigh", n=Inf, sort.by="P")
tt$ensembl <- rownames(tt); tt$symbol <- sym[tt$ensembl]
tt <- tt[tt$ensembl != STK39, ]   # remove the grouping gene itself
fwrite(tt, "results/DE_STK39_high_vs_low_limma_voom.csv")
sig <- subset(tt, adj.P.Val < 0.05 & abs(logFC) > 1)
de <- list(n_genes_tested=nrow(tt), n_up=sum(sig$logFC>0), n_down=sum(sig$logFC<0),
           n_fdr05_any=sum(tt$adj.P.Val<0.05), n_high=sum(tm$group=="High"), n_low=sum(tm$group=="Low"))
print(de)
tt$cls <- with(tt, ifelse(adj.P.Val<0.05 & logFC>1,"Up", ifelse(adj.P.Val<0.05 & logFC< -1,"Down","NS")))
lab <- head(sig[order(sig$adj.P.Val),],20)
p2 <- ggplot(tt, aes(logFC, -log10(adj.P.Val), colour=cls)) + geom_point(size=.6, alpha=.6) +
  scale_colour_manual(values=c(Up="#E64B35",Down="#3C5488",NS="grey70")) +
  geom_vline(xintercept=c(-1,1), lty=2) + geom_hline(yintercept=-log10(.05), lty=2) +
  geom_text_repel(data=lab, aes(label=symbol), size=3, colour="black", max.overlaps=30) +
  labs(title="STK39-high vs STK39-low tumors", subtitle=sprintf("limma-voom; FDR<0.05 & |log2FC|>1: %d up, %d down", de$n_up, de$n_down),
       x="log2 fold change", y="-log10 FDR", colour=NULL) + theme_classic(13)
ggsave("figures/fig2_volcano.pdf", p2, width=6.5, height=5.5); ggsave("figures/fig2_volcano.png", p2, width=6.5, height=5.5, dpi=300)
top <- head(tt[order(tt$adj.P.Val),], 50)
hm <- v$E[top$ensembl, order(tm$STK39)]; rownames(hm) <- top$symbol
hm <- t(scale(t(hm))); hm[hm>2] <- 2; hm[hm< -2] <- -2
ann <- data.frame(STK39=tm$group[order(tm$STK39)], Project=tm$project[order(tm$STK39)], row.names=colnames(hm))
for (ext in c("pdf","png")) pheatmap(hm, cluster_cols=FALSE, show_colnames=FALSE, annotation_col=ann, fontsize_row=6,
  color=colorRampPalette(c("#3C5488","white","#E64B35"))(100), main="Top 50 DE genes (z-score)",
  filename=sprintf("figures/fig3_heatmap.%s",ext), width=7, height=8)

## ---- 5. Enrichment ----
universe <- na.omit(mapIds(org.Hs.eg.db, tt$ensembl, "ENTREZID", "ENSEMBL"))
enr <- list()
for (dir in c("Up","Down")) {
  g <- na.omit(mapIds(org.Hs.eg.db, tt$ensembl[tt$cls==dir], "ENTREZID", "ENSEMBL"))
  if (length(g) < 10) next
  ego <- enrichGO(g, OrgDb=org.Hs.eg.db, ont="BP", universe=universe, readable=TRUE)
  enr[[paste0("GO_",dir)]] <- ego
  ek <- tryCatch(enrichKEGG(g, organism="hsa", universe=universe), error=function(e){message("KEGG failed: ",conditionMessage(e)); NULL})
  if (!is.null(ek)) enr[[paste0("KEGG_",dir)]] <- setReadable(ek, org.Hs.eg.db, "ENTREZID")
}
enr_n <- sapply(enr, function(e) sum(as.data.frame(e)$p.adjust < 0.05))
print(enr_n)
for (n in names(enr)) {
  fwrite(as.data.frame(enr[[n]]), sprintf("results/enrich_%s.csv", n))
  # dot plots are drawn in scripts/03_enrich_plots.R
}

## ---- 6. Survival ----
sv <- merge(tm, surv[, .(sample, OS.time, OS)], by="sample")
sv <- merge(sv, clin[, .(sample, gender, age, stage)], by="sample", all.x=TRUE)
sv <- sv[!is.na(OS.time) & OS.time > 0]
sv[, stage_grp := fifelse(grepl("^Stage (III|IV)", stage), "III-IV", fifelse(grepl("^Stage (I|II)", stage), "I-II", NA_character_))]
sv[, STK39_z := as.numeric(scale(STK39))]
fitkm <- survfit(Surv(OS.time/30.44, OS) ~ group, data=sv)
lr <- survdiff(Surv(OS.time, OS) ~ group, data=sv)
p_lr <- 1 - pchisq(lr$chisq, 1)
km <- ggsurvplot(fitkm, data=sv, pval=TRUE, risk.table=TRUE, palette=c("#3C5488","#E64B35"),
  legend.labs=c("STK39 low","STK39 high"), xlab="Months", title="Overall survival, TCGA COAD+READ")
for (ext in c("pdf","png")) { if (ext=="pdf") pdf("figures/fig5_KM.pdf", 6.5, 6.5, onefile=FALSE) else png("figures/fig5_KM.png", 6.5, 6.5, units="in", res=300); print(km); dev.off() }
cox_uni <- coxph(Surv(OS.time, OS) ~ STK39_z, data=sv)
cox_grp <- coxph(Surv(OS.time, OS) ~ group, data=sv)
cox_mv  <- coxph(Surv(OS.time, OS) ~ STK39_z + age + gender + stage_grp + project, data=sv)
zph <- cox.zph(cox_mv)
sumcox <- function(m, term) { s <- summary(m); c(HR=s$conf.int[term,1], lo=s$conf.int[term,3], hi=s$conf.int[term,4], p=s$coefficients[term,5], n=m$n, events=m$nevent) }
cox <- rbind(univariate_continuous=sumcox(cox_uni,"STK39_z"), univariate_high_vs_low=sumcox(cox_grp,"groupHigh"),
             multivariable_continuous=sumcox(cox_mv,"STK39_z"))
print(cox); print(zph)
capture.output(summary(cox_mv), zph, file="results/cox_multivariable_full.txt")
write.csv(cox, "results/cox_summary.csv")
fwrite(meta, "results/sample_STK39_expression.csv")

saveRDS(list(tn=tn, de=de, enr_n=enr_n, cox=cox, p_logrank=p_lr, n_surv=nrow(sv), events=sum(sv$OS),
  top_up=head(sig[sig$logFC>0,][order(sig$adj.P.Val[sig$logFC>0]),"symbol"],15),
  top_down=head(sig[sig$logFC<0,][order(sig$adj.P.Val[sig$logFC<0]),"symbol"],15),
  enr_top=lapply(enr, function(e) head(as.data.frame(e)[,c("Description","p.adjust","Count")],8)),
  zph=zph$table), "results/summary.rds")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("DONE\n")

# Enrichment dot plots (ORA + GSEA) drawn from saved result tables, with wrapped labels
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(stringr)})
setwd(normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE))),"..")))
lab <- function(x, max=60, w=40) str_wrap(str_trunc(x, max), w)
save2 <- function(p, n, w, h) { ggsave(sprintf("figures/%s.pdf",n), p, width=w, height=h)
  ggsave(sprintf("figures/%s.png",n), p, width=w, height=h, dpi=300) }
th <- theme_bw(12) + theme(axis.text.y=element_text(size=10, lineheight=.85), plot.margin=margin(10,15,10,10),
  plot.title=element_text(size=13, face="bold"), legend.position="right")
# GSEA: top 8 positive + top 8 negative by |NES|, faceted by direction
for (n in c("GSEA_GO_BP","GSEA_KEGG")) {
  f <- sprintf("results/%s.csv", n); if (!file.exists(f)) next
  d <- fread(f)[p.adjust < 0.05]; if (!nrow(d)) next
  d[, dir := factor(ifelse(NES > 0, "Higher in STK39-high", "Higher in STK39-low"), levels=c("Higher in STK39-high","Higher in STK39-low"))]
  d <- rbind(d[NES>0][order(-NES)][1:min(8,.N)], d[NES<0][order(NES)][1:min(8,.N)])
  d[, term := lab(Description)]; d[, term := factor(term, levels=unique(term[order(NES)]))]
  p <- ggplot(d, aes(NES, term, size=setSize, colour=p.adjust)) + geom_point() +
    geom_vline(xintercept=0, lty=2, colour="grey50") +
    facet_grid(dir ~ ., scales="free_y", space="free_y") +
    scale_colour_gradient(low="#E64B35", high="#3C5488", name="FDR") + scale_size(range=c(2,7), name="Set size") +
    labs(x="Normalised enrichment score (NES)", y=NULL, title=sprintf("%s (top 8 per direction)", gsub("_"," ",n))) +
    th + theme(strip.text.y=element_text(size=10))
  save2(p, paste0("fig4_", n), 9, 9)
}
# ORA: all significant terms up to 10
for (n in c("GO_Up","GO_Down","KEGG_Up","KEGG_Down")) {
  f <- sprintf("results/enrich_%s.csv", n); if (!file.exists(f)) next
  d <- fread(f)[p.adjust < 0.05]; if (!nrow(d)) next
  d <- d[order(p.adjust)][1:min(10,.N)]
  d[, gr := sapply(strsplit(GeneRatio,"/"), function(z) as.numeric(z[1])/as.numeric(z[2]))]
  d[, term := lab(Description)]; d[, term := factor(term, levels=term[order(gr)])]
  p <- ggplot(d, aes(gr, term, size=Count, colour=p.adjust)) + geom_point() +
    scale_colour_gradient(low="#E64B35", high="#3C5488", name="FDR", n.breaks=3, labels=function(x) signif(x,2),
      guide=guide_colourbar(barwidth=unit(4,"cm"), title.vjust=.8)) + scale_size(range=c(3,7), breaks=scales::breaks_pretty(3)) +
    scale_x_continuous(expand=expansion(mult=.25)) +
    labs(x="Gene ratio", y=NULL, title=sprintf("ORA: %s", sub("_"," (",n)) |> paste0(")"), subtitle="DE genes, STK39-high vs low") + th +
    theme(legend.position="bottom", legend.box="horizontal", plot.title.position="plot")
  save2(p, paste0("fig4_", n), 7.5, max(4.5, 2.2 + 0.6*nrow(d)))
}
cat("DONE\n")

#!/usr/bin/env bash
set -euo pipefail; cd "$(dirname "$0")"
bash scripts/00_download.sh
Rscript scripts/01_analysis.R      # expression, DE, ORA, survival (~2 min, ~4 GB RAM)
Rscript scripts/02_gsea.R          # GSEA GO/KEGG (~5 min; KEGG needs internet)
Rscript scripts/03_enrich_plots.R  # enrichment dot plots
Rscript -e 'rmarkdown::render("report.Rmd", output_file="report.html")'

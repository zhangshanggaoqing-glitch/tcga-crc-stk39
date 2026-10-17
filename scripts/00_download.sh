#!/usr/bin/env bash
# Download public TCGA COAD/READ data from the UCSC Xena GDC hub
set -euo pipefail; cd "$(dirname "$0")/../data"
B=https://gdc-hub.s3.us-east-1.amazonaws.com/download
for p in COAD READ; do for f in star_counts.tsv.gz survival.tsv.gz clinical.tsv.gz; do curl -fLO "$B/TCGA-$p.$f"; done; done
curl -fLO "$B/gencode.v36.annotation.gtf.gene.probemap"

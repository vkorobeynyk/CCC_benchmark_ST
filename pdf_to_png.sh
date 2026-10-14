#!/bin/bash
# Convert every PDF in output/**/results/ to one PNG per page.
# output/spatialScattering/results/CosMx_HFC_diagnostic_plots.pdf
#   -> output_individual_images/spatialScattering/results/CosMx_HFC_diagnostic_plots/page-1.png, page-2.png, ...

set -euo pipefail

SRC="output"
DEST="output_individual_images"
DPI=300

command -v pdftoppm >/dev/null || { echo "pdftoppm not found (install poppler-utils)"; exit 1; }

find "$SRC" -type f -path "*/results/*" -name "*.pdf" -print0 |
while IFS= read -r -d '' pdf; do
  rel="${pdf#"$SRC"/}"                 # spatialScattering/results/CosMx_HFC_diagnostic_plots.pdf
  outdir="$DEST/${rel%.pdf}"           # output_individual_images/spatialScattering/results/CosMx_HFC_diagnostic_plots
  mkdir -p "$outdir"
  echo "Converting $pdf -> $outdir/"
  pdftoppm -png -r "$DPI" "$pdf" "$outdir/page"
done

echo "Done."

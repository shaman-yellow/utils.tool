pandoc --metadata reference-section-title="Reference" \
  --metadata link-citations=true \
  --filter pandoc-fignos \
  --filter pandoc-tablenos \
  --citeproc \
  --csl=/mnt/data/wizard/Documents/article/styles.csl \
  --bibliography=/mnt/data/wizard/Documents/article/library.bib \
  -f markdown+tex_math_dollars --mathjax $1 -o $1.tex
## pandoc -o custom-reference.docx --print-default-data-file reference.docx

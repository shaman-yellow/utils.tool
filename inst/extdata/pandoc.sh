pandoc --metadata reference-section-title="Reference" \
  --reference-doc custom-reference.docx \
  --metadata link-citations=true \
  --filter pandoc-fignos \
  --filter pandoc-tablenos \
  --citeproc \
  --csl=../styles.csl \
  --bibliography=../library.bib \
  --pdf-engine=xelatex \
  -f markdown+tex_math_dollars --mathjax $1 -o $2
## pandoc -o custom-reference.docx --print-default-data-file reference.docx

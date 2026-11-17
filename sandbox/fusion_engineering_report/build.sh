#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error main.tex > build/latex_pass1.log
bibtex main > build/bibtex.log
pdflatex -interaction=nonstopmode -halt-on-error main.tex > build/latex_pass2.log
pdflatex -interaction=nonstopmode -halt-on-error main.tex > build/latex_pass3.log
cp main.pdf fusion_nuclear_data_sensitivity_report.pdf
printf '%s\n' 'Built fusion_nuclear_data_sensitivity_report.pdf'

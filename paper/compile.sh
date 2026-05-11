#!/bin/sh
# compile.sh — Compile paper.tex with pdflatex + bibtex
# Created: 2026-05-09
TEXBIN="$HOME/Library/TinyTeX/bin/universal-darwin"
cd "$(dirname "$0")"

echo "=== Pass 1 ==="
$TEXBIN/pdflatex -interaction=nonstopmode paper.tex

echo "=== BibTeX ==="
$TEXBIN/bibtex paper

echo "=== Pass 2 ==="
$TEXBIN/pdflatex -interaction=nonstopmode paper.tex

echo "=== Pass 3 ==="
$TEXBIN/pdflatex -interaction=nonstopmode paper.tex

echo "=== Done ==="
grep "Output written" paper.log

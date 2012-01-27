#
# Makefile for guhdi
# 20120127
#

papers= B1 B2

default: B2.bbl

B2: B2.pdf
B1: B1.pdf

complete: bib


all: $(addsuffix .pdf,$(basename ${papers}))

bib: $(addsuffix .bbl,$(basename ${papers}))

%.pdf: %.tex $(wildcard *.tex)
	pdflatex $<
	pdflatex $<

%.bbl: %.pdf *.bib
	bibtex $(basename $<)
	pdflatex $(basename $<)
	pdflatex $(basename $<)

clean:
	rm -f *.bbl *.blg *.aux *.out *.log \
$(addsuffix .pdf,$(basename ${papers}))

cleanmax: clean
	rm -f *~ lablst.pdf


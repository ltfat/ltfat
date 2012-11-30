PDFNAME = ltfat.pdf

all:	$(PDFNAME)

%.pdf: %.tex Makefile $(shell ls *.bib)
	pslatex $*
	-bibtex $*
	pslatex $*
	pslatex $*
	dvips -t A4 $*.dvi
	ps2pdf $*.ps $*.pdf


%.pdf: %.lyx Makefile $(shell ls *.bib)
	lyx -f --export pdf2 $<

%.tex: %.lyx Makefile $(shell ls *.bib)
	lyx -f --export latex $<

texclean:	
	-rm *.aux *.bbl *.blg *.log *.dvi *.out

clean:	texclean
	-rm $(PDFNAME)

.PHONY: all texclean clean

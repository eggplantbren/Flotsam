FILENAMES = ms
PDFS = $(foreach f, $(FILENAMES), $(f).pdf)

# The default build rule
%.pdf: %.tex
	latex $*
	latex $*
	dvipdf $*

default: ms.pdf

clean:
	rm -f *.aux *.log *.dvi *.pdf


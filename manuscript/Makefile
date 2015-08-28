OUTNAME = valet-$(shell date +%F).pdf

all: valet.pdf

valet.pdf: valet.tex
	pdflatex valet.tex; bibtex valet; pdflatex valet.tex; pdflatex valet.tex;
	mv valet.pdf ${OUTNAME}

clean:
	rm -rf valet*.pdf
	rm -f valet.log
	rm -f valet.dvi
	rm -f valet.aux
	rm -f valet.bbl
	rm -f valet.blg
	rm -f valet.bcf

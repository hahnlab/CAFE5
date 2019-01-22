.PHONY: CAFE5.pdf all clean

all: CAFE5.pdf

CAFE5.pdf: CAFE5.tex
	latexmk -f -bibtex -pdf -pdflatex="pdflatex -interaction=nonstopmode" CAFE5.tex

clean:
	latexmk -CA


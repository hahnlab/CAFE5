.PHONY: CAFE5.pdf all clean

all: CAFE5.pdf

CAFE5.pdf: CAFE5.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make CAFE5.tex

clean:
	latexmk -CA


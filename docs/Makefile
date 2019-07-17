.PHONY: CAFE5.pdf all clean

all: CAFE5.pdf

CAFE5.pdf: CAFE5.tex README.tex
	latexmk -f -bibtex -pdf -pdflatex="pdflatex -interaction=nonstopmode" CAFE5.tex

README.tex: ../README.md
	pandoc $? -o $@

clean:
	latexmk -CA
	rm -f README.tex


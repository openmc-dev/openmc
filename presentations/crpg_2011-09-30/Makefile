all: slides.tex
	mkdir -p build
	pdflatex --output-directory=build slides.tex > /dev/null
	pdflatex --output-directory=build slides.tex > /dev/null

clean:
	rm -rf build

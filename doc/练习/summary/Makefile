
main = main
all_file = ${main}.tex 

pdf: ${main}.pdf

${main}.pdf: ${all_file}
	xelatex ${main}.tex
	bibtex  ${main}.aux
	xelatex ${main}.tex
	xelatex ${main}.tex
	evince  ${main}.pdf&

view:
	evince ${main}.pdf&

.PHONY:clean  
clean:
	-rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *~ *.tid *.tms *.pdf *.bak *.blg *.bbl


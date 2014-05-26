SOURCES=mcs.w common.w csg.w error.w main.w matrix.w mem.w sim.w vector.w world.w heap.w table.w types.w cuda.w postfix.w bstack.w

mcs: mcs.c
	gcc -g -Wall -o mcs mcs.c -lm

mcs.c: ${SOURCES}
	ctangle mcs.w

mcs.pdf: mcs.tex
	pdflatex mcs.tex

mcs.tex: ${SOURCES}
	cweave mcs.w

clean:
	rm -f *~ mcs.tex *.aux *.out *.ps *.dvi *.pdf *.scn *.log *.toc *.idx *.c mcs


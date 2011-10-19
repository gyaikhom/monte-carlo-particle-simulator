SOURCES=mcs.w common.w csg.w error.w main.w matrix.w mem.w sim.w vector.w world.w heap.w table.w

mcs: mcs.c
	gcc -g -Wall -o mcs mcs.c -lm

mcs.c: ${SOURCES}
	ctangle mcs.w

mcs.pdf: mcs.ps
	pstopdf mcs.ps

mcs.ps: mcs.tex
	tex mcs.tex
	dvips mcs.dvi

mcs.tex: ${SOURCES}
	cweave mcs.w

clean:
	rm -f *~ mcs.tex *.ps *.dvi *.pdf *.scn *.log *.toc *.idx *.c mcs


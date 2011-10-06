SOURCES=mcs.w common.w csg.w error.w main.w matrix.w mem.w sim.w vector.w world.w

test_input: test_input.c
	gcc -g -Wall -o test_input test_input.c -lm

test_input.c: ${SOURCES}
	ctangle mcs.w

test_containment: test_containment.c
	gcc -g -Wall -o test_containment test_containment.c -lm

test_containment.c: ${SOURCES}
	ctangle mcs.w

mcs.c: ${SOURCES}
	ctangle mcs.w

mcs.pdf: mcs.tex
	pdftex mcs.tex

mcs.tex: ${SOURCES}
	cweave mcs.w

clean:
	rm -f *~ mcs.tex *.scn *.log *.toc *.idx *.c mcs test_*


test_input: test_input.c
	gcc -g -Wall -o test_input test_input.c -lm

test_input.c:  mcs.w
	ctangle mcs.w

test_containment: test_containment.c
	gcc -g -Wall -o test_containment test_containment.c -lm

test_containment.c:  mcs.w
	ctangle mcs.w

mcs.c: mcs.w
	ctangle mcs.w

mcs.pdf: mcs.tex
	pdftex mcs.tex

mcs.tex: mcs.w
	cweave mcs.w

clean:
	rm -f *~ *.tex *.scn *pdf *.log *.toc *.idx *.c mcs test_*


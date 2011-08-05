mcs: mcs.c
	gcc -o mcs mcs.c

mcs.c: mcs.w
	ctangle mcs.w

mcs.pdf: mcs.tex
	pdftex mcs.tex

mcs.tex: mcs.w
	cweave mcs.w

clean:
	rm -f *~ *.tex *.scn *pdf *.log *.toc *.idx *.c


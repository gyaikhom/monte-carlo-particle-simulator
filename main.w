@** The main program.

@p
@h
@<Include preamble for applications@>@/

@ @<Include preamble for applications@>=
@<Include system libraries@>@/
@<Type definitions@>@/
@<Global variables@>@/
@<Forward declare functions@>@/
@<Global functions@>@/

@ @<Parse command line arguments@>=

@ @<Process input files@>=

@ @<Create physics tables@>=

@ @<Create and process events@>=

@ @<Clean up the system@>=
mem_free(mem);
mem_free(mem_p);

@ @<Include system libraries@>=
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

@ @<Global variables@>=
Vector positive_xaxis_unit_vector = { 1.0, 0.0, 0.0, 1.0 };

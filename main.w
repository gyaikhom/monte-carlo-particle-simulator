@** The main program.

@p
@h
@<Include preamble for applications@>@/
int main(int argc, char *argv[])
{
	@<Test geometry table generation@>;
	cuda_mem_typed_alloc(10,uint32_t,mem_gpu);
	cuda_mem_typed_alloc(10,float,mem_gpu);
	@<Clean up the system@>;
	return 0;
}

@ @<Include preamble for applications@>=
@<Include system libraries@>@/
@<Type definitions@>@/
@<Forward declare functions@>@/
@<Global variables@>@/
@<Global functions@>@/

@ @<Parse command line arguments@>=

@ @<Process input files@>=

@ @<Create physics tables@>=

@ @<Create and process events@>=

@ @<Clean up the system@>=
mem_free(mem_phase_one);
cuda_mem_free(mem_gpu); /* must be done before freeing phase two memory */
mem_free(mem_phase_two);

@ @<Include system libraries@>=
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cuda_runtime_api.h"

@ @<Global variables@>=
Vector positive_xaxis_unit_vector = { 1.0, 0.0, 0.0, 1.0 };


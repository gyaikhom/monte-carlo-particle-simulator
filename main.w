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
for (i = 0; i < num_events; i++) {
	@<Simulate $i$th event@>;
}

@ @<Clean up the system@>=
destroy_csg_tree(csg_tree.root);

@ @<Include system libraries@>=
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

@ @<Global variables@>=
uint32_t i, j, k, l, m, n; /* counters */
Event *current_event = NULL;
Vertex *current_vertex = NULL;
Particle *current_particle = NULL;
uint32_t num_events = 5;
uint32_t num_vertices = 5;
uint32_t num_particles = 5;
Particle_gun particle_gun;
Particle_stack particle_stack;
Vector temp_vector;
Vector positive_xaxis_unit_vector = { 1.0, 0.0, 0.0, 1.0 };

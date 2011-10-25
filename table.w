@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Preparing the shared tables.

@<Type definitions@>=
typedef struct {
    @<Counters used during geometry table generation@>;
    @<Define shared primitives table@>;
    @<Define shared solids table@>;
    @<Define shared subcuboids table@>;
} GeometryTable;
GeometryTable geotab;

@ @<Counters used during geometry table generation@>=
uint32_t ip; /* index within primitive table */
uint32_t is; /* index within solid table */
uint32_t ic; /* index within subcuboid table */
uint32_t ipf; /* index within postfix expression buffer */

@
@d primitives_table_size 1024 /* maximum number of primitives */
@<Define shared primitives table@>=
struct {
    Matrix a, i; /* affine and inverse */
    Primitive p; /* primtive data */
} p[primitives_table_size]; /* primitives table */

@
@d solids_table_size 512 /* maximum number of CSG solids */
@d postfix_buffer_size 4096 /* size of the postfix buffer */
@<Define shared solids table@>=
struct {
    uint32_t s, e; /* postfix expression start and end indices */
} s[solids_table_size]; /* solids table */
int32_t pf[postfix_buffer_size]; /* postfix expression buffer */

@
@d subcuboids_table_size 256
@d max_solids_per_subcuboid 100
@<Define shared subcuboids table@>=
struct {
    uint32_t n; /* number of solids inside subcuboid */
    uint32_t s[max_solids_per_subcuboid]; /* indices of solids */
} c[subcuboids_table_size]; /* subcuboids table */

@ @<Global functions@>=
void fill_geotab_with_csg(GeometryTable *g, CSG_Node *n) {
    if (NULL == n) return;
    if (SOLID == n->op) {
        matrix_copy(g->p[g->ip].a, n->affine);
	matrix_copy(g->p[g->ip].i, n->inverse);
	g->p[g->ip].p = *(n->leaf.p);
	g->pf[g->ipf++] = g->ip++;
        return;
    }
    fill_geotab_with_csg(g, n->internal.left);
    switch(n->op) {
    case UNION:
        g->pf[g->ipf++] = BOOLEAN_UNION;
        break;
    case INTERSECTION:
        g->pf[g->ipf++] = BOOLEAN_INTERSECTION;
        break;
    case DIFFERENCE:
        g->pf[g->ipf++] = BOOLEAN_DIFFERENCE;
        break;
    default: ;
    }
    fill_geotab_with_csg(g, n->internal.right);
}

@ @<Global functions@>=
void create_geotab(GeometryTable *g)
{
    int i, j;
    CSG_Node *s;
    @<Initialise the geometry table@>;
    for (i = 0; i < MAX_CSG_SOLIDS; ++i) {
        s = forest_of_solids.t[i].s;
        @<Fill table entries for this solid@>;
	@<Add solid to subcuboid if contained by subcuboid@>;
    }
}

@ @<Initialise the geometry table@>=
g->ip = g->is = g->ic = g->ipf = 0;
for (i = 0; i < subcuboids_table_size; ++i) g->c[i].n = 0;

@ @<Fill table entries for this solid@>=
if (NULL == s) continue;
g->s[g->is].s = g->ipf;
fill_geotab_with_csg(g, s);
g->s[g->is].e = g->ipf;

@ @<Add solid to subcuboid if contained by subcuboid@>=
for (j = 0; j < MAX_SUBCUBOIDS; ++j) {
    if (no_intersection_bb(subcuboids[j].bb, s->bb)) continue;
    g->c[j].s[g->c[j].n++] = g->is++;
}

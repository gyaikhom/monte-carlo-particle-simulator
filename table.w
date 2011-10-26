@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Preparing the shared tables.

\bigskip

\centerline{\epsfig{file=figures/geometry-solids,scale=1}}

\bigskip

@<Type definitions@>=
typedef struct {
    @<Sizes of the geometry table components@>;
    @<Counters used during geometry table generation@>;
    struct primitives_table_item *p;
    struct solids_table_item *s;
    struct subcuboids_table_item *c;
    int32_t *pb; /* pointer to the postfix expression buffer */
    int32_t *sb; /* pointer to the solid indices buffer */
} GeometryTable;
GeometryTable geotab;

@

\bigskip

\centerline{\epsfig{file=figures/geometry-primitives-table,scale=1}}

\bigskip

@<Type definitions@>=
struct primitives_table_item {
    Matrix a, i; /* accumulated affine and inverse transformations */
    Primitive p; /* primtive data */
};

@

\bigskip

\centerline{\epsfig{file=figures/geometry-solids-table,scale=1}}

\bigskip

@<Type definitions@>=
struct solids_table_item {
    uint32_t s, c; /* start index and item count in postfix expression buffer */
};

@

\bigskip

\centerline{\epsfig{file=figures/geometry-subcuboids-table,scale=1}}

\bigskip

@<Shared subcuboids table@>=
struct subcuboids_table_item {
    uint32_t s, c;  /* start index and item count in solid indices buffer */
};

@ @<Sizes of the geometry table components@>=
uint32_t np; /* number of entries in primitives table */
uint32_t ns; /* number of entries in solids table */
uint32_t nc; /* number of entries in subcuboids table */
uint32_t npb; /* number of entries in postfix expression buffer */
uint32_t nsb; /* number of entries in solid indices buffer */

@ @<Counters used during geometry table generation@>=
uint32_t ip; /* index within primitives table */
uint32_t is; /* index within solids table */
uint32_t ic; /* index within subcuboids table */
uint32_t ipb; /* index within postfix expression buffer */
uint32_t isb; /* index within solid indics buffer */

@ @<Global functions@>=
void fill_geotab_with_csg(GeometryTable *g, CSG_Node *n) {
    if (NULL == g || NULL == n) return;
    if (is_primitive(n->op)) {
        matrix_copy(g->p[g->ip].a, n->affine);
	matrix_copy(g->p[g->ip].i, n->inverse);
	g->p[g->ip].p = *(n->leaf.p);
	g->pb[(g->ipb)++] = (g->ip)++;
        return;
    }
    fill_geotab_with_csg(g, n->internal.left);
    fill_geotab_with_csg(g, n->internal.right);
    switch(BIT_MASK_NODE & n->op) {
    case UNION:
        g->pb[g->ipb++] = BOOLEAN_UNION;
        break;
    case INTERSECTION:
        g->pb[g->ipb++] = BOOLEAN_INTERSECTION;
        break;
    case DIFFERENCE:
        g->pb[g->ipb++] = BOOLEAN_DIFFERENCE;
        break;
    default: ;
    }
}

@ @<Global functions@>=
void create_geotab(GeometryTable *g)
{
    int i;
    CSG_Node *s;
    @<Initialise the geometry table@>;
    for (i = 0; i < forest_of_solids.n; ++i) {
        s = forest_of_solids.s[i];
        @<Fill table entries for this solid@>;
    }
}

@ @<Initialise the geometry table@>=
g->ip = g->is = g->ic = g->ipb = g->isb = g->nc = g->nsb = 0;
g->np = nodes_repo->stat[PRIMITIVE];
g->npb = nodes_repo->stat[PRIMITIVE] +
       nodes_repo->stat[UNION] +
       nodes_repo->stat[INTERSECTION] +
       nodes_repo->stat[DIFFERENCE];
g->ns = forest_of_solids.n;
g->p = mem_typed_alloc(g->np, struct primitives_table_item, mem_phase_two);
g->pb = mem_typed_alloc(g->npb, int32_t, mem_phase_two); 
g->s = mem_typed_alloc(g->ns, struct solids_table_item, mem_phase_two);

@ @<Fill table entries for this solid@>=
if (NULL == s) continue;
g->s[g->is].s = g->ipb;
fill_geotab_with_csg(g, s);
g->s[g->is].c = g->ipb - g->s[g->is].s;
++(g->is);

@ @<Add solid to subcuboid if contained by subcuboid@>=
for (j = 0; j < MAX_SUBCUBOIDS; ++j) {
    if (no_intersection_bb(subcuboids[j].bb, s->bb)) continue;
    g->c[j].s[g->c[j].n++] = g->is++;
}

@ @<Global functions@>=
void print_geotab(FILE *f, GeometryTable *g)
{
	uint32_t i, j;
	fprintf(f, "G E O M E T R Y  T A B L E\nPrimitives table:\n");
	for (i = 0; i < g->np; ++i) {
	    fprintf(f, "Affine:\n");
	    matrix_print(f, g->p[i].a, 4, 4, 0);
	    fprintf(f, "\nInverse:\n");
	    matrix_print(f, g->p[i].i, 4, 4, 0);
	    for (j = 0; j < 80; ++j) fprintf(f, "-");
	    fprintf(f, "\n");
	}
	fprintf(f, "Solids table:\n");
	for (i = 0; i < g->ns; ++i) fprintf(f, "%u %u\n", g->s[i].s, g->s[i].c);
	fprintf(f, "Postfix buffer:\n");
	for (i = 0; i < g->npb; ++i) fprintf(f, "%d ", g->pb[i]);
	fprintf(f, "\n");
}

@ @<Test geometry table generation@>=
{
        if (false == read_geometry("input.dat")) exit(1);
	print_geom_statistics(stdout);
	create_geotab(&geotab);
	print_geotab(stdout, &geotab);
}
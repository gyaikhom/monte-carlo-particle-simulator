@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Preparing the shared tables.

\bigskip

\centerline{\epsfig{file=figures/geometry-solids,scale=1}}

\bigskip

@<Type definitions@>=
struct geomtab_struct {
    @<Sizes of the geometry table components@>;
    @<Counters used during geometry table generation@>;
    @<Subcuboids related data structures inside the geometry table@>;
    struct primitives_table_item *p;
    struct solids_table_item *s;
    int32_t *pb; /* pointer to the postfix expression buffer */
    uint32_t *sb; /* pointer to the solid indices buffer */
    BoundingBox sw; /* the simulation world cuboid */
    uint32_t l, m, n; /* divisions along $x$, $y$ and $z$ axes */
} geotab;

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
void fill_geotab_csg_table(GeometryTable *g, CSG_Node *n) {
    if (NULL == g || NULL == n) return;
    if (is_primitive(n->op)) {
        matrix_copy(g->p[g->ip].a, n->affine);
	matrix_copy(g->p[g->ip].i, n->inverse);
	g->p[g->ip].p = *(n->leaf.p);
	g->pb[(g->ipb)++] = (g->ip)++;
        return;
    }
    fill_geotab_csg_table(g, n->internal.left);
    fill_geotab_csg_table(g, n->internal.right);
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
bool fill_geotab_subcuboids_table(GeometryTable *g)
{
    Area t, r; /* temporary paged solids indices buffer */
    uint32_t m = 3; /* maximum number of items per page */
    uint32_t c = 0; /* number of items in current page */
    CSG_Node *s = NULL;
    uint32_t *sb; /* current solid indices buffer page */
    uint32_t i, j;
    mem_init(t);
    mem_init(r);
    @<Fill in the subcuboids table and create a paged solid indices buffer@>;
    @<Finalise solid indices buffer by moving paged data to contiguous memory@>;
    return true;

exit_error: mem_free(t);
    return false;
}

@ @<Fill in the subcuboids table and create a paged solid indices buffer@>=
sb = mem_typed_alloc(m, uint32_t, t);
for (i = 0; i < num_subcuboids; ++i) {
    g->ctab[i].s = g->isb;
    for (j = 0; j < forest_of_solids.n; ++j) {
        s = forest_of_solids.s[j];
        if (no_intersection_bb(g->ctab[i].bb, s->bb)) continue;
        if (c == m) {
           sb = mem_typed_alloc(m, uint32_t, t);
           if (NULL == sb) goto exit_error;
	   c = 0;
	}
    	sb[c++] = j;
	++g->isb;
    }
    g->ctab[i].c = g->isb - g->ctab[i].s;
}

@ @<Finalise solid indices buffer by moving paged data to contiguous memory@>=
g->nsb = g->isb;
if (*t) {
    g->sb = mem_typed_alloc(g->nsb, uint32_t, mem_phase_two);
    if (NULL == g->sb) goto exit_error;
    @<Transfer the last page of the solid indices buffer@>;
    @<Transfer the remaining pages of the solid indices buffer@>;
}

@ Note that the memory pages are maintained inside the memory area as
a reverse linked list.

@<Transfer the last page of the solid indices buffer@>=
i = g->nsb - c; /* index in contiguous memory for the last page */
memcpy(&(g->sb[i]), (*t)->first, sizeof(uint32_t) * c);
*r = (*t)->next;
free((*t)->first);
*t = *r;
i -= m;

@ @<Transfer the remaining pages of the solid indices buffer@>=
while (*t) {
    *r = (*t)->next;
    memcpy(&(g->sb[i]), (*t)->first, sizeof(uint32_t) * m);
    free((*t)->first);
    *t = *r;
    i -= m;
}

@ @<Global functions@>=
void create_geotab(GeometryTable *g)
{
    uint32_t i;
    CSG_Node *s;
    @<Initialise the geometry table@>;
    @<Build tables and search trees for managing the subcuboids@>;
    fill_geotab_subcuboids_table(g);
    @<Fill in the primitives table, the solids table and postfix buffer@>;
    @<Check if there are stray node@>;
}

@ @<Initialise the geometry table@>=
g->ip = g->is = g->ic = g->ipb = g->isb = g->nsb = 0;
g->sw = sim_world;
g->nc = num_subcuboids;
g->l = div_subcuboids[0];
g->m = div_subcuboids[1];
g->n = div_subcuboids[2];
g->np = nodes_repo->stat[PRIMITIVE];
g->npb = nodes_repo->stat[PRIMITIVE] +
       nodes_repo->stat[UNION] +
       nodes_repo->stat[INTERSECTION] +
       nodes_repo->stat[DIFFERENCE];
g->ns = forest_of_solids.n;
g->p = mem_typed_alloc(g->np, struct primitives_table_item, mem_phase_two);
g->pb = mem_typed_alloc(g->npb, int32_t, mem_phase_two); 
g->s = mem_typed_alloc(g->ns, struct solids_table_item, mem_phase_two);

@ @<Fill in the primitives table, the solids table and postfix buffer@>=
for (i = 0; i < forest_of_solids.n; ++i) {
    s = forest_of_solids.s[i];
    @<Fill table entries for this solid@>;
}

@ @<Fill table entries for this solid@>=
if (NULL == s) continue;
g->s[g->is].s = g->ipb;
fill_geotab_csg_table(g, s);
g->s[g->is].c = g->ipb - g->s[g->is].s;
++(g->is);

@ @<Check if there are stray node@>=
i = g->npb - g->ipb;
if (i) {
   fprintf(stderr, "! There are %u stray nodes that are not in any of the solids:\n", i);
   for (i = 0; i < MAX_CSG_NODES; ++i) {
       s = nodes_repo->table[i];
       if (NULL == s || is_used(s)) continue;
       fprintf(stderr, "\t\"%s\" at line %u\n", s->name, get_line(s));
   }
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
	fprintf(f, "\nSubcuboids table:\n");
	for (i = 0; i < g->nc; ++i) fprintf(f, "%u %u\n", g->ctab[i].s, g->ctab[i].c);
	fprintf(f, "Solid indices buffer:\n");
	for (i = 0; i < g->nsb; ++i) fprintf(f, "%u ", g->sb[i]);
	fprintf(f, "\n");
	print_subcuboid_search_trees(stdout, geotab.ctree);
	print_neighbour_table(stdout, &geotab);
	print_subcuboids_table(stdout, &geotab);
}

@ @<Test geometry table generation@>=
{
        if (false == read_geometry("input.dat")) exit(1);
	print_geom_statistics(stdout);
	create_geotab(&geotab);
	print_geotab(stdout, &geotab);
}
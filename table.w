@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@** Preparing the shared tables.

All of the geometry information must be translated into tables before
they can be used during the particle simulations. The aim is to
simplify the data structures so that they are efficient, and also
allow using devices where management of complex data structures are
prohibitive. For instance, GPU memory management relies on the host
application; hence, it is complicated to manage data structures that
uses pointers.

\.{MCS} maintains all of the simplified data structures as tables
indside a container of tables, known as |geotab|. This is used as a
temporary container during the translation and simplification
process. Once it is ready, the compact form is transferred to a stage
two memory on the devices that will run the simulations.

To store information defining the forest of solids and the subcuboids
containing them, we could have opted for a straightforward approach
where each of the geometry components are stored separately. However,
most of the components in |geotab| are related to one another during
processing. For instance, searching for the solid which contains a
particle requires starting at the simulation world cuboid, and going
deeper through the subcuboids, and the CSG tree until we find a
containing primitive. Hence, to take advantage of these relationships,
the components are stored interlinked using a hierarchcal
structure. The tabular representation of this hierarchy is designed so
that the cache utilisation is efficient.

@ Assume that our simulation world is divided into four subcuboids as
shown below, and that it contains three silids $A$, $B$, and $C$. The
CSG tree for each of the solids are shown on the left-hand
side.

\bigskip

\centerline{\epsfig{file=figures/geometry-solids,scale=1}}

\bigskip

According to this diagram, the solid containment scenario is: (0:
B), (1: B, C), (2: A, B), (3: B). This is stored in a compact form
using the subcuboids table, which we saw in the previous section. This
assumes that the solids $A$, $B$ and $C$ have been assigned the
indices 0, 1, and 2 respectively.

@ We have already seen part of this hierarchical representation when we
discussed division of the simulation world into subcuboids, which is
reproduced in the following diagram. This captures the link:
Simulation world $\rightarrow$ Subcuboids Table $\rightarrow$ Solid
Indices Buffer.

\bigskip

\centerline{\epsfig{file=figures/geometry-subcuboids-table,scale=1}}

\bigskip

This hierarchical linkage can be expanded until we reach the lowest
level CSG primitives: Solid Indices $\rightarrow$ Solids Table
$\rightarrow$ Postfix Expression Buffer $\rightarrow$ Primitives
Table. We shall now discuss this extension.

@ To simulate a particle, we start by finding the subcuboid that
contains the particle. Then, we go to the row in the subcuboids table
and retrieve the start index and solids count. Using these values, we
retrieve the solids from the solid indices buffer.

For each of the solids indexed by values within the specified range in
the buffer, we check if that solid contains the particle. To do this,
we then use the {\it solids table}, which captures information
concerning the solids, i.e., the corresponding CSG tree.

@ A solids table is a compact form that stores the forest of CSG
trees. Instead of storing each of the solids separately, we store the
CSG trees as postfix expressions inside a common {\it postfix
expression buffer}. Now, as in the case with the subcuboids table, we
store the corresponding start indices and the expression length inside
the solids table. This is shown in the following example:

\bigskip

\centerline{\epsfig{file=figures/geometry-solids-table,scale=1}}

@<Type definitions@>=
struct solids_table_item {
    uint32_t s, c; /* start index and item count in postfix expression buffer */
    BoundingBox bb; /* bounding box for the solid */
};

@ The postfix expression buffer stores indices and operators of a CSG
tree. All of the negative integers are operators, where -1 represents
a boolean difference, -2 a boolean intersection, and -3 a boolean
union. All of the positive integers are indices to the {\it primitives
table}, which store information concerning the parameters specific to
the primitive instances.

\bigskip

\centerline{\epsfig{file=figures/geometry-primitives-table,scale=1}}

@<Type definitions@>=
struct primitives_table_item {
    Matrix a, i; /* accumulated affine and inverse transformations */
    Primitive p; /* primitive data */
};

@ Once we reach the primitives table, we can use the affine
transformation matrices (|a| or |i|) to transform the particle's
position vector, and test containment inside the primitive by using
the primitive's parameters available in |p|. This containment testing
inside primitives is carried out as we evaluate the boolean postfix
expression for a given solid.

For instance, if we are testing containment inside the solid $B$, we
will first evaluate containment inside the primitives using rows 4 and
5 in the primitives table, which should return boolean values, and
then calculate the boolean difference of the results.

@ The collection of tables that captures the hierarchical relationship
is stored inside the following data structure, which is
defined as the type |GeometryTable|.
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
    if (is_primitive(n)) {
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

@ A solid with a long postfix expression will take longer to evaluate,
compared to shorter expressions. Hence, we use |compare_solids(a,b)|
to sort the solids table in ascending order using the postfix expression
length as the comparison key.
@<Global functions@>=
static int compare_solids(const void *a, const void *b)
{
	struct solids_table_item *ap = (struct solids_table_item *) a;
	struct solids_table_item *bp = (struct solids_table_item *) b;
	return (ap->c - bp->c);
}

@ @<Global functions@>=
bool fill_geotab_subcuboids_table(GeometryTable *g)
{
    Area t, r; /* temporary paged solids indices buffer */
    uint32_t m = 3; /* maximum number of items per page */
    uint32_t c = 0; /* number of items in current page */
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
if (NULL == sb) goto exit_error;
qsort(g->s, g->ns, sizeof(struct solids_table_item), compare_solids);
for (i = 0; i < g->nc; ++i) {
    g->ctab[i].s = g->isb;
    for (j = 0; j < g->ns; ++j) {
        if (no_intersection_bb(g->ctab[i].bb, g->s[j].bb)) continue;
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
    @<Fill in the primitives table, the solids table and postfix buffer@>;
    @<Check if there are stray node@>;
    fill_geotab_subcuboids_table(g);
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
g->s[g->is].bb = s->bb;
g->s[g->is].s = g->ipb;
fill_geotab_csg_table(g, s);
g->s[g->is].c = g->ipb - g->s[g->is].s;
++(g->is);

@ @<Check if there are stray node@>=
i = g->npb - g->ipb;
if (i) {
   uint32_t j, p = 0;
   fprintf(stderr, "! There are %u stray nodes that are not in any of the solids:\n", i);
   for (j = 0; i && j < MAX_CSG_NODES; ++j) {
       s = nodes_repo->table[j];
       if (NULL == s || is_inuse(s)) continue;
       if (is_primitive(s)) {
           ++p;
       	   fprintf(stderr, "\tPrimitive \"%s\" at line %u\n", s->name, get_line(s));
       } else {
           fprintf(stderr, "\tOperator \"%s\" at line %u\n", s->name, get_line(s));
       }
       --i;
   }
   g->np -= p; /* correct the number of active primitives */
   g->npb = g->ipb; /* correct the number of items in postfix buffer */
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
        if (false == read_geometry("test/test_gpu_table.data")) exit(1);
	print_geom_statistics(stdout);
	create_geotab(&geotab);
	print_geotab(stdout, &geotab);
}

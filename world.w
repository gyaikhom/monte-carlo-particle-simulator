@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@*1 The simulation world.
The simulation world is defined by a {\sl cuboid},
which is composed of solids. Only particles inside this cuboid are
tracked during a simulation. To exploit parallelism during
simulations, {\tt MCS} logically divides this cuboid into smaller
cuboids of equal shape and size, each referred to as a {\sl
subcuboid}. A set of subcuboids decomposes the simulation world into
disjointed volumes where particles can be tracked independently of
other particles in other subcuboids.

The simulation world must be axis-aligned with the world coordinate
frame. It is decomposed into subcuboids $S = \{s_{ijk} : 1 \le i \le l,
1 \le j \le m, 1 \le k \le n\}$, where each of the cuboid's three
dimensions has been divided into $l$, $m$ and $n$ equal parts along
the $x$, $y$ and $z$ axes, respectively. This division can be carried
out in any combinations of $l$, $m$ and $n$, as long as they are all
whole numbers. The other conditions are that all of the subcuboids
resulting from a division must all have the same size and shape
(i.e., must be a cuboid themselves), must be disjointed, and their
union must exactly produce the original cuboid.

For instance, in the following, the two divisions on the left are
valid, whereas the right is invalid:

\bigskip

\centerline{\epsfig{file=figures/subdivision,scale=0.3}\hfil \epsfig{file=figures/nonuniform,scale=0.3}}

\bigskip

We do not use space-partitioning data structures such as the {\sl
octree} or the {\sl $kd$-tree} to divide the simulation world because
it is important that the subcuboids are of the same size and
shape. With unequal partitions, we must carry out a lot more
calculations in order to determine to which neighbouring subcuboid a
particle has escaped. Even worse, this computational overhead is
nonuniformly spread differing from one particle to the next, depending
on their positions. If we were simulating a few thousand particles,
the computational overhead may be considered acceptable. However,
since {\tt MCS} is expected to simulate millions, or even billions of
particles, the overhead becomes significant and therefore
unacceptable. By using subcuboids of the same size and shape, we can 
avoid all of this overhead by using an efficient lookup-table. We
shall discuss this in the following sections. Note, however, that we
could still use an octree or a $kd$-tree to accelerate the search for
a solid inside a subcuboid.

@ After the simulation world has been divided, all of the subcuboids
are stored in a table, with entries as shown below. We shall refer to
this as the {\sl subcuboids table}. The solid indices buffer stores
indices of solids that are inside a given subcuboid. Hence, for a
given subcuboid, the index |s| gives the starting location of the
solids list within this buffer, and |c| gives the number of items in
that list. For instance, the subcuboid with index 2 contains two
solids, whose solid indices are 0 and 1. What these solid indices mean
will become clear when we discuss the solids table in later sections.

\bigskip

\centerline{\epsfig{file=figures/geometry-subcuboids-table,scale=1}}

@<Type definitions@>=
struct subcuboids_table_item {
    uint32_t s, c;  /* start index and item count in solid indices buffer */
    BoundingBox bb;
};

@ @<Subcuboids related data structures inside the geometry table@>=
struct subcuboids_table_item *ctab;

@ Function |build_subcuboids_table(bb,l,m,n)| builds the subcuboids
table by filling in the bounding box information for each of the
subcuboids. The dimension of each subcuboid is determined from the
size of the original cuboid |bb|, and the number of equal divisions
in each of the dimensions $x$, $y$ and $z$ as given respectively by
|l|, |m| and |n|.
@<Global functions@>=
typedef struct geomtab_struct GeometryTable; /* forward-declaration */
bool build_subcuboids_table(GeometryTable *g)
{
	uint32_t i, j, k, c, t;
	double dx, dy, dz, x[2], y[2], z;
	g->ctab = mem_typed_alloc(g->nc, struct subcuboids_table_item,
	mem_phase_two);
	if (NULL == g->ctab) return false;
	t = g->m * g->n;
	dx = (g->sw.u[0] - g->sw.l[0]) / g->l; /* determine subcuboid size using upper and lower bounds */
	dy = (g->sw.u[1] - g->sw.l[1]) / g->m; /* of the cuboid that represents the simulation world */
	dz = (g->sw.u[2] - g->sw.l[2]) / g->n;
        x[0] = g->sw.l[0];
	for (i = 0; i < g->l; ++i) {
            x[1] = x[0] + dx;
            y[0] = g->sw.l[1];
	    for (j = 0; j < g->m; ++j) {
	        y[1] = y[0] + dy;
                z = g->sw.l[2];
	        for (k = 0; k < g->n; ++k) {
	            c = i * t + j * g->n + k;
	            g->ctab[c].bb.l[0] = x[0];
		    g->ctab[c].bb.l[1] = y[0];
		    g->ctab[c].bb.l[2] = z;
	            z += dz;
	            g->ctab[c].bb.u[0] = x[1];
	            g->ctab[c].bb.u[1] = y[1];
                    g->ctab[c].bb.u[2] = z;
                }
	        y[0] += dy;
            }
            x[0] += dx;
     	}
	return true;
}

@ Function |print_subcuboids_table(f,g)| prints the subcuboids table
from the geometry table |g| to the I/O stream pointed to by |f|.
@<Global functions@>=
void print_subcuboids_table(FILE *f, GeometryTable *g)
{
	uint32_t i;
	fprintf(f, "Subcuboids Table:\n");
	for (i = 0; i < g->nc; ++i)
	    fprintf(f, "%3u [%10.3lf, %10.3lf, %10.3lf : %10.3lf, %10.3lf, %10.3lf]\n",
	    i, g->ctab[i].bb.l[0], g->ctab[i].bb.l[1],
	    g->ctab[i].bb.l[2], g->ctab[i].bb.u[0],
	    g->ctab[i].bb.u[1], g->ctab[i].bb.u[2]);
}


@*2 Finding the subcuboid containing a particle.
Once the simulation world has been divided into subcuboids, how do we
find the subcuboid that contains a given particle? This question 
arises when we are generating primary particles using a particle
gun. The only way to determine the containing subcuboid is to search
throughout the simulation world by going through all of the
subcuboids. Again, for a few thousand particles the computational cost
for an exhaustive search would be acceptable; however, since {\tt MCS}
is expected to simulate millions, we must devise a better strategy.

First of all, we only search for particles that are inside the
simulation world. Then, we take advantage of the following four design
conditions:

$$\vcenter{\halign{# & #\hfil \cr
1) & all of the subcuboids are of the same size and shape,\cr
2) & their union exactly defines the simulation world,\cr
3) & no two subcuboids intersect, and\cr
4) & the division into subcuboids is immutable throughout the entire simulation.\cr
}}$$

This means that, we can decompose a single search for the subcuboid in
three-dimensional space into three separate, but faster, searches in
one-dimensional space. At the end of the three searches, we simply
combine the separate results to produce the index of the subcuboid
that contains the particle.

@<Subcuboids related data structures inside the geometry table@>=
double *ctree; /* stores three binary search trees */

@ Every node in the binary search tree stores the upper bound of each
subcuboid in the selected axis, except for the upper bound that
coincides with the upper bound of the simulation world. Hence, if we
divide the simulation world into $n$ subcuboids in a given dimension,
we require $n - 1$ nodes for the search tree in that
dimension. Furthermore, since the division into subcuboids is
immutable for the entire simulation, we could improve the search by
using a {\sl complete binary tree}, instead of using an arbitrary
binary search tree. With a complete binary tree the number of
comparisons per search with $n$ subcuboids per dimension is bound by
$O(\lceil \lg n \rceil)$. Finally, we can take advantage of the fact
that a complete binary search tree can be stored using a
one-dimensional array. This not only reduces the space requirements
for storing the binary tree, but also simplifies the tree traversal
mechanism: Instead of using pointers, we use {\sl bit shift} operators.

@ Function |build_complete_tree(t,n,i)| recursively builds a complete
binary search tree |t| so that, during an {\sl inorder} traversal, the
values in the internal nodes will produce the sequence $s = \{s_i, 0
\le i < n\}$, where $s_0 = \delta$ and $s_j = s_{j-1} + \delta, 0 < j <
n$. The increment $\delta$ is the dimension of each subcuboid in the
selected axes. The leaf nodes, on the other hand, store partial
indices of the subcuboids that contain the value range defined by the
internal nodes in the chosen axis. Later on, partial indices from the
three axes can be combined to produce a complete subcuboid index.
@<Global functions@>=
static double build_complete_tree_val = 0.0; /* upper bound stored in
internal node */
static double build_complete_tree_inc = 0.0; /* size of the subcuboid */
static int build_complete_leaf_idx = 0; /* partial subcuboid index
stored in leaf */
void build_complete_tree(double *t, unsigned long n, unsigned long i)
{
	unsigned long j;
	if (i < n) {
		j = i << 1;
		build_complete_tree(t, n, j);
		build_complete_tree_val += build_complete_tree_inc;
		t[i] = build_complete_tree_val;
		build_complete_tree(t, n, j + 1);
	} else *(unsigned long *) &t[i] = build_complete_leaf_idx++;
}

@ Function |build_subcuboid_binary_search_tree(t,n,u,l)| builds a complete
binary search tree |t| where the upper and lower bounds of the
simulation world are respectively |u| and |l| units in the selected
axis. On the selected axis, the simulation world has been divided
into |n| equal subcuboids.
@<Global functions@>=
void build_subcuboid_binary_search_tree(double *t, unsigned long n, double u, double l)
{
	build_complete_tree_val = 0.0;
	build_complete_tree_inc = (u - l) / n;
	build_complete_leaf_idx = 0;
	build_complete_tree(t, n, 1);
}

@ Function |build_subcuboid_search_trees(g)| builds three complete binary
search trees using the simulation world data inside the geometry table
|g|. The simulation world is bound by the bounding box |bb|, and is
divided into |l|, |m| and |n| equal parts along the $x$, $y$ and $z$
axes. There will be a total of $|l| \times |m| \times |n|$ subcuboids,
which this function assumes is less than |MAX_SUBCUBOIDS|.

@<Global functions@>=
bool build_subcuboid_search_trees(GeometryTable *g)
{
	double *x, *y, *z;
	uint32_t nx, ny, nz;
	@<Allocate memory for the three subcuboid search trees@>;
	build_subcuboid_binary_search_tree(x, g->l, g->sw.u[0], g->sw.l[0]);
	build_subcuboid_binary_search_tree(y, g->m, g->sw.u[1], g->sw.l[1]);
	build_subcuboid_binary_search_tree(z, g->n, g->sw.u[2], g->sw.l[2]);
	@<Finalise the subcuboid search trees@>;
	return true;
}

@ This allocates a single array to hold all of the three trees. Each
search tree with |n| subcuboids only requires |n - 1| internal nodes
and |n| leaf nodes. However, we allocate one extra element per tree
because the root only begins at the second element (index 1). Instead
of wasting this first element, we use it to store the number of
nodes. Hence, to store a binary subcuboid search tree with $k$
subcuboids, we allocate an array with $2k$ elements. The same is
calculated and allocated for each of the other two axes.
@<Allocate memory for the three subcuboid search trees@>=
nx = 2 * g->l;
ny = 2 * g->m;
nz = 2 * g->n;
g->ctree = mem_typed_alloc(nx + ny + nz, double, mem_phase_two);
if (NULL == g->ctree) return false;
x = g->ctree;
y = x + nx;
z = y + ny;

@ Don't forget to store the number of nodes as the first element of
the tree array.
@<Finalise the subcuboid search trees@>=
*(unsigned long *)x = nx;
*(unsigned long *)y = ny;
*(unsigned long *)z = nz;

@ Function |print_subcuboid_search_trees(f,t)| prints the three
complete binary search trees in |t| to the I/O stream that is pointed
to by |f|.
@<Global functions@>=
void print_subcuboid_search_trees(FILE *f, double *t)
{
	unsigned long i, j, k, size;
	for (i = 0; i < 3; ++i) {
	    size = *(unsigned long *) t;
	    k = size / 2;
	    fprintf(f, "tree with %lu nodes:\n\tinternal: ", size);
	    for (j = 1; j < k; ++j) fprintf(f, "%lf ", t[j]);
    	    fprintf(f, "\n\tleaves: ");
	    for (; j < size; ++j) fprintf(f, "%lu ", *(unsigned long *) &t[j]);
	    fprintf(f, "\n");
   	    t += size; /* move to next tree array */
	}
}

@ Function |find_subcuboid(t,v)| finds the subcuboid that contains a
point with position vector |v| by searching the three binary search
trees pointed to by |t|. 

Assume that the simulation world is divided into $l \times m
\times n$ subcuboids. After the three one-dimensional searches, let
$i$, $j$ and $k$ represent the search results along the $x$, $y$ and
$z$ axes respectively. The index of the subcuboid is then given by
$(i \times m \times n + j \times n + k)$.

@<Global functions@>=
unsigned long find_subcuboid(double *t, Vector v)
{
	unsigned long i, j, c[3], size[3];
	for (i = 0; i < 3; ++i) @<Search for subcuboid in the current tree@>;
	return (c[0] * size[1] * size[2] + c[1] * size[2] + c[2]);
}

@ We start at the root and travel left or right, depending on the
value at the node, until we reach a leaf, where we have already stored
the index of the containing subcuboid.
@<Search for subcuboid in the current tree@>=
{
size[i] = *(unsigned long *) t / 2;
j = 1;
while (j < size[i]) {
    if (v[i] < t[j]) j <<= 1; /* left subtree */
    else j = (j << 1) + 1; /* right subtree */
}
c[i] = *(unsigned long *) &t[j]; /* index found */
t += 2 * size[i]; /* move to the next binary search tree */
}

@*2 Relationship between subcuboids.
Once inside a simulation, we must determine at each step of
the particle's trajectory which subcuboid currently encloses the
particle. This is important because in order to apply physics
processes to a particle, we must first determine the physical
properties of the material with which the particle will interact in
that step. To determine the material properties, however, we must
know which solid is currently enclosing the particle.

For primary particles that are only starting their trajectory, we
must use the function |find_subcuboid(t,v)| to find the containing
subcuboid. However, for particles already in simulation, we can devise
a better strategy which utilises the current particle position
relative to its position in the previous trajectory step. To do this,
assume that the planes that divide the simulation world also divide
the void outside the simulation world cuboid. Then every subcuboid is
surrounded by 26 neighbouring subcuboids, where each neighbour is
either part of the simulation world or belongs to the void. This is
shown below:

\bigskip

\centerline{\epsfig{file=figures/subcuboids,scale=0.35}}

\bigskip

For a given subcuboid, we can define all of its 26 neighbours using
the six faces of the subcuboid. Let $f_i, 0 \le i < 6$ represent the
six faces so that the intervals $[f_0, f_1]$, $[f_2, f_3]$, and
$[f_4, f_5]$ respectively define the region enclosed by the subcuboid
along the $x$, $y$ and $z$ axes. Now, for every point $(x, y, z)$ in
the three-dimensional space, we can define a bit-field
$s = \langle s_5, s_4, s_3, s_2, s_1, s_0 \rangle$, where $s_0$ is the
least significant bit. In the following sections, we shall refer to
this bit-field as the {\sl $s$-field}.@^$s$-field@>

All of the bits in $s$ are zero, except in the following cases:

$s_0 = 1$ if $x < f_0$, $s_1 = 1$ if $x > f_1$,

$s_2 = 1$ if $y < f_2$, $s_3 = 1$ if $y > f_3$,

$s_4 = 1$ if $z < f_4$, and $s_4 = 1$ if $z > f_5$.

For instance, points inside the subcuboid will have $s =
\langle000000\rangle$;  whereas, points inside the top neighbour will
have $s = \langle001000\rangle$. As we can see, the bit pairs $(s_0,
s_1)$, $(s_2, s_3)$, and $(s_4, s_5)$ will never be both
nonzero. Furthermore, we can represent the faces using a bounding box,
where the points $(f_0, f_2, f_4)$ and $(f_1, f_3, f_5)$ are
respectively the lower and upper bounds.

@ Function |update_sfield(s,bb,v)| updates the $s$-field |s| for a
particle currently at the position vector |v| using the subcuboid
boundary defined by the bounding box |bb|. 
@<Global functions@>=
void update_sfield(uint8_t *s, BoundingBox *bb, Vector v)
{
	*s = 0x0;
	if (v[0] < bb->l[0]) *s |= 0x1; /* $s_0$ */
	else if (v[0] > bb->u[0]) *s |= 0x2; /* $s_1$ */
	if (v[1] < bb->l[1]) *s |= 0x4; /* $s_2$ */
	else if (v[1] > bb->u[1]) *s |= 0x8; /* $s_3$ */
	if (v[2] < bb->l[2]) *s |= 0x16; /* $s_4$ */
	else if (v[2] > bb->u[2]) *s |= 0x32; /* $s_5$ */
}

@ While a particle is being tracked, we must check if it continues to
exists inside the same subcuboid, or if it has escaped to a
neighbouring subcuboid. We shall use the $s$-field to efficiently
determine the relevant case for every particle after applying a
trajectory step.

Every particle in simulation maintains an $s$-field, which gives its
location relative to the current subcuboid enclosing it. Before
applying a trajectory step for the first time, the $s$-field is
initialised once to $\langle000000\rangle$. After applying 
a trajectory step, which could have changed the particle's location, we
update the $s$-field using |update_sfield(s,bb,v)|. After subsequent
application of a trajectory step, we simply determine the next
effective subcuboid using the particle's current $s$-field, its
containing subcuboid, and a neighbourhood table.

The subcuboid neighbourhood table is an $m \times 26$ two-dimensional
array where each row represents one of the $m$ subcuboids, and every
element of the row points to one of the 26 neighbouring
subcuboids. During lookup, the relevant $s$-field determines the
appropriate column to choose. For instance, if a particle escapes to
one of the 26 neighbouring subcuboids, say the top neighbour, it will
have the $s$-field $\langle001000\rangle$. To find the subcuboid
containing the particle in the next step, we use this as the column
during table lookup. So, if the particle is inside the subcuboid with
index 10 in the current step, and the $s$-field is
$\langle001000\rangle$, the index of the containing subcuboid for the
next step is given by the table entry at row 10 and column 8.

@ If we were to allocate the subcuboid neighbourhood lookup-table
using direct-mapping, the table will waste $16 \times
|MAX_SUBCUBOIDS|$ table elements. This is because, in order to make
every valid $s$ value indexable, we must allocate |MAX_SFIELD| table
elements per row. However, out of this row only 26 elements actually
contain a valid index to one of the neighbouring subcuboids. This is a
consequence of the fact that the bit pairs $(s_0, s_1)$, $(s_2, s_3)$,
and $(s_4, s_5)$ will never be both nonzero. Direct mapping,
therefore, is space inefficient.

\bigskip

\centerline{\epsfig{file=figures/neighbour-direct,scale=1}}

\bigskip

If we were to use a two-tiered mapping, however, we only waste space
equivalent to 16 table elements. This works by allocating an index
lookup-table, which maps the first |MAX_SFIELD| $s$ values to an index
that points to one of the 26 valid neighbouring subcuboids.

\bigskip

\centerline{\epsfig{file=figures/neighbour-tiered,scale=1}}

@ We now allocate a two-tiered neighbourhood table.
@d MAX_SUBCUBOIDS 1024
@d NUM_NEIGHBOURS 26 /* number of cuboid neighbours */
@d MAX_SFIELD 42 /* maximum $s$-field value: $\langle101010\rangle$ */
@<Subcuboids related data structures inside the geometry table@>=
int8_t iltab[MAX_SFIELD + 1]; /* index lookup table */
uint32_t **ntab; /* each row containing |NUM_NEIGHBOURS| entries */

@ Function |get_neighbour(g,s,i)| returns the next effective subcuboid
for a particle, where it was previously in subcuboid |i| and the most
recent application of a trajectory step updated the $s$-field to
|s|. It uses the neighbourhood table in |g|.

@d subcuboid_lookup(g,i,s) ((g)->ntab[(i)][(g)->iltab[s]])
@<Global functions@>=
uint32_t get_neighbour(GeometryTable *g, uint32_t i, uint8_t s)
{
	if (s) return subcuboid_lookup(g,i,s);
	return i; /* particle continues to exist in the same subcuboid */
}

@ Function |build_neighbour_table(g)| builds the subcuboid
neighbourhood lookup-table inside the geometry table |g|. The cuboid
which represents the simulation world has been divided into |l|, |m|
and |n| equal parts along the $x$, $y$ and $z$ axes. There will be a
total of $|l| \times |m| \times |n|$ subcuboids, which this function
assumes is less than |MAX_SUBCUBOIDS|.

@d subcuboid_assign(g,r,s,v) (g)->ntab[(r)][(g)->iltab[(int)(s)]] = (v)
@<Global functions@>=
bool build_neighbour_table(GeometryTable *g)
{
    uint32_t i, j, k, x, y, z, t = g->m * g->n;
    uint32_t r; /* row for subcuboid currently being filled in */
    uint8_t s, xb, yb, zb; /* $s$-field and extracted bit pairs */
    int8_t q[] = {-1,  0,  1, -1,  2,  3,  4, -1,  5,  6,
     7, -1, -1, -1, -1, -1,  8,  9, 10, -1,
    11, 12, 13, -1, 14, 15, 16, -1, -1, -1,
    -1, -1, 17, 18, 19, -1, 20, 21, 22, -1,
    23, 24, 25}; /* index lookup table data */
    g->ntab = mem_typed_alloc2d(g->nc, NUM_NEIGHBOURS, uint32_t,
    mem_phase_two);
    if (NULL == g->ntab) return false;
    memcpy(g->iltab, q, MAX_SFIELD + 1);
    for (i = 0; i < g->l; ++i)
        for (j = 0; j < g->m; ++j)
	    for (k = 0; k < g->n; ++k) { /* set neighbours */
                r = i * t + j * g->n + k; /* row index for the current subcuboid */
		for (s = 1; s <= MAX_SFIELD; ++s) {
    		    x = i; y = j; z = k;
    		    if (((xb = s & 0x3) == 0x3) ||
        	        ((yb = s & 0xC) == 0xC) ||
			((zb = s & 0x30) == 0x30)) continue;
                    else @<Calculate the neighbour using the axes bit-pairs@>;
                }
            }
    return true;
}

@ @<Calculate the neighbour using the axes bit-pairs@>=
{
@<Adjust index of the neighbour cuboid along the $x$-axis@>;
@<Adjust index of the neighbour cuboid along the $y$-axis@>;
@<Adjust index of the neighbour cuboid along the $z$-axis@>;
subcuboid_assign(g, r, s, x * t + y * g->n + z);
}

@ @<Adjust index of the neighbour cuboid along the $x$-axis@>=
if (xb == 0x1) {
   if (--x < 0) @<Neighbour subcuboid is outside the simulation world@>;
} else if (xb == 0x2) {
   if (++x == g->l) @<Neighbour subcuboid is outside the simulation world@>;
}

@ @<Adjust index of the neighbour cuboid along the $y$-axis@>=
if (yb == 0x4) {
   if (--y < 0) @<Neighbour subcuboid is outside the simulation world@>;
} else if (yb == 0x8) {
   if (++y == g->m) @<Neighbour subcuboid is outside the simulation world@>;
}

@ @<Adjust index of the neighbour cuboid along the $z$-axis@>=
if (zb == 0x10) {
   if (--z < 0) @<Neighbour subcuboid is outside the simulation world@>;
} else if (zb == 0x20) {
   if (++z == g->n) @<Neighbour subcuboid is outside the simulation world@>;
}

@ All of the subcuboids that belong to the void outside the simulation
world are given the same subcuboid index |OUTSIDE_WORLD|.
@d OUTSIDE_WORLD (MAX_SUBCUBOIDS + 1)
@<Neighbour subcuboid is outside the simulation world@>=
{
    subcuboid_assign(g,r,s, OUTSIDE_WORLD);
    continue;
}

@ @<Build tables and search trees for managing the subcuboids@>=
build_subcuboid_search_trees(g);
build_neighbour_table(g);
build_subcuboids_table(g);

@ Function |print_neighbour_table(f,g)| prints the subcuboid
neighbourhood table in |g| to the I/O stream pointed to by |f|.
@<Global functions@>=
void print_neighbour_table(FILE *f, GeometryTable *g)
{
	uint32_t i, j, c, k = MAX_SUBCUBOIDS + 1;
	fprintf(f, "Subcuboid neighbourhood table:\n");
	for (i = 0; i < g->nc; ++i) {
	    for (j = 0; j < NUM_NEIGHBOURS; ++j) {
	        c = g->ntab[i][j];
	    	if (c < k) fprintf(f, "%3d ", c);
		else fprintf(f, " .  ");
	    }
	    fprintf(f, "\n");
	}	
}

@ The usage of the neighbourhood table makes the assumption that
particles escaping a subcuboid will be found in one of the 26
surrounding subcuboids. However, this assumption may not apply if
the distance that a particle travels after a step application exceeds
the dimension of the subcuboid in one of the three axes. In these
cases, we must fall-back to finding the subcuboid using
|find_subcuboid(t,v)|. In most cases, the subcuboid dimensions will be
larger than the distance travelled, however, it is important to note
the exception. After every step, we must check this distance.

@ The following code segment tests the functionalities provided by
this sections.

@<Test subcuboid functionalities@>=
{
	uint32_t i, j, k, e, r;
	double dx, dy, dz;
	Vector v;
	do {
	       printf("Give upper bound (x, y, z) of the simulation world:\n");
	       scanf("%lf %lf %lf", &geotab.sw.u[0], &geotab.sw.u[1], &geotab.sw.u[2]);
	       printf("Give lower bound (x, y, z) of the simulation world:\n");
	       scanf("%lf %lf %lf", &geotab.sw.l[0], &geotab.sw.l[1], &geotab.sw.l[2]);
	} while (geotab.sw.u[0] < geotab.sw.l[0] ||
	         geotab.sw.u[1] < geotab.sw.l[1] ||
	         geotab.sw.u[2] < geotab.sw.l[2]);
        do {
	       printf("Give number of divisions on the x, y, and z axes:\n");
	       scanf("%u %u %u", &geotab.l, &geotab.m, &geotab.n);
	       geotab.nc = geotab.l * geotab.m * geotab.n;
	} while (geotab.nc > MAX_SUBCUBOIDS);	
	build_subcuboid_search_trees(&geotab);
	print_subcuboid_search_trees(stdout, geotab.ctree);
	@<Search subcuboid for all points in each of the subcuboids@>;

	build_neighbour_table(&geotab);
	build_subcuboids_table(&geotab);
	print_neighbour_table(stdout, &geotab);
	print_subcuboids_table(stdout, &geotab);
}

@ @<Search subcuboid for all points in each of the subcuboids@>=
dx = (geotab.sw.u[0] - geotab.sw.l[0]) / (double) geotab.l;
dy = (geotab.sw.u[1] - geotab.sw.l[1]) / (double) geotab.m;
dz = (geotab.sw.u[2] - geotab.sw.l[2]) / (double) geotab.n;
v[0] = geotab.sw.l[0];
for (i = 0; i < geotab.l; ++i, v[0] += dx) {
    v[1] = geotab.sw.l[1];
    for (j = 0; j < geotab.m; ++j, v[1] += dy) {
        v[2] = geotab.sw.l[2];
	for (k = 0; k < geotab.n; ++k, v[2] += dz) {
	    e = i * geotab.m * geotab.n + j * geotab.n + k;
	    r = find_subcuboid(geotab.ctree, v);
	    if (e != r) goto test_subcuboid_search_failed;
        }
    }
}
printf("Test was a success...\n");
goto test_subcuboid_search_done;
test_subcuboid_search_failed:
       printf("Failure at (%u, %u, %u): %u instead of %u for (%lf, %lf, %lf)\n",
       i, j, k, r, e, v[0], v[1], v[2]);
test_subcuboid_search_done:

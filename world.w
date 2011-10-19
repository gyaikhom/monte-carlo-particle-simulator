@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

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

For instance, the following are valid divisions of
the simulation world:

\bigskip

\centerline{\epsfig{file=figures/subdivision,scale=0.3}}

\bigskip

We do not use space-partitioning data structures such as the {\sl
octree} or the {\sl $kd$-tree} to divide the simulation world because
it is important that the subcuboids are of the same size and
shape. With unequal partitions, we must carry out a lot more
calculations in order to determine to which neighbouring subcuboid a
particle has escaped. Even worse, this computational overhead is
nonuniformly spread differing from one particle to the next, depending
on their positions. The following figure demonstrates.

\bigskip

\centerline{\epsfig{file=figures/nonuniform,scale=0.3}}

\bigskip

If we were simulating a few thousand particles, the computational
 overhead may be considered acceptable. However, since {\tt
MCS} is expected to simulate millions, or even billions of
particles, the overhead becomes significant and therefore
unacceptable. By using subcuboids of the same size and shape, we can
avoid all of this overhead by using an efficient lookup-table. We
shall discuss this in the following sections.

To accelerate particle tracking inside a subcuboid, however, we could
still take advantage of an octree or a $kd$-tree to accelerate the
search for a solid inside a subcuboid.

@ After the simulation world has been divided, all of the subcuboids
are stored in the array |subcuboids|. We shall refer to this as the
{\sl subcuboids table}.
@<Type definitions@>=
struct {
    BoundingBox bb;
} subcuboids[MAX_SUBCUBOIDS]; /* subcuboids table */
uint32_t num_subcuboids = 0;

@ Function |build_subcuboids_table(bb,l,m,n)| builds the subcuboids
table by filling in the bounding box information for each of the
subcuboids. The dimension of each subcuboid is determined from the
size of the original cuboid, and the number of equal divisions carried
out in each of the dimensions.
@<Global functions@>=
void build_subcuboids_table(BoundingBox *bb, uint32_t l, uint32_t m, uint32_t n)
{
	uint32_t i, j, k, c, t = m * n;
	double dx, dy, dz, x[2], y[2], z;
	dx = (bb->u[0] - bb->l[0]) / l;
	dy = (bb->u[1] - bb->l[1]) / m;
	dz = (bb->u[2] - bb->l[2]) / n;
        x[0] = bb->l[0];
	for (i = 0; i < l; ++i) {
            x[1] = x[0] + dx;
            y[0] = bb->l[1];
	    for (j = 0; j < m; ++j) {
	        y[1] = y[0] + dy;
                z = bb->l[2];
	        for (k = 0; k < n; ++k) {
	            c = i * t + j * n + k;
	            subcuboids[c].bb.l[0] = x[0];
		    subcuboids[c].bb.l[1] = y[0];
		    subcuboids[c].bb.l[2] = z;
	            z += dz;
	            subcuboids[c].bb.u[0] = x[1];
	            subcuboids[c].bb.u[1] = y[1];
                    subcuboids[c].bb.u[2] = z;
                }
	        y[0] += dy;
            }
            x[0] += dx;
     	}
	num_subcuboids = t * l;
}

@ Function |print_subcuboids_table(f)| prints the subcuboids table to
the I/O stream pointed to by |f|.
@<Global functions@>=
void print_subcuboids_table(FILE *f)
{
	uint32_t i;
	fprintf(f, "Subcuboids Table:\n");
	for (i = 0; i < num_subcuboids; ++i)
	    fprintf(f, "%3u [%8.3lf, %8.3lf, %8.3lf : %8.3lf, %8.3lf, %8.3lf]\n",
	    i, subcuboids[i].bb.l[0], subcuboids[i].bb.l[1],
	    subcuboids[i].bb.l[2], subcuboids[i].bb.u[0],
	    subcuboids[i].bb.u[1], subcuboids[i].bb.u[2]);
}


@*2 Finding the subcuboid containing a particle.
Once the simulation world has been divided into subcuboids, how do we
find the subcuboid that contains a given particle? This question 
arises when we are generating primary particles using a particle
gun. The only way to determine the containing subcuboid is to search
throughout the simulation world by going through all of the
subcuboids. Again, for a few thousand particles the computational cost
for an exhaustive search could be acceptable; however, since {\tt MCS}
is expected to simulate millions, we must devise a better strategy.

First of all, we only search for particles that are inside the
simulation world. Then, we take advantage of the following four design
conditions:

$$\vcenter{\halign{# & #\hfil \cr
1) & all of the subcuboids are of the same size and shape,\cr
2) & their union exactly defines the simulation world,\cr
3) & no two subcuboids intersect, and\cr
4) & the division into subcuboids is immutable for the entire simulation.\cr
}}$$

This means that, we can decompose a single search for the subcuboid in
three-dimensional space into three separate, but faster, searches in
one-dimensional space. At the end of the three searches, we simply
combine the results to give the index of the subcuboid.

@<Global variables@>=
double *subcuboid_search_tree; /* stores three binary search trees */

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
mechanism: Instead of using pointers, we use {\sl binary shift} operators.

@ Function |build_complete_tree(t,n,i)| recursively builds a complete
binary search tree |t| so that, during an {\sl inorder} traversal, the
values in the internal nodes will produce the sequence $s = \{s_i, 0
\le i < n\}$, where $s_0 = \delta$ and $s_j = s_{j-1} + \delta, 0 < j <
n$. The increment $\delta$ is the dimension of each subcuboid in the
selected axes. The leaf nodes, on the other hand, store indices of the
subcuboids that contains the value range defined by the internal nodes.
@<Global functions@>=
static double build_complete_tree_val = 0.0;
static double build_complete_tree_inc = 0.0;
static int build_complete_leaf_idx = 0;
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

@ Function |build_subcuboid_search_tree(t,n,u,l)| builds a complete
binary search tree |t| where the upper and lower bounds of the
simulation world are respectively |u| and |l| units in the selected
axis. In the selected dimension, the simulation world has been divided
into |n| equal subcuboids.
@<Global functions@>=
void build_subcuboid_search_tree(double *t, unsigned long n, double u, double l)
{
	build_complete_tree_val = 0.0;
	build_complete_tree_inc = (u - l) / n;
	build_complete_leaf_idx = 0;
	build_complete_tree(t, n, 1);
}

@ Function |build_subcuboid_trees(bb, l,m,n)| builds three complete
binary search trees. The simulation world is bound by the bounding box
|bb|, and is divided into |l|, |m| and |n| equal parts along the $x$,
$y$ and $z$ axes. There will be a total of $|l| \times |m| \times |n|$
subcuboids, which this function assumes is less than
|MAX_SUBCUBOIDS|.

@<Global functions@>=
bool build_subcuboid_trees(BoundingBox *bb, unsigned long l, unsigned long m, unsigned long n)
{
	double *x, *y, *z;
	@<Allocate memory for the three cuboid search trees@>;
	build_subcuboid_search_tree(x, l, bb->u[0], bb->l[0]);
	build_subcuboid_search_tree(y, m, bb->u[1], bb->l[1]);
	build_subcuboid_search_tree(z, n, bb->u[2], bb->l[2]);
	@<Finalise the subcuboid search trees@>;
	return true;
}

@ This allocates a single array to hold all of the three trees. Each
search tree with |n| subcuboids only requires |n - 1| internal nodes
and |n| leaf nodes. However, we allocate one extra element per tree
because the root only begins at the second element (index 1). Instead
of wasting this first element, we use it to store the number of nodes.
@<Allocate memory for the three cuboid search trees@>=
subcuboid_search_tree = mem_typed_alloc(2 * (l + m + n), double, mem_p);
if (NULL == subcuboid_search_tree) return false;
x = subcuboid_search_tree;
y = x + 2 * l;
z = y + 2 * m;

@ Don't forget to store the number of nodes as the first element of
the tree array.
@<Finalise the subcuboid search trees@>=
*(unsigned long *)x = 2 * l;
*(unsigned long *)y = 2 * m;
*(unsigned long *)z = 2 * n;

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
properties of the material with which the particle will interacting at
that given step. However, to determine the materal properties, we must
know which solid is currently enclosing the particle. We therefore
start at the subcuboid.

For primary particles, we must use the function
|find_subcuboid(t,v)|to find the subcuboid; however, for particles in
simulation, we can devise a better strategy. Assume that the planes
that divide the simulation world also divides the void outside the
simulation. Then every subcuboid is surrounded by 26 other 
subcuboids, where each is either part of the simulation world or
belongs to the void, as shown below:

\bigskip

\centerline{\epsfig{file=figures/subcuboids,scale=0.4}}

\bigskip

For a given subcuboid, we can define all of its 26 neighbours using
the six faces of the subcuboid. Let $f_i, 0 \le i < 6$ represent the
six faces so that the intervals $[f_0, f_1]$, $[f_2, f_3]$, and
$[f_4, f_5]$ respectively define the region inside the subcuboid
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
exists inside the same subcuboid, or if it has exited to a
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

The neighbourhood table is an $m \times 26$ two-dimensional array where each
row represents one of the $m$ subcuboids, and every element of the row
points to one of the 26 neighbouring subcuboids. Which column index is
chosen is determined by the relevant $s$-field.

@ If we were to allocate the lookup-table for direct-mapping, the
table will waste $16 \times |MAX_SUBCUBOIDS|$ table elements. This is
because, in order to make every valid $s$ value indexable, we must
allocate |MAX_SFIELD| table elements per row. However, out of this row
only 26 elements actually contain a valid index to one of the
neighbouring subcuboids. Direct mapping is therefore space inefficient.


\bigskip

\centerline{\epsfig{file=figures/neighbour,scale=0.5}}

\bigskip

If we were to use a two-tiered mapping, however, we only waste space
equivalent to 16 table elements. This works by allocating an index
lookup-table, which maps the first |MAX_SFIELD| $s$ values to an index
that points to one of the 26 valid neighbouring subcuboids.

@d MAX_SUBCUBOIDS 1024
@d NUM_NEIGHBOURS 26 /* number of cuboid neighbours */
@d MAX_SFIELD 42 /* maximum $s$-field value: $\langle101010\rangle$ */
@<Global variables@>=
int8_t neighbour_idx_table[MAX_SFIELD + 1] = {
-1,  0,  1, -1,  2,  3,  4, -1,  5,  6,
 7, -1, -1, -1, -1, -1,  8,  9, 10, -1,
11, 12, 13, -1, 14, 15, 16, -1, -1, -1,
-1, -1, 17, 18, 19, -1, 20, 21, 22, -1,
23, 24, 25
};
uint32_t neighbour_table[MAX_SUBCUBOIDS][NUM_NEIGHBOURS];

@ Function |get_neighbour(s,i)| returns the next effective subcuboid
for a particle, where it was previously in subcuboid |i| and the most
recent application of a trajectory step updated the $s$-field to |s|.

@d cuboid_lookup(i,s) neighbour_table[(i)][neighbour_idx_table[s]]
@<Global functions@>=
uint32_t get_neighbour(uint32_t i, uint8_t s)
{
	if (s) return cuboid_lookup(i,s);
	return i; /* particle continues to exist in the same subcuboid */
}

@ Function |build_neighbour_table(l,m,n)| builds the subcuboid
neighbourhood lookup-table when the cuboid representing the simulation
world was divided into |l|, |m| and |n| equal parts along the $x$, $y$
and $z$ axes. There will be a total of $|l| \times |m| \times |n|$
subcuboids, which this function assumes is less than |MAX_SUBCUBOIDS|.

@d cuboid_assign(r,s,v) neighbour_table[(r)][neighbour_idx_table[(int)(s)]] = (v)
@<Global functions@>=
void build_neighbour_table(uint32_t l, uint32_t m, uint32_t n)
{
    uint32_t i, j, k, x, y, z, t = m * n;
    uint32_t r; /* row for subcuboid currently being filled in */
    uint8_t s, xb, yb, zb; /* $s$-field and extracted bit pairs */
    for (i = 0; i < l; ++i)
        for (j = 0; j < m; ++j)
	    for (k = 0; k < n; ++k) { /* set neighbours */
                r = i * t + j * n + k; /* row index for the current subcuboid */
		for (s = 1; s <= MAX_SFIELD; ++s) {
    		    x = i; y = j; z = k;
    		    if (((xb = s & 0x3) == 0x3) ||
        	        ((yb = s & 0xC) == 0xC) ||
			((zb = s & 0x30) == 0x30)) continue;
                    else @<Calculate the neighbour using the axes bit-pairs@>;
                }
            }
}

@ @<Calculate the neighbour using the axes bit-pairs@>=
{
@<Adjust index of the neighbour cuboid along the $x$-axis@>;
@<Adjust index of the neighbour cuboid along the $y$-axis@>;
@<Adjust index of the neighbour cuboid along the $z$-axis@>;
cuboid_assign(r,s, x * t + y * n + z);
}

@ @<Adjust index of the neighbour cuboid along the $x$-axis@>=
if (xb == 0x1) {
   if (--x < 0) @<Neighbour subcuboid is outside the simulation world@>;
} else if (xb == 0x2) {
   if (++x == l) @<Neighbour subcuboid is outside the simulation world@>;
}

@ @<Adjust index of the neighbour cuboid along the $y$-axis@>=
if (yb == 0x4) {
   if (--y < 0) @<Neighbour subcuboid is outside the simulation world@>;
} else if (yb == 0x8) {
   if (++y == m) @<Neighbour subcuboid is outside the simulation world@>;
}

@ @<Adjust index of the neighbour cuboid along the $z$-axis@>=
if (zb == 0x10) {
   if (--z < 0) @<Neighbour subcuboid is outside the simulation world@>;
} else if (zb == 0x20) {
   if (++z == n) @<Neighbour subcuboid is outside the simulation world@>;
}

@ All of the subcuboids that belong to the void outside the simulation
world are given the same subcuboid index |OUTSIDE_WORLD|.
@d OUTSIDE_WORLD (MAX_SUBCUBOIDS + 1)
@<Neighbour subcuboid is outside the simulation world@>=
{
    cuboid_assign(r,s, OUTSIDE_WORLD);
    continue;
}

@ Function |print_neighbour_table(f)| prints the subcuboid
neighbourhood table to the I/O stream pointed to by |f|.
@<Global functions@>=
void print_neighbour_table(FILE *f)
{
	uint32_t i, j, c, k = MAX_SUBCUBOIDS + 1;
	fprintf(f, "Subcuboid neighbourhood table:\n");
	for (i = 0; i < num_subcuboids; ++i) {
	    for (j = 0; j < NUM_NEIGHBOURS; ++j) {
	        c = neighbour_table[i][j];
	    	if (c < k) fprintf(f, "%3d ", c);
		else fprintf(f, " .  ");
	    }
	    fprintf(f, "\n");
	}	
}

@ The following code segment tests the functionalities provided by
this sections.

@<Test subcuboid functionalities@>=
{
	unsigned long l, m, n, i, j, k, e, r;
	double dx, dy, dz;
	Vector v;
	BoundingBox bb = {{0.0, 0.0, 0.0, 1.0},{0.0, 0.0, 0.0, 1.0}};
	do {
	       printf("Give upper bound (x, y, z) of the simulation world:\n");
	       scanf("%lf %lf %lf", &bb.u[0], &bb.u[1], &bb.u[2]);
	       printf("Give lower bound (x, y, z) of the simulation world:\n");
	       scanf("%lf %lf %lf", &bb.l[0], &bb.l[1], &bb.l[2]);
	} while (bb.u[0] < bb.l[0] ||
	         bb.u[1] < bb.l[1] ||
	         bb.u[2] < bb.l[2]);
        do {
	       printf("Give number of divisions on the x, y, and z axes:\n");
	       scanf("%lu %lu %lu", &l, &m, &n);
	} while (l * m * n > MAX_SUBCUBOIDS);	
	build_subcuboid_trees(&bb, l, m, n);
	print_subcuboid_search_trees(stdout, subcuboid_search_tree);
	@<Search subcuboid for all points in each of the subcuboids@>;

	build_neighbour_table(l, m, n);
	build_subcuboids_table(&bb, l, m, n);
	print_neighbour_table(stdout);
	print_subcuboids_table(stdout);
}

@ @<Search subcuboid for all points in each of the subcuboids@>=
dx = (bb.u[0] - bb.l[0]) / (double) l;
dy = (bb.u[1] - bb.l[1]) / (double) m;
dz = (bb.u[2] - bb.l[2]) / (double) n;
v[0] = bb.l[0];
for (i = 0; i < l; ++i, v[0] += dx) {
    v[1] = bb.l[1];
    for (j = 0; j < m; ++j, v[1] += dy) {
        v[2] = bb.l[2];
	for (k = 0; k < n; ++k, v[2] += dz) {
	    e = i * m * n + j * n + k;
	    r = find_subcuboid(subcuboid_search_tree, v);
	    if (e != r) goto test_subcuboid_search_failed;
        }
    }
}
printf("Test was a success...\n");
goto test_subcuboid_search_done;
test_subcuboid_search_failed:
       printf("Failure at (%lu, %lu, %lu): %lu instead of %lu for (%lf, %lf, %lf)\n",
       i, j, k, r, e, v[0], v[1], v[2]);
test_subcuboid_search_done:

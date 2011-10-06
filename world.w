@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@*1 Cuboid. The simulation world is defined by a {\sl cuboid}, which is
composed of solids. Only particles inside this cuboid are tracked
during a simulation, and particles outside this cuboid are marked as
``particles outside the simulation world''. To exploit parallelism
during a simulation, {\tt MCS} logically divides this cuboid into
smaller cuboids of equal shape and size. Each of these cuboids are
referred to as a {\sl subcuboid}. A set of subcuboids decomposes the
simulation world into disjointed volumes where particles can be
tracked independently of other particles in other subcuboids.

The simulation world cuboid must be axis-aligned with the world
coordinate frame. It is decomposed into subcuboids by dividing
along each of the three axes. Division can be carried out in any
combinations, however, all of the subcuboids resulting from a division
must all have the same size and shape (i.e., must be a cuboid), and
all of the union of all the subcuboids must exactly produce the
original cuboid that defines the simulation world. For instance, the
following are valid divisions of the simulation world:

\centerline{PUT FIGURE HERE}

We do not use space-partitioning data structures such as the {\sl
octree} or the {\sl $kd$-tree} here because it is important that the
subcuboids are of the same size and shape. With unequal partitions, we
must carry out a lot more calculations in order to determine to which
neighbour a particle has escaped. Even worse is that this computational
costs is non-uniform across the particles, as the following figure
demonstrates.

\centerline{PUT FIGURE HERE}

If we were simulating a few thousand particles, the overhead
computational cost may be considered acceptable. However, since {\tt
MCS} is expected to simulate close to millions, or even billions of
particles, the overhead becomes significant and therefore
unacceptable. By using subcuboids of the same size and shape, we can
avoid all of the computations (including the overhead) by using an
efficient lookup-table, for which the cost of finding a neighbour is
equivalent to two memory accesses for all particles throughout. We
shall discuss this in the following sections.

To accelerate particle tracking inside a subcuboid, however, we could
still take advantage of an octree or a $kd$-tree to accelerate the
search for the solid which contains a particle.

@*2 Relationship between subcuboids.
Once we are in a simulation, we need to determine at each step of the
trajectory which subcuboid currently contains the particle. This is
important because different subcuboids have different solid
configurations, and this determines the materials for the
physics processes.

If we assume that the planes that divide the simulation world cuboid
also divides the void outside this cuboid, then every subcuboid is
surrounded by other subcuboids (for internal subcuboids) or void
cuboids (for subcuboids that are on the boundary of the simulation
world). Either way, every subcuboid is surrounded by 26 logical
subcuboids. These are shown below:

\centerline{PUT FIGURE HERE}

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

Every particle has an $s$-field, which gives its location relative to
the current subcuboid containing it. Before a trajectory step is
applied for the first time, the $s$-field for a particle is
initialised to $\langle000000\rangle$ once. After applying a
trajectory step, which could have changed the particle's location, we
update the $s$-field using |update_sfield(s,bb,v)|. After subsequent
application of a trajectory step, we simply determine the next
effective subcuboid using the particle's current $s$-field, its
containing subcuboid, and a mapping table.

The mapping table is an $m \times 26$ two-dimensional array where each
row represents one of the $m$ subcuboids, and every element of the row
points to one of the 26 neighbouring subcuboids. Which column index is
chosen is determined by the relevant $s$-field.

@ If we were to allocate the lookup-table for direct-mapping, the
table will waste $16 \times |MAX_SUBCUBOIDS|$ table elements. This is
because, in order to make every valid $s$ value indexable, we must
allocate |MAX_SFIELD| table elements per row. However, out of this row
only 26 elements actually contain a valid index to one of the
neighbouring cuboids. Direct mapping is therefore space inefficient.

If we were to use a two-tiered mapping, however, we only waste space
equivalent to 16 table elements. This works by allocating an index
lookup-table, which maps the first |MAX_SFIELD| $s$ values to an index
that points to one of the 26 valid neighbouring cuboids.

@f int8_t int
@d MAX_SUBCUBOIDS 64
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
uint32_t num_subcuboids = 0;

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

@f uint8_t int
@d cuboid_assign(r,s,v) neighbour_table[(r)][neighbour_idx_table[(int)(s)]] = (v)
@<Global functions@>=
void build_neighbour_table(uint32_t l, uint32_t m, uint32_t n)
{
	int i, j, k; /* indices along $x$, $y$ and $z$ axes */
	int r; /* row for subcuboid currently being filled in */
	int x, y, z; /* indices along $x$, $y$ and $z$ axes for neighbour */
	int xb, yb, zb; /* extracted bit pairs from $s$-field */
	int t = m * n; /* cached value */
	uint8_t s; /* index for the entire row */
	num_subcuboids = t * l;
	for (i = 0; i < l; ++i)
	for (j = 0; j < m; ++j)
	for (k = 0; k < n; ++k)
	        @<Set neighbours for the current subcuboid@>;
}

@ @<Set neighbours for the current subcuboid@>=
{
r = i * t + j * n + k; /* row index for the current subcuboid */
for (s = 1; s <= MAX_SFIELD; ++s) {
    @<Initialise neighbour indices with those from the current subcuboid@>;
    if (((xb = s & 0x3) == 0x3) ||
        ((yb = s & 0xC) == 0xC) ||
	((zb = s & 0x30) == 0x30))
	continue;
    else @<Calculate the neighbour using the axes bit-pairs@>;
}
}

@ @<Initialise neighbour indices with those from the current subcuboid@>=
x = i; y = j; z = k;

@
@d OUTSIDE_WORLD (MAX_SUBCUBOIDS + 1)
@<Neighbour subcuboid is outside the simulation world@>=
{
    cuboid_assign(r,s, OUTSIDE_WORLD);
    continue;
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


@ @<Global functions@>=
void print_neighbour_table(FILE *f)
{
	int i, j, c, k = MAX_SUBCUBOIDS + 1;
	for (i = 0; i < MAX_SUBCUBOIDS; ++i) {
	    for (j = 0; j < NUM_NEIGHBOURS; ++j) {
	        c = neighbour_table[i][j];
	    	if (c < k) fprintf(f, "%3d ", c);
		else fprintf(f, " .  ");
	    }
	    fprintf(f, "\n");
	}	
}

@*1 Finding the subcuboid containing a particle.
The previous sections gives us an efficient way to determine where a
particle will go to if we know the containing subcuboid and the
particle's $s$-field. Now, how do we know where a particle is if we do
not know either the containing subcuboid or the particle's $s$-field?
This question arises when we are generating primary particles using a
particle gun. The only way to determine the containing subcuboid is to
search throughout the simulation world. Again, for a few thousand
particles, the cost for an exhaustive search could be
acceptable; however, {\tt MCS} is expected to simulate
millions. Hence, we must devise a better strategy.

@ First of all, we only search for particles that are inside the
simulation world. The, we take advantage of four design condition:

$$\vcenter{\halign{# & #\hfil \cr
1) & all of the subcuboids are of the same size and shape,\cr
2) & their union exactly defines the simulation world,\cr
3) & no two subcuboids intersect, and\cr
4) & the division into subcuboids is immutable for the entire simulation.\cr
}}$$

This means that, we can decompose a single search for the subcuboid
containing a given particle in the three-dimensional space into three
separate, but faster, searches in one-dimension. At the end of the
three searches, we simply combine the results to give the index of the
subcuboid.

@<Global variables@>=
double *cuboid_search_tree; /* stores three binary search trees */

@ If we us assume that the simulation world is divided into $l \times m
\times n$ subcuboids. After the three one-dimensional searches, let
$i$, $j$ and $k$ represent the search results along the $x$, $y$ and
$z$ respectively. Then the index of the subcuboid is given by
$(i \times m \times n + j \times n + k)$. The search on each axes is
carried out using a binary search tree.

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

@ Function |build_tree_subcuboid_search(t,n,u,l)| builds a perfect
binary search tree |t| where the upper and lower bounds of the
simulation world are respectively |u| and |l| units in the selected
axis. In the selected dimension, the simulation world has been divided
into |n| equal subcuboids. This function requires that $|n| = 2^k$ for
some positive integer $k$.
@<Global functions@>=
void build_tree_subcuboid_search(double *t, uint32_t n, double u, double l)
{
	int i, j, p;
	double d;
	@<Initialise the first two elements@>;
	@<Build the tree top-down by starting at the root@>;
	@<Now displace the tree to the correct position@>;
}

@ The first initialisation is required because subsequent calculations
use their parent node, which is calculated using a {\sl right
shift}. Although the root of the tree is |t[1]|, |t[0]| is the parent
of node |t[1]| according to the shift calculation, hence, it must be
initialised. The second initialisation ensures that the tree being
built is complete, no matter what the number of divisions |n| happens
to be.
@<Initialise the first two elements@>=
t[0] = 0.0;
t[1] = n / 2.0;

@ For each node, we first calculate the left child. If this is outside
the allowed node index, we are done. Otherwise, we calculate the
node's parent value, and use its displacement from the parent to
calculate the appropriate displacements for its children.
@<Build the tree top-down by starting at the root@>=
for (i = 1; i < n; ++i) {
    j = i << 1; /* left child */
    if (j < n) {
        p =  i % 2 ? (i - 1) >> 1 : i >> 1;
	@<Fill in the children with appropriate displacement@>;
    } else break;
}

@ @<Fill in the children with appropriate displacement@>=
d = (t[i] - t[p]) / 2.0;
if (d < 0.0) d *= -1;
t[j] = t[i] - d;
t[j + 1] = t[i] + d;

@ The tree array, as it currently is, assumes that the lower bound for
the simulation space is at the origin, and that the size of each
subcuboid is 1 units. We must therefore displace the value in the
tree nodes using the correct lower bound and subcuboid size in the
selected direction.
@<Now displace the tree to the correct position@>=
d = (u - l) / n;
for (i = 1; i < n; ++i) t[i] = (t[i] * d) + l;

@ Function |build_cuboid_trees(bb, l,m,n)| builds three complete
binary search trees. In each of these trees, every node is the lower
bound of the subcuboids produced after dividing the simulation world,
represented by the bounding box |bb|, into |l|, |m| and |n| equal
parts along the $x$, $y$ and $z$ axes. There will be a total of $|l|
\times |m| \times |n|$ subcuboids, which this function assumes is less
than |MAX_SUBCUBOIDS|.

@<Global functions@>=
bool build_subcuboid_trees(BoundingBox *bb, uint32_t l, uint32_t m, uint32_t n)
{
	double *x, *y, *z;
	@<Allocate memory for the three cuboid search trees@>;
	build_tree_subcuboid_search(x, l, bb->u[0], bb->l[0]);
	build_tree_subcuboid_search(y, m, bb->u[1], bb->l[1]);
	build_tree_subcuboid_search(z, n, bb->u[2], bb->l[2]);
	@<Finalise the subcuboid search trees@>;
	return true;
}

@ This allocates a single array to hold all of the three trees. Each
search tree with $n$ subcuboids only requires $n - 1$ tree
nodes. However, we allocate one extra element per tree because the 
root only begins at the second element (index 1). Instead of wasting
the first element, we use it to store the tree size.
@<Allocate memory for the three cuboid search trees@>=
cuboid_search_tree = mem_typed_alloc(l + m + n, double, mem_p);
if (cuboid_search_tree == NULL) return false;
x = cuboid_search_tree;
y = x + l;
z = y + m;

@ Don't forget to store the number of nodes as the first element of
the tree array.
@<Finalise the subcuboid search trees@>=
*(long *)x = l - 1;
*(long *)y = m - 1;
*(long *)z = n - 1;

@ @<Global functions@>=
void print_subcuboid_search_trees(FILE *f, double *t)
{
	long i, j, size;
	for (i = 0; i < 3; ++i) {
	    size = *(long *) t + 1;
	    fprintf(f, "tree[%ld]: ", size - 1);
	    for (j = 1; j < size; ++j) fprintf(f, "%lf ", t[j]);
	    fprintf(f, "\n");
   	    t += size; /* move to next tree array */
	}
}

@ Function |find_subcuboid(v)| finds the subcuboid that contains a
point with position vector |v| by searching the three binary search
trees pointed to by |t|. 
@<Global functions@>=
uint32_t find_subcuboid(double *t, Vector v)
{
	uint32_t i, j, f, c[3], size[3];
	for (i = 0; i < 3; ++i) {
	    @<Prepare the current subcuboid search tree@>;
	    @<Search for subcuboid in the current tree@>;
	    @<Move to the next subcuboid search tree@>;
	}
	return (c[0] * size[1] * size[2] + c[1] * size[2] + c[2]);
}
@ @<Prepare the current subcuboid search tree@>=
size[i] = *(long *) t + 1;
j = 1;

@ We start at the root and travel left or right, depending on the
value at the node, until we reach a leaf. Each of this left-right
decision is recorded in the bit-field |f|. Once we have reached the
leaf, the value of |f| gives the index of the subcuboid in the
selected dimension, which we copy to |c[i]| for later use in combining
the indices along the three axes.
@<Search for subcuboid in the current tree@>=
f = 0x0;
while (j < size[i]) {
    f <<= 1;
    if (v[i] > t[j]) {
        f |= 0x1;
	j = (j << 1) + 1; /* right subtree */
    } else j <<= 1; /* left subtree */
}
c[i] = f;

@ @<Move to the next subcuboid search tree@>=
t += size[i];

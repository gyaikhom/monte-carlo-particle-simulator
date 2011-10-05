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

The simulation world cuboid must be axis-aligned. They are decomposed
into subcuboids by dividing along each of the three axes. Division can be
carried out in any combinations, however, all of the subcuboids resulting
from a division must all have the same size and shape (i.e., must be a
cuboid), and all of the union of all the subcuboids must exactly produce
the original cuboid that defines the simulation world. For instance,
the following are valid divisions of the simulation world:

\centerline{PUT FIGURE HERE}

@*2 Relationship between subcuboids.
If we assume that the planes that divided the simulation world cuboid
also divided the void outside this cuboid, then every subcuboid is
surrounded by other subcuboids (for internal subcuboids) or void
cuboids (for subcuboids that are on the boundary of the simulation
world). In either case, every subcuboid is surrounded by 26 logical
subcuboids.

For a given subcuboid, we can define all of its 26 neighbours using
the six faces of the subcuboid. Let $f_i, 0 \le i < 6$ represent the
six faces so that the intervals $[f_0, f_1]$, $[f_2, f_3]$, and
$[f_4, f_5]$ respectively define the region inside the subcuboid
along the $x$, $y$ and $z$ axes. Now, for every point $(x, y, z)$ in
the three-dimensional space, we can define a bit-field
$s = \langle s_0, s_1, s_2, s_3, s_4, s_5 \rangle$, where $s_0$ is the
least significant bit. In the following sections, we shall refer to
the bit-field represented by $s$ as the {\sl $s$-field}.@^$s$-field@>

All of the bits in $s$ are zero, except in the following cases:

$s_0 = 1$ if $x < f_0$, $s_1 = 1$ if $x > f_1$,

$s_2 = 1$ if $y < f_2$, $s_3 = 1$ if $y > f_3$,

$s_4 = 1$ if $z < f_4$, and $s_4 = 1$ if $z > f_5$.

For instance, points inside the subcuboid will have $s =
\langle000000\rangle$;  whereas, points inside the top neighbour will
have $s = \langle000100\rangle$. All of the 26 neighbours are shown
below:

\centerline{PUT FIGURE HERE}

As we can see, the bit pairs $(s_0, s_1)$, $(s_2, s_3)$, and $(s_4,
s_5)$ will never be both nonzero. Furthermore, we can represent the
faces using a bounding box, where the points $(f_0, f_2, f_4)$ and
$(f_1, f_3, f_5)$ are respectively the lower and upper bounds.

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
applied for the first time, the $s$-field for a particle is set to
$\langle000000\rangle$. After applying a step, which could have
changed the particle's location, we update the $s$-field using the six
faces of the current subcuboid. During subsequent application of a
trajectory step, we simply use a mapping table to retrieve the
effective subcuboid for a given particle by using its current
$s$-field and the previous subcuboid.

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

@d cuboid_lookup(i,s) neighbour_table[(i)][neighbour_idx_table[(int)(s)]]
@<Global functions@>=
uint32_t get_neighbour(uint32_t i, uint8_t *s)
{
	if (s) return cuboid_lookup(i,s);
	return i; /* particle continues to exist in the same subcuboid */
}

@ Function |build_neighbour_table(l,m,n)| builds the subcuboid
neighbourhood lookup-table when the cuboid representing the simulation
world was divided into |l|, |m| and |n| equal parts along the $x$, $y$
and $z$ axes. There will be a total of $|l| \times |m| \times |n|$
subcuboids.

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
	if (num_subcuboids > MAX_SUBCUBOIDS) {
	    num_subcuboids = 0;
	    return;
	}
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

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

This table will waste $16 \times |MAX_SUBCUBOIDS|$ table elements,
since not all $s$ values are valid.

@d MAX_SUBCUBOIDS 64
@d OUTSIDE_WORLD (MAX_SUBCUBOIDS + 1)
@d MAX_SFIELD 42 /* maximum $s$-field value: $\langle101010\rangle$ */
@<Global variables@>=
uint32_t neighbour_table[MAX_SUBCUBOIDS][MAX_SFIELD];
uint32_t num_subcuboids = 0;

@ Function |get_neighbour(s,i)| returns the next effective subcuboid
for a particle, where it was previously in subcuboid |i| and the most
recent application of a trajectory step updated the $s$-field to |s|.
@<Global functions@>=
uint32_t get_neighbour(uint32_t i, uint8_t *s)
{
	if (s) return neighbour_table[i][(int)s - 1];
	return i; /* particle continues to exist in the same subcuboid */
}

@ Function |build_neighbour_table(l,m,n)| builds the subcuboid
neighbourhood lookup-table when the cuboid representing the simulation
world was divided into |l|, |m| and |n| equal parts along the $x$, $y$
and $z$ axes. There will be a total of $|l| \times |m| \times |n|$
subcuboids.
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
	for (k = 0; k < n; ++k) {
	        @<Set neighbours for the current subcuboid@>;
	}
}

@ @<Set neighbours for the current subcuboid@>=
r = i * t + j * n + k; /* row index for the current subcuboid */
for (s = 1; s <= MAX_SFIELD; ++s) {
    x = i; y = j; z = k;
    if (((xb = s & 0x3) == 0x3) ||
        ((yb = s & 0xC) == 0xC) ||
	((zb = s & 0x30) == 0x30)) {
	neighbour_table[r][s - 1] = OUTSIDE_WORLD;
        continue;
    } else {
    printf("%d %d\n", s, s % 26);

        @<Calculate the neighbour using the axes bits@>;
    }
}

@ @<Calculate the neighbour using the axes bits@>=
@<Adjust index of the neighbour cuboid along the $x$-axis@>;
@<Adjust index of the neighbour cuboid along the $y$-axis@>;
@<Adjust index of the neighbour cuboid along the $z$-axis@>;
neighbour_table[r][s - 1] = x * t + y * n + z;

@ @<Adjust index of the neighbour cuboid along the $x$-axis@>=
if (xb == 0x1) {
   if (--x < 0) @<Neighbour subcuboid is outside the simulation world@>;
} else if (xb == 0x2) {
   if (++x == l) @<Neighbour subcuboid is outside the simulation world@>;
}

@ @<Neighbour subcuboid is outside the simulation world@>=
{
    neighbour_table[r][s - 1] = OUTSIDE_WORLD;
    continue;
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
	int i, s, k = MAX_SUBCUBOIDS + 1;
	for (i = 0; i < MAX_SUBCUBOIDS; ++i) {
	    for (s = 0; s < MAX_SFIELD; ++s)
	    	if (neighbour_table[i][s] < k)
		fprintf(f, "%3d ", neighbour_table[i][s]);
		else fprintf(f, " .  ");
	    fprintf(f, "\n");
	}	
}

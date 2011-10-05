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
applied for the first time, the $s$-field for a prticle is set to
$\langle000000\rangle$. After applying a step, which could have
changed the particles location, we update the $s$-field using the six
faces of the containing subcuboid. During subsequent application of a
trajectory step, we simply use a mapping table to retrieve the
effective subcuboid for a given particle by using the index of the
previous subcuboid and its current $s$-field.


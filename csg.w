@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@** Constructive Solid Geometry.
During a simulation event, the tracking of a particle's trajectory
happens inside a closed three-dimensional simulation world composed of
solids. Each of these solids is defined by a {\sl binary tree}
data-structure, which is built recursively bottom-up from primitive
convex solids through affine transformations and boolean combination
of existing solids. A simulation world may be composed of several
non-intersecting solids, and hence, it is best represented by a {\sl
forest} of binary trees. Once the forest is ready, we pre-process it
to generate a compact data-structure that is optimised for efficient
processing during particle tracking.

We use {\sl Constructive Solid Geometry} to represent solids. Most of
the concepts used in this implementation are based on the books {\sl
Geometric and Solid Modelling: An Introduction} by Christoph
M. Hoffman [Morgan Kaufmann, ({\bf 1989})] and {\sl An Integrated
Introduction to Computer Graphics and Geometric Modelling} by Ronald
Goldman [CRC Press ({\bf 2009})].

@*1 Primitive solids.
{\tt MCS} currently supports four types of {\sl primitive solids}: the
{\sl parallelepiped}@^parallelepiped@>, the {\sl sphere}@^sphere@>,
the {\sl cylinder}@^cylinder@>, and the {\sl torus}@^torus@>.

@<Type definitions@>=
typedef enum {
	BLOCK = 0, SPHERE, CYLINDER, TORUS
} Primitive_type;

@ The standard primitives listed above are generic; solid instances of
these primitives can take different shapes and form. To use solids of a
given primitive type, a new instance of the generic primitive
must first be created. This must then be initialised to the required
specification by filling in the relevant parameters.

An instance of a primitive solid defines a convex solid volumes. A solid
is convex if any line segment between any two points on the surface of
the solid only intersects the solid at those two points. We represent
an instance of a primitive solid using the |Primitive| data type. 

{\bf NOTE:}
In each of the following sections, we specify only the relevant
details in relation to the context of the section. Additional details
will be incorporated when we discuss other aspects of the solids. For
instance, at the moment we are only concerned with the geometry of the
solids, and not with their material properties. Hence, the
details concerning their material properties will be added later on,
when we discuss materials.

@<Type definitions@>=
typedef struct {
	@<Information common to all primitives@>;
	union {
       	        @<Data for primitive block@>;
       		@<Data for primitive sphere@>;
       		@<Data for primitive cylinder@>;
       		@<Data for primitive torus@>;
	};
} Primitive;

@ The first field of all primitive data stores a primitive type. This
is used while deciding the manner in which a solid instance must be
processed: based on this type, we choose the appropriate \CEE/ union
field.

@<Information common to all primitives@>=
Primitive_type type;

@ In addition to the {\sl global coordinate frame} defined by the
simulation world, each initialised primitive also defines a {\sl local
coordinate frame}. The origin of this coordinate frame is used by the
geometry construction algorithm to conceptually place a primitive
inside the simulation world. By design, a primitive solid is created
in such a way that its origin always coincides with the origin of the
world coordinate frame. Any translation or transformation henceforth
is recorded separately in a binary tree. This will become clearer once
we reach the section on {\sl Constructive Solid Geometry Tree}.
@^global coordinate frame@>@^local coordinate frame@>
@^Constructive Solid Geometry Tree@>

@*2 Parallelepiped.
For simplicity, we shall refer to a parallelepiped as
a {\sl block}@^block@>. A primitive block stores the following
information.

@<Data for primitive block@>=
struct {
       @<Information that defines a primitive block@>;
} b;

@ The geometry of a block is defined by its length, width and
height. The origin of the block's local coordinate frame is defined by
its {\sl centroid}@^centroid@>, and the orthogonal axes incident on
this origin, i.e., the $x$, $y$ and $z$ axes of the local coordinate
frame, are aligned respectively in parallel to its length, height and
width. For efficient containment testing, we store half-lengths,
instead of storing the full length, width, or height.

@<Information that defines a primitive block@>=
double length, height, width; /* half-lengths */

@*2 Sphere.
A primitive sphere stores the following information.

@<Data for primitive sphere@>=
struct {
       @<Information that defines a primitive sphere@>;
} s;

@ The geometry of a sphere is defined by its radius, and the origin of
its local coordinate frame is defined by the sphere's {\it
center}@^center@>.

@<Information that defines a primitive sphere@>=
double radius;

@*2 Cylinder.
A primitive cylinder stores the following information.

@<Data for primitive cylinder@>=
struct {
       @<Information that defines a primitive cylinder@>;
} c;

@ The geometry of a cylinder is defined by its base radius and its
height.  The origin of the block's local coordinate frame is defined by
its {\sl centroid}@^centroid@>, and the orthogonal axes incident on
this origin, i.e., the $x$, $y$ and $z$ axes of the local coordinate
frame, are aligned so that $y$ axis is perpendicular to the base.

@<Information that defines a primitive cylinder@>=
double radius, height;

@*2 Torus.
A primitive torus stores the following information.

@<Data for primitive torus@>=
struct {
       @<Information that defines a primitive torus@>;
} t;

@ The origin of the torus' local coordinate frame is defined by
its {\sl center}@^center@>, which coincides with the origin of the
world coordinate frame. We also assume that they are radially
symmetrical to the $y$-axis of the world coordinate frame. Under these
condition, a torus may be defined parametrically as follows:

$$\vcenter{\halign{\hfil $#$ & $#$ & $#$ \hfil \cr
    x(\phi, \theta) & = & (R + r \cos{\theta}) \cos{\phi} \cr
    y(\phi, \theta) & = & (R + r \cos{\theta}) \sin{\phi} \cr
    z(\phi, \theta) & = & r \sin{\theta} \cr}}$$

\noindent where, the length $R$ gives the distance from the center of
the tube to the center of the torus, and the length $r$ gives the
radius of the tube. They are referred to as the {\sl
major radius}@^major radius@> and the {\sl minor radius}@^minor radius@>,
respecticely; and their ratio is referred to as the {\sl aspect
ratio}@^aspect ratio@>. Parameters $\phi$ and $\theta$ are
angles in radians, where $0 \le \phi, \theta < 2\pi$.

Parameter $\phi$ is the angle subtended by the center of the tube
on the $xz$-plane due to a counter-clockwise rotation about the
$y$-axis, where the rotation begins at the positive
$x$-axis. Parameter $\theta$, on the other hand, is the angle
subtended by the surface of the tube on a cross-section
of the tube, where the cross section is defined by a plane
perpendicular to the $xz$-plane which passes through the center of the
torus and the center of the tube, at the point where the angle
subtended by this plane on the $xz$-plane is $\phi$ radians. 

@<Information that defines a primitive torus@>=
double major, minor;
double phi, phi_start, theta, theta_start; /* subtended angle, and
start angle (in degrees) */

@ To create a primitive solid, we use the global memory area
|mem|. This memory area stores all of the primitive solids in blocks
of items of type |Primitive|. Hence, we need a pointer that points to
the slot that is currently available inside the active |Primitive|
block, and a pointer to check if the available slot is indeed
valid. The pointers |next_primitive| and |bad_primitive| keep track of
this information.

The |next_primitive| pointer is updated after the allocation of a slot
from the current |Primitive| block, or when a new |Primitive| block is
allocated using |mem_alloc(n,s)|; the |bad_primitive| is updated only
when a new block is allocated.
@<Global variables@>=
static Primitive *next_primitive = NULL;
static Primitive *bad_primitive = NULL;

@ To store a primitive solid, we check if the available slot is
valid, i.e., |next_primitive != bad_primitive|. If valid, we use this
slot; otherwise, we allocate a new block that can accommodate
|primitives_per_block| primitive solids, and use the first slot in
this block.

@d primitives_per_block 195 /* 16384 bytes per block */
@<Global functions@>=
Primitive *create_primitive_solid() {
        Primitive *slot = next_primitive;
	if (slot == bad_primitive) {
	   slot = mem_typed_alloc(primitives_per_block, Primitive, mem_phase_one);
	   if (slot == NULL) return NULL;
	   else {
	   	next_primitive = slot + 1;
		bad_primitive = slot + primitives_per_block;
	   }
	} else next_primitive++;
        return slot;
}

@*1 Constructive Solid Geometry Tree.
All solids in the simulation world are built from instances of
primitive solids by using volumetric boolean operators. These
operators work on the volume defined by the primitive solids. There
are three volumetric boolean operators: {\sl union}@^union@>,
{\sl intersection}@^intersection@> and {\sl
difference}@^difference@>. They respectively define the volumetric
union, intersection and difference of two solid operands.

Further to the boolean operators, affine transformation operators are
defined for translation, rotation and scaling. They allow us to change
the shape, form, orientation and location of its solid operand
relative to the world coordinate frame. Affine transformation
operators are unary, whereas boolean operators are binary.

@<Type definitions@>=
typedef enum {
	PRIMITIVE = 0, UNION, INTERSECTION, DIFFERENCE,
	TRANSLATE, ROTATE, SCALE, PARAMETER
} csg_node_type;

@ A CSG solid is represented by a binary tree where the root
represents the solid that is being created. The primitive solids are
stored at the leaves, and all of the internal nodes either represent
affine transformations or boolean combinations. This tree is generally
referred to as the {\sl Constructive Solid Geometry Tree}, or {\it CSG
Tree} for short. @^Constructive Solid Geometry Tree@>@^CSG Tree@>

The performance of many of the algorithms in the following sections
can be improved significantly if each of the CSG tree nodes are
enclosed inside an axis-aligned parallelipiped, commonly referred to
as the solid's {\sl bounding box}. This is defined by a pair of three
dimensional coordinates $(\lambda, \tau)$, where $\lambda$ and $\tau$
are respectively the lower and upper bounds of all the points that are
inside the solid.

@<Type definitions@>=
typedef struct {
	Vector l, u; /* lower and upper bounds */
} BoundingBox;

@ @<Type definitions@>=
typedef struct csg_node_struct CSG_Node;
struct csg_node_struct {
       @<Information common to all CSG nodes@>;
       union {
               @<Information stored in CSG leaf node@>;
	       @<Information stored in CSG internal node@>;
       };
};

@ Each node stores a node name, that uniquely identifies the
node in the CSG tree. They are supplied by the user through the
geometry input file. Node names are used during registration of
solids.

@d MAX_LEN_NODE_NAME 64
@<Information common to all CSG nodes@>=
char name[MAX_LEN_NODE_NAME];

@ To differentiate between node types, each node stores a
CSG node type. This field determines if a node is a primitive
solid node, a parameter node, or an intermediate solid node. Since we
only require eight values (i.e., |csg_node_type|), it will be a waste
to use more than three bits. Instead, we will use a bit field:

\medskip

\centerline{\epsfig{file=figures/node-type-bit-field,scale=1}}

\medskip

We make the three least significant bits to represent the node type. We
are choosing the least significant bits because it will be used more
often than the rest of the bit fields, which only come into effect
during error diagnosis. We then use the fourth least significant bit
to mark if the node is currently used in a solid CSG tree. Because the
CSG tree is built bottom-up, these is a possibility that the user has
made a mistake and added stray CSG nodes. Finally, the rest of the 28
bits are used to record the line number in the geometry file which
gave the command to create the node. This will be displayed to the
user during error diagnostics (e.g., to display a warning to the user
where in the geomtry file an unused CSG node was created).

@<Information common to all CSG nodes@>=
uint32_t op;

@ The following macros set the bit fields.
@d BIT_MASK_NODE 0x7
@d BIT_MASK_UBIT 0x8
@d BIT_MASK_LINE 0xfffffff0
@d set_ubit(n) ((n)->op |= BIT_MASK_UBIT)
@d set_line(n) ((n)->op |= (input_file_current_line << 4))
@d get_line(n) ((n)->op >> 4)

@ The following macros check the bit fields.
@d is_inuse(n) (BIT_MASK_UBIT & (n)->op)
@d is_primitive(n) (!(BIT_MASK_NODE & (n)->op)) /* special because |PRIMITIVE == 0x0| */
@d is_parameter(n) (PARAMETER == (BIT_MASK_NODE & (n)->op))
@d is_union(n) (UNION == (BIT_MASK_NODE & (n)->op))
@d is_intersection(n) (INTERSECTION == (BIT_MASK_NODE & (n)->op))
@d is_difference(n) (DIFFERENCE == (BIT_MASK_NODE & (n)->op))
@d is_translate(n) (TRANSLATE == (BIT_MASK_NODE & (n)->op))
@d is_rotate(n) (ROTATE == (BIT_MASK_NODE & (n)->op))
@d is_scale(n) (SCALE == (BIT_MASK_NODE & (n)->op))

@ All nodes store a pointer to its parent node. This is used for
backtracking during tree traversal.
@<Information common to all CSG nodes@>=
CSG_Node *parent; /* pointer to parent node */

@ The internal nodes of a CSG tree store operators. If the operator
stored is binary (i.e., the set-theoretic operators) both left and
right children point to solids. If the operator is unary, the left
child points to a solid, and the right child stores the operator
parameters (e.g., the displacement if we are translating a solid, or
the scaling factor if we are scaling a solid, etc.).

@<Information stored in CSG internal node@>=
struct {
       CSG_Node *left, *right; /* pointers to the left and right subtrees */
} internal;

@ CSG leaf nodes that point to solids, or unary operator nodes that
store parameters, require additional data. A node that corresponds to a
unary operator stores operator specific parameters, whereas a node that
corresponds to a solid stores a pointer to a primitive solid. A leaf
node, however, never points to an {\sl intermediate
solid}@^intermediate solid@>, which is defined as a complex solid
composed of several primitive solids using a CSG tree of operators. In
fact, each internal node in a CSG tree defines an intermediate solid.

@<Information stored in CSG leaf node@>=
union {
       @<Translation parameters@>;
       @<Rotation parameters@>;
       @<Scaling parameters@>;
       Primitive *p;
} leaf;

@ To create a new CSG node, we use the global memory area |mem|. This
memory area stores all of the CSG nodes in blocks of items of type
|CSG_Node|. Hence, we need a pointer that points to the slot that is
currently available inside the active |CSG_Node| block, and a pointer
to check if the available slot is indeed valid. The pointers
|next_csg_node| and |bad_csg_node| keep track of this information.

The |next_csg_node| pointer is updated after the allocation of a slot
from the current |CSG_Node| block, or when a new |CSG_Node| block is
allocated using |mem_alloc(n,s)|; the |bad_csg_node| is updated only
when a new block is allocated.
@<Global variables@>=
static CSG_Node *next_csg_node = NULL;
static CSG_Node *bad_csg_node = NULL;

@ To create a new CSG node, we check if the available slot is
valid, i.e., |next_csg_node != bad_csg_node|. If valid, we use this
slot; otherwise, we allocate a new block that can accommodate
|csg_node_per_block| CSG nodes, and use the first slot in this block.

@d csg_nodes_per_block 128 /* 65536 bytes per block */
@<Global functions@>=
CSG_Node *create_csg_node() {
        CSG_Node *slot = next_csg_node;
	if (slot == bad_csg_node) {
	   slot = mem_typed_alloc(csg_nodes_per_block, CSG_Node, mem_phase_one);
	   if (slot == NULL) return NULL;
	   else {
	   	next_csg_node = slot + 1;
		bad_csg_node = slot + csg_nodes_per_block;
	   }
	} else next_csg_node++;
	@<Initialise affine matrices to the identity matrix@>;
        return slot;
}

@*2 Nodes repository.
In the first phase, when the input files are processed, we must build
the CSG trees which correspond to all of the solids bottom-up. This
requires maintaining the CSG nodes already defined, which may be
disconnected from the others nodes, so that subsequent commands may
combine them using CSG operators. Hence, we use a nodes repository to
stores all of the nodes already defined. This is implemented using a
hash table, where the name of the node is used as the hash key.

@d MAX_CSG_NODES 65536 /* maximum number of CSG nodes allowed */
@<Type definitions@>=
typedef struct {
       uint16_t stat[7]; /* geometry statistics */
       CSG_Node *root; /* pointer to the CSG root */
       CSG_Node *table[MAX_CSG_NODES]; /* hash table of nodes */
} NodesRepository;
NodesRepository *nodes_repo = NULL;

@ The nodes repository is only used in phase one, when the data from
the input files are processed at the beginning. In phase two, this
information will be reduced to a compact representation which requires
less memory. After phase two is complete, we must free the resources
allocated to the nodes repository as it consumes quite a lot of
memory. Of course, since the nodes repository is allocated in the
memory area |mem_phase_one|, it will be freed automatically when we
free phase one memory area at the end of phase two.

Note that the macro |mem_typed_alloc()| initialises all of the fields
to zero, since it uses the |calloc()| system call to allocate the
memory. Hence, we do not need to initiliase the fields separately.

@<Global functions@>=
bool create_nodes_repository()
{
    nodes_repo = mem_typed_alloc(1, NodesRepository, mem_phase_one);
    if (NULL == nodes_repo) return false;
    return true;
}

@ Function |print_geom_statistics(f)| prints the geometry statistics
in the nodes repository to the I/O stream pointed to by |f|.
Inside the nodes repository, we only maintain the overall statistics
concerning the number of primitives, operators and solids. If we wish
to obtain statistics for a specific solid, it may be derived at
by traversing the CSG tree that corresponds to the selected solid. The
array stores in order the number of: primitives, unions,
intersections, differences, translations, rotations and scalings.

@<Global functions@>=
void print_geom_statistics(FILE *f)
{
    if (NULL == nodes_repo) return;
    fprintf(f, "Geometry statistics\n\tprimitive: %u\n\tunion: %u\n\tintersection: %u\n\tdifference: %u\n\ttranslation: %u\n\trotation: %u\n\tscaling: %u\n",
    nodes_repo->stat[0], nodes_repo->stat[1],
    nodes_repo->stat[2], nodes_repo->stat[3],
    nodes_repo->stat[4], nodes_repo->stat[5],
    nodes_repo->stat[6]);
}

@ The hash code for a string $c_1 c_2 \ldots c_l$ of length $l$ is a
 nonlinear function of the characters. We borrow the hash function
 described by Donald E. Knuth in the book {\sl The Stanford
 GraphBase}~[Addison-Wesley ({\bf 1993})] to calculate $h$. As noted
 therein, this hash function only works with ASCII encoded strings.

@d HASH_MULT 314159 /* random multiplier */
@d HASH_PRIME 516595003 /* the 27182818th prime; which is less than $2^{29}$ */

@<Global functions@>=
long hash(const char *name, long ubound)
{
        register long h;
        const char *t = name;
        for (h = 0; *t; t++) {
                h += (h ^ (h >> 1)) + HASH_MULT * (unsigned char) *t;
	        while (h >= HASH_PRIME) h -= HASH_PRIME;
        }
        return (h % ubound);
}

@ Function |register_csg_node(n,s)| registers into the nodes
repository the CSG node pointed to by |n| using the unique name
|s|. If any of the inputs are invalid, or if the supplied name is
already in use, it returns |false| to notify registration failure; 
otherwise, it returns |true| to indicate success. 

@<Global functions@>=
bool register_csg_node(CSG_Node *n, char *s) {
	 long h;
	 if (NULL == n || NULL == s) return false;
	 h = hash(s, MAX_CSG_NODES);
	 if (NULL == nodes_repo->table[h]) {
                 strcpy(n->name, s);
	         nodes_repo->table[h] = n;
		 nodes_repo->root = n; /* set as current root */
		 return true;
	 }
	 return false;
}

@ Function |find_csg_node(s)| finds in the nodes repository the CSG
node named |s|. If the input is invalid, or the node does not exists,
NULL is returned; otherwise, a pointer to the node is returned.
@<Global functions@>=
CSG_Node *find_csg_node(char *s) {
	 if (NULL == s) return NULL;
	 return nodes_repo->table[hash(s, MAX_CSG_NODES)];
}

@*2 The geometry input file.
The geometry of the solids and their placement and orientation within
the world is specified in the input file. The grammar for this input
file is very simple. The file consist of several commands, where an
entire line of text is used to specify a specific command. Each
command has the following format:

\smallskip

$\langle command \rangle \langle parameters \rangle \langle newline \rangle$

\smallskip

Here, $\langle command \rangle$ is a single character code, which
defines its intended action. The commands and their intended actions
are as follows: 

\smallskip

$$\vcenter{\halign{\hfil {\tt #} & # \hfil \cr
B, S, C, T & Create a primitive block, sphere, cylinder, or torus.\cr
u, i, d & Carry out a union, intersection, or difference.\cr
t, r, s & Translate, rotate, or scale the solid.\cr
+ & Registers the current CSG tree as solid.\cr
* & Defines the simulation world as an axis-aligned bounding box.\cr
\%\ & Begin comment line. Stop at the first newline character.\cr
}}$$

\smallskip

To support this intended action, the user must supply all of the
required parameters in the $\langle parameters \rangle$
field. Finally, every command must be terminated by a $\langle newline
\rangle$ character. The following is an example:

\bigskip

{\tt
\% Define primitive solids

T ("Torus A" 0.0 359.999999 10.0 2.0)

C ("Cylinder A" 10.0 20.0)

C ("Cylinder B" 10.0 20.0)

\

\% Operation on primitive solids

u ("U1" "Torus A" "Cylinder A")

d ("D1" "U1" "Cylinder B")

t ("T1" "D1" 10.0 50.0 20.0)

s ("S1" "T1" 10.0 20.0 30.0)

+ ("SOLID") \% register CSG tree as a solid

* (-100.0 -100.0 -100.0 100.0 100.0 100.0) \% define the simulation world
}

\bigskip

@ Function |read_geometry(n)| reads the geometry data by opening the
file named |n|.
@<Global functions@>=
bool read_geometry(const char *n)
{
        @<Local variables: |read_geometry()|@>;
	@<Open input geometry file and initialise nodes repository@>;
	while (EOF != (c = fgetc(f))) {
		@<Discard comments, white spaces and empty lines@>;
	        @<Process input command@>;
	}
	fclose(f);
	return true;
	
	@<Handle geometry file errors@>;
	return false;
}

@ @<Local variables: |read_geometry()|@>=
FILE *f;
char c;

@ We use the temporary variable |input_file_name| to store the name of
the file that is currently being processed, so that if we wish, we
might generate the name of the file at runtime (e.g., by appending
prefix and suffix strings to a base filename).

@<Open input geometry file and initialise nodes repository@>=
input_file_current_line = 1;
strcpy(input_file_name, n); /* generate filename */
f = fopen(input_file_name, "r");
if (NULL == f) {
        fprintf(stderr, "Failed to open input geometry file: %s\n",
	        input_file_name);
	return 1;
}
if (false == create_nodes_repository()) {
        fprintf(stderr, "Failed to initialise nodes repository\n");
	return 1;
}

@ To improve readability, we allow comments, empty lines and
indentation of commands.

@<Discard comments, white spaces and empty lines@>=
@<Discard comment lines@>;
@<Discard indentations@>;
@<Discard empty lines@>;

@ Comments begin with the `{\tt \%}' character, and end after the
end-of-line character.
@<Discard comment lines@>=
if (c == '%') {
        while (EOF != (c = fgetc(f)))
                if ('\n' == c) break; /* gobble comments */
	if (EOF == c) break; /* done reading input file */
}

@ @<Discard indentations@>=
if (' ' == c || '\t' == c) continue;

@ @<Discard empty lines@>=
if ('\n' == c) {
        ++input_file_current_line;
	continue;
}

@ @<Process input command@>=
switch(c) {
case 'B': @<Read block geometry@>;@+ break;
case 'S': @<Read sphere geometry@>;@+ break;
case 'C': @<Read cylinder geometry@>;@+ break;
case 'T': @<Read torus geometry@>;@+ break;
case 'u': @<Read union operation@>;@+ break;
case 'i': @<Read intersection operation@>;@+ break;
case 'd': @<Read difference operation@>;@+ break;
case 't': @<Read translation operation@>;@+ break;
case 'r': @<Read rotation operation@>;@+ break;
case 's': @<Read scaling operation@>;@+ break;
case '+': @<Read solid registration@>;@+ break;
case '*': @<Read simulation world specification@>;@+ break;
default:@/
	fprintf(stderr, "%s[%u] Invalid command '%c' in input file\n", input_file_name, input_file_current_line, c);
        goto error_invalid_file;
}

@ When we cannot recover from an error (e.g., incorrect input file),
we must exit the system after cleaning up the resources that were
allocated by previous commands. Furthermore, the system must also
alert the user about the error. This section defines all of the exit
points and the corresponding error messages.

@<Handle geometry file errors@>=
@<Alert error while reading file@>;
@<Alert failure to create primitive solid@>;
@<Alert failure to create operator node@>;
@<Alert failure to create parameter node@>;
@<Alert solid already exists@>;
@<Alert solid does not exists@>;
@<Alert invalid division of simulation world@>;
error_invalid_file:@/
@<Cleanup resources allocated to invalid geometry@>;
fclose(f);

@ @<Alert error while reading file@>=
failed_read_exit_after_cleanup:@/
fprintf(stderr, "Failed to read file '%s' at line %u\n"
        "\tInvalid formatting of parameters\n",
	input_file_name, input_file_current_line);
goto error_invalid_file;

@ @<Alert failure to create primitive solid@>=
create_primitive_failed_exit_after_cleanup:@/
fprintf(stderr, "%s[%u] Failed to create primitive solid\n",
	input_file_name, input_file_current_line);
goto error_invalid_file;

@ @<Alert failure to create operator node@>=
create_operator_failed_exit_after_cleanup:@/
fprintf(stderr, "%s[%u] failed to create internal node\n",
	input_file_name, input_file_current_line);
goto error_invalid_file;

@ @<Alert failure to create parameter node@>=
create_parameter_failed_exit_after_cleanup:@/
fprintf(stderr, "%s[%u] failed to create leaf node\n",
input_file_name, input_file_current_line);
goto error_invalid_file;

@ @<Alert solid already exists@>=
solid_exists_exit_after_cleanup:@/
fprintf(stderr, "%s[%u] Invalid geometry specification... "@/
"Solid named '%s' already exists\n", input_file_name,
input_file_current_line, solid_name);
goto error_invalid_file;

@ @<Alert solid does not exists@>=
no_solid_exists_exit_after_cleanup:@/
fprintf(stderr, "%s[%u] Invalid geometry specification... "@/
	"Solid named '%s' does not exists\n", input_file_name,
	input_file_current_line, solid_name);
goto error_invalid_file;

@ @<Alert invalid division of simulation world@>=
invalid_subcuboid_division_exit_after_cleanup:@/
fprintf(stderr, "%s[%u] Invalid geometry specification... "@/
	"Number of subcuboids %u exceeds allowed %u\n",
	input_file_name, input_file_current_line,
	num_subcuboids, MAX_SUBCUBOIDS);
goto error_invalid_file;

@ @<Cleanup resources allocated to invalid geometry@>=
fprintf(stderr, "Cleaning up resources...\n");

@ While reading in the operations, we use the following variables to
store temporary values. Since each operation is defined independently
of other operations before and after, these variables can be shared by
all of the operation commands without confusion.

@d MAX_LEN_FILENAME 256
@<Global variables@>=
double op_x, op_y, op_z, op_theta;
char op_target[MAX_LEN_NODE_NAME], op_solid[MAX_LEN_NODE_NAME];
char op_left[MAX_LEN_NODE_NAME], op_right[MAX_LEN_NODE_NAME];
char *solid_name; /* points to name of solid while alerting error */
CSG_Node *target_solid, *left_solid, *right_solid;

@ The parameters that are required to initialise a standard primitive
is read from an input file. These parameters must be supplied using
specific formatting requirements. The parameters for the block
geometry are supplied in the following format.

\smallskip

(``name" $length$ $width$ $height$)

\smallskip

For instance, the specification {\tt ("Block A" 10.0 5.0 15.0)}
initialises a new block named ``Block A" with centroid at (0.0, 0.0,
0.0) in the world coordinate frame with length 10.0, width 5.0 and
height 15.0. It is important that all of these lengths are specified
using the same unit of measurement. We shall discuss in later sections
how a unit of measurement is chosen.

@<Read block geometry@>=
@<Create a new primitive solid@>;
@<Initialise primitive block with relevant data@>;
@<Register the primitive solid@>;

@ @<Create a new primitive solid@>=
p = create_primitive_solid();
if (NULL == p) {
        @<Exit after cleanup: failed to create primitive solid@>;
}

@ To increase clarity and to reduce code size, we keep recurring error
messages in common code-blocks, and jump to these blocks as
required. When it comes to implementing error handling, such as this,
where we wish to terminate the application after cleanup, I prefer
{\tt goto} statements. This is where bending the rules of {\it
Structured Programming}@^Structured Programming@> leads to cleaner
source code and smaller object code.

@<Exit after cleanup: failed to create primitive solid@>=
goto create_primitive_failed_exit_after_cleanup;

@ We use the following variables while reading input data from a file.

@<Global variables@>=
char input_file_name[MAX_LEN_FILENAME]; /* current input file */
uint32_t input_file_current_line; /* current line in current input file */
int read_count; /* returned by fscanf() */
Primitive *p; /* used during initialisation of primitive solids */
CSG_Node *internal_node, *leaf_node, *temp_node;

@ @<Initialise primitive block with relevant data@>=
read_count = fscanf(f, "(\"%[^\"]\" %lf %lf %lf)",
       op_solid, &p->b.length, &p->b.width, &p->b.height);
if (EOF == read_count || 4 != read_count) {
        @<Exit after cleanup: failed to read from file@>;
}
p->type = BLOCK;
@<Prepare block for containment testing@>;
++(nodes_repo->stat[PRIMITIVE]);

@ To improve containment testing, we precalculate some of the values
that are used frequently. We first half the dimensions, and then
calculate the containment range for the block.

@<Prepare block for containment testing@>=
@<Half the length, height and width of the block@>;
@<Calculate containment range for the block@>;

@ When checking containment, we require the halves of the block
length, height and width. Hence, although the dimensions were
specified in full, we shall store their halves internally.

@<Half the length, height and width of the block@>=
p->b.length /= 2.0;
p->b.height /= 2.0;
p->b.width /= 2.0;

@ @<Exit after cleanup: failed to read from file@>=
goto failed_read_exit_after_cleanup;

@ @<Register the primitive solid@>=
leaf_node = create_csg_node();
if (NULL == leaf_node) {
        @<Exit after cleanup: failed to create leaf node@>;
} else {
        leaf_node->op = PRIMITIVE; /* this is a primitive solid leaf node */
	leaf_node->leaf.p = p;
	leaf_node->parent = NULL;
	set_line(leaf_node);
	register_csg_node(leaf_node, op_solid);
}

@ The parameters for the sphere geometry are supplied in the following
format:

\smallskip

(``name" $radius$)

\smallskip

For instance, the specification {\tt ("Sphere A" 10.0)} initialises a
new solid sphere named ``Sphere A" located at (0.0, 0.0, 0.0) in the
world coordinate frame with a radius of 10.0. Again, we assume that
the radius is specified in the chosen unit of measurement, yet to be
discussed.

@<Read sphere geometry@>=
@<Create a new primitive solid@>;
@<Initialise primitive sphere with relevant data@>;
@<Register the primitive solid@>;

@ @<Initialise primitive sphere with relevant data@>=
read_count = fscanf(f, "(\"%[^\"]\" %lf)",
       op_solid, &p->s.radius);
if (EOF == read_count || 2 != read_count) {
        @<Exit after cleanup: failed to read from file@>;
}
p->type = SPHERE;
++(nodes_repo->stat[PRIMITIVE]);


@ The parameters for the cylinder geometry are supplied in the following
format:

\smallskip

(``name" $radius$ $height$)

\smallskip

For instance, the specification {\tt ("Cylinder A" 10.0 20.0)}
initialises a new solid cylinder named ``Cylinder A" 
located at (0.0, 0.0, 0.0) in the world coordinate frame with
base radius 10.0 and height 20.0.

@<Read cylinder geometry@>=
@<Create a new primitive solid@>;
@<Initialise primitive cylinder with relevant data@>;
@<Register the primitive solid@>;

@ @<Initialise primitive cylinder with relevant data@>=
read_count = fscanf(f, "(\"%[^\"]\" %lf %lf)",
       op_solid, &p->c.radius, &p->c.height);
if (EOF == read_count || 3 != read_count) {
        @<Exit after cleanup: failed to read from file@>;
}
p->type = CYLINDER;
@<Prepare cylinder for containment testing@>;
++(nodes_repo->stat[PRIMITIVE]);

@ @<Prepare cylinder for containment testing@>=
@<Half the height the cylinder@>;
@<Calculate containment range for the cylinder@>;

@ When checking containment, we require half of the cylinder
|height|, because the origin of the cylinder is given by its
centroid. Hence, although the |height| were specified in full, we
shall store their halves internally.

@<Half the height the cylinder@>=
p->c.height /= 2.0;

 
@ The parameters for the torus geometry are supplied in the following
format:

\smallskip

(``name" $phi$ $phi\_start$ $theta$ $theta\_start$ $major$ $minor$)

\smallskip

For instance, the specification {\tt ("Torus A" 180.0 45.0 45.0 15.0
10.0 2.0)} initialises a partial torus named ``Torus A" located at
(0.0, 0.0, 0.0) in the world coordinate frame with major radius 10.0
and minor 2.0. Its start and end cross sections subtend an angle of
$180^{\circ}$, and begins it rotation at $45^{\circ}$ from the
positive $x$-axis. Furthermore, the tube subtends and angle of
$45^{\circ}$, and begins its rotation about the center of the tube
starting at $15^{\circ}$ from the $xz$-plane.

Note here that, although mathematically the parametric angles $\phi,
\theta < 2\pi$ are defined as radians, their computational values
|phi| and |theta| (and their start angles |phi_start| and
|theta_start|) are supplied in degrees. This avoids conversion 
between radians and degrees while carrying out containment
testing. Furthermore, we allow variables |phi| and |theta| to take
values of $360^{\circ}$, which signifie complete rotations.

@<Read torus geometry@>=
@<Create a new primitive solid@>;
@<Initialise primitive torus with relevant data@>;
@<Register the primitive solid@>;

@ @<Initialise primitive torus with relevant data@>=
read_count = fscanf(f, "(\"%[^\"]\" %lf %lf %lf %lf %lf %lf)",
       op_solid, &p->t.phi, &p->t.phi_start, &p->t.theta,
       &p->t.theta_start, &p->t.major, &p->t.minor);
if (EOF == read_count || 7 != read_count) {
        @<Exit after cleanup: failed to read from file@>;
}
p->type = TORUS;
@<Prepare torus for containment testing@>;
++(nodes_repo->stat[PRIMITIVE]);

@ @<Prepare torus for containment testing@>=
@<Calculate radial containment range for the torus@>;
@<Calculate end angles for the torus@>;

@*2 Union.

The union of two solids is defined as the volume which consist of
points that are {\bf inside either of} these solids. In the geometry
input file, the union of two solids are defined using the following
format: $$\vcenter{(``target" ``left solid" ``right solid")}$$

Here, ``left solid" and ``right solid" are the names of the
union operands, and the result of the union is to be stored using the
name ``target". For the union command to be valid, no solid with the
name ``target" must already exists within the system. Whereas, both
operands must already exists within the system and could represent
primitive solids, or intermediate solids that are defined by a CSG
sub-tree. They are not required to share space (i.e., the solids may
be detached from one another). Because the union operator is {\it
commutative}@^commutative operator@>, the order of the operands is
unimportant. However, for performance considerations, it is advisable
to order the solids so that the solid on the left requires the least
amount of computation to determine inclusion (i.e., determination of
whether a point exists inside the solid).

For instance, the union specification {\tt ("U1" "Cylinder A" "Torus
A")} finds the union of two solids named ``Cylinder A" and ``Torus A"
and stores the result using the name ``U1".

@ @<Read union operation@>=
read_count = fscanf(f, "(\"%[^\"]\" \"%[^\"]\" \"%[^\"]\")",
	   op_target, op_left, op_right);
if (EOF == read_count || 3 != read_count)
        @<Exit after cleanup: failed to read from file@>;
@<Check that the target does not already exists@>;
@<Create union operator node@>;

@ @<Check that the target does not already exists@>=
target_solid = find_csg_node(op_target);
if (NULL != target_solid) {
        solid_name = op_target;
        @<Exit after cleanup: solid already exists@>;
}

@ @<Exit after cleanup: solid already exists@>=
goto solid_exists_exit_after_cleanup;

@ We are now ready to create a union operator. But first, we must
ensure that the solids specified as the operands already exists within
the system. If they are, we create a new operator node and make its
left and right subtrees point to these existing solids.

@<Create union operator node@>=
@<Find solid that corresponds to the left-hand operand@>;
@<Find solid that corresponds to the right-hand operand@>;
@<Create new union operator node@>;

@ @<Find solid that corresponds to the left-hand operand@>=
left_solid = find_csg_node(op_left);
if (NULL == left_solid) {
        solid_name = op_left;
        @<Exit after cleanup: solid does not exists@>;
}

@ @<Exit after cleanup: solid does not exists@>=
goto no_solid_exists_exit_after_cleanup;

@ @<Find solid that corresponds to the right-hand operand@>=
right_solid = find_csg_node(op_right);
if (NULL == right_solid) {
        solid_name = op_right;
        @<Exit after cleanup: solid does not exists@>;
}

@ When a new union operator node is created, we are actually creating
a new {\sl intermediate solid}@^intermediate solid@>. Hence, this new
solid must be registered with the system, so that they may be used by
later commands. 

@<Create new union operator node@>=
internal_node = create_csg_node();
if (NULL == internal_node) {
        @<Exit after cleanup: failed to create internal node@>;
} else {
        internal_node->op = UNION; /* this is an internal node */
	@<Set line number of node and usage bits of children@>;
	@<Set parents, and pointers to intermediate solids@>;
	register_csg_node(internal_node, op_target); /* register operator
	node using the target name */
	++(nodes_repo->stat[UNION]);
}

@ @<Set line number of node and usage bits of children@>=
set_line(internal_node);
set_ubit(left_solid);
set_ubit(right_solid);

@ @<Set parents, and pointers to intermediate solids@>=
internal_node->internal.left = left_solid;
internal_node->internal.right = right_solid;
internal_node->parent = NULL;
left_solid->parent = internal_node;
right_solid->parent = internal_node;

@ @<Exit after cleanup: failed to create internal node@>=
goto create_operator_failed_exit_after_cleanup;

@*2 Intersection.

The intersection of two solids is defined as the volume which consist of
points that are {\bf both inside} these solids. In the geometry input file,
the intersection is specified using the same format as that of union:
$$\vcenter{(``target" ``left solid" ``right solid")}$$

Here, ``left solid" and ``right solid" are the names of the
intersection operands, and the result of the intersection is to be
stored using the name ``target". For the intersection command to be
valid, no solid with the name ``target" must already exists within the
system. Whereas, both operands must already exists within the system
and could represent primitive solids, or intermediate solids that are
defined by a CSG sub-tree. Because the intersection operator is {\it
commutative}@^commutative operator@>, the order of the operands is
unimportant.

For instance, the intersection specification {\tt ("I1" "Cylinder A" "Torus
A")} finds the intersection of two solids named ``Cylinder A" and ``Torus A"
and stores the result using the name ``I1".

@ @<Read intersection operation@>=
read_count = fscanf(f, "(\"%[^\"]\" \"%[^\"]\" \"%[^\"]\")",
	   op_target, op_left, op_right);
if (EOF == read_count || 3 != read_count)
        @<Exit after cleanup: failed to read from file@>;
@<Check that the target does not already exists@>;
@<Create intersection operator node@>;

@ We are now ready to create a intersection operator. But first, we
must ensure that the solids specified as the operands already exists
within the system. If they are, we create a new operator node and make
its left and right subtrees point to these existing solids.

@<Create intersection operator node@>=
@<Find solid that corresponds to the left-hand operand@>;
@<Find solid that corresponds to the right-hand operand@>;
@<Create new intersection operator node@>;

@ When a new intersection operator node is created, we are actually
creating a new {\sl intermediate solid}@^intermediate solid@>. Hence,
this new solid must be registered with the system, so that they may be
used by later commands. 

@<Create new intersection operator node@>=
internal_node = create_csg_node();
if (NULL == internal_node) {
        @<Exit after cleanup: failed to create internal node@>;
} else {
        internal_node->op = INTERSECTION; /* this is an internal node */
	@<Set line number of node and usage bits of children@>;
        @<Set parents, and pointers to intermediate solids@>;
	register_csg_node(internal_node, op_target); /* register operator
	node using the target name */
	++(nodes_repo->stat[INTERSECTION]);
}

@*2 Difference.

The difference between two solids is defined as the volume which
consist of points that are {\bf inside the first} solid, but {\bf not
inside the second}. In the geometry input file, the difference
is specified using the same format as that of union and intersection:
$$\vcenter{(``target" ``left solid" ``right solid")}$$

Here, ``left solid" and ``right solid" are the names of the operands,
and the difference is to be stored using the name ``target". For the
difference command to be valid, no solid with the name ``target" must
already exists within the system. Whereas, both operands must already
exists within the system and could represent primitive solids, or
intermediate solids that are defined by a CSG sub-tree. Because the
difference operator is {\sl non-commutative}@^non-commutative operator@>,
the order of the operands is important. Hence, the
difference consists of all the points that are inside the left-hand
operand, but not inside the right-hand operand.

For instance, the difference specification {\tt ("D1" "Cylinder A"
"Torus A")} finds the difference by subtracting ``Torus A" from
``Cylinder A", and stores the result using the name ``D1".

@ @<Read difference operation@>=
read_count = fscanf(f, "(\"%[^\"]\" \"%[^\"]\" \"%[^\"]\")",
	   op_target, op_left, op_right);
if (EOF == read_count || 3 != read_count)
        @<Exit after cleanup: failed to read from file@>;
@<Check that the target does not already exists@>;
@<Create difference operator node@>;

@ We are now ready to create a difference operator. But first, we
must ensure that the solids specified as the operands already exists
within the system. If they are, we create a new operator node and make
its left and right subtrees point to these existing solids.

@<Create difference operator node@>=
@<Find solid that corresponds to the left-hand operand@>;
@<Find solid that corresponds to the right-hand operand@>;
@<Create new difference operator node@>;

@ When a new difference operator node is created, we are actually
creating a new {\sl intermediate solid}@^intermediate solid@>. Hence,
this new solid must be registered with the system, so that they may be
used by later commands. 

@<Create new difference operator node@>=
internal_node = create_csg_node();
if (NULL == internal_node) {
        @<Exit after cleanup: failed to create internal node@>;
} else {
        internal_node->op = DIFFERENCE; /* this is an internal node */
	@<Set line number of node and usage bits of children@>;
	@<Set parents, and pointers to intermediate solids@>;
	register_csg_node(internal_node, op_target); /* register operator
	node using the target name */
	++(nodes_repo->stat[DIFFERENCE]);
}

@*2 Translation.

Translation relocates a solid by displacing the solid's origin in the
$x$, $y$ and $z$ axes with respect to the world coordinate frame. The
paramaters for a translation operator is specified as a |displacement|
vector. The unit of displacement depends on the unit chosen by the
user when specifying the CSG tree; however, this unit does not matter
within the system since all measurement units are normalised
internally, so that all measurements in a given category can be
combined readily with other measurements in the same category.

@<Translation parameters@>=
struct {
       Vector displacement; /* displacement of solid's origin */
} t;

@ @<Read translation operation@>=
@<Read translation parameters@>;
@<Find the target solid for the operation@>;
@<Create translation parameter node@>; 
@<Create and register translation operator node@>;

@ The translation parameters are specified in the geometry input file
using the following format:

\smallskip

(``target'' ``solid'' $dx$ $dy$ $dz$)

\smallskip

\noindent where, ``target'' is the name of the target solid which is
obtained by translating ``solid''. The lengths $dx$, $dy$ and $dz$ are
the respective displacements along the $x$, $y$ and $z$ axes in the
world coordinate frame. For instance, the translation specification
{\tt ("T1" "Solid" 10.0 50.0 20.0)} means, translate the solid
associated with the name ``Solid'' by displacing its origin by 10.0
units along the $x$ axis, 50.0 units along the $y$ axis, and 20.0
units along the $z$ axis, and register this intermediate solid as
``T1''.

@<Read translation parameters@>=
read_count = fscanf(f, "(\"%[^\"]\" \"%[^\"]\" %lf %lf %lf)",
	   op_target, op_solid, &op_x, &op_y, &op_z);
if (EOF == read_count || 5 != read_count)
        @<Exit after cleanup: failed to read from file@>;


@ In order for a translation to be applicable, the supplied target
solid must already exists within the system, either by prior
definition as a primitive solid, or as an intermediate solid
defined by a CSG subtree.

@<Find the target solid for the operation@>=
target_solid = find_csg_node(op_solid);
if (NULL == target_solid) {
        solid_name = op_solid;
        @<Exit after cleanup: solid does not exists@>;
}

@ We now have a valid translation, so create the parameter leaf
node which will become the right-child of the translation operator
node.

@<Create translation parameter node@>=
leaf_node = create_csg_node();
if (NULL == leaf_node) {
        @<Exit after cleanup: failed to create leaf node@>;
} else {
        leaf_node->op = PARAMETER; /* this is a parameter leaf node */
	leaf_node->leaf.t.displacement[0] = op_x;
	leaf_node->leaf.t.displacement[1] = op_y;
	leaf_node->leaf.t.displacement[2] = op_z;
        @<Set up the matrix for translation@>;
}

@ @<Exit after cleanup: failed to create leaf node@>=
goto create_parameter_failed_exit_after_cleanup;

@ Finally, create the translation operator node, and attach the target
solid and the parameter node.

@<Create and register translation operator node@>=
internal_node = create_csg_node();
if (NULL == internal_node) {
        @<Exit after cleanup: failed to create internal node@>;
} else {
        internal_node->op = TRANSLATE; /* this is an operator internal node */
	@<Set line number of node and usage bit of transformed node@>;
        @<Set target, parameters and parent pointers@>;
	register_csg_node(internal_node, op_target);
	++(nodes_repo->stat[TRANSLATE]);
}

@ @<Set line number of node and usage bit of transformed node@>=
set_line(internal_node);
set_ubit(target_solid);

@ @<Set target, parameters and parent pointers@>=
internal_node->internal.left = target_solid;
internal_node->internal.right = leaf_node;
internal_node->parent = NULL;
target_solid->parent = internal_node;
leaf_node->parent = internal_node;


@*2 Rotation.

We rotate a solid relative to an axis. The angle of rotation
$\theta$ is specified in {\sl radians}@^radians@>. If $\theta < 0$, we
have {\sl clockwise rotation}@^clockwise rotation@>; if $\theta > 0$,
we have {\sl counter-clockwise rotation}@^counter-clockwise rotation@>.
When $\theta = 0$, no rotation is applied. In order to
apply a rotation, we also need an {\sl axis of rotation}@^axis of rotation@>.
This is specified as a unit vector@^unit vector@> named |axis|, which
is defined relative to the origin of the world coordinate frame.

@<Rotation parameters@>=
struct {
       double theta; /* angle of rotation in radians */
       Vector axis; /* axis of rotation (unit vector from origin of world
       coordinate frame) */
} r;

@ @<Read rotation operation@>=
@<Read rotation parameters@>;
@<Find the target solid for the operation@>;
@<Create rotation parameter node@>; 
@<Create and register rotation operator node@>;

@ The rotation parameters are specified in the geometry input file
using the following format:

\smallskip

(``target'' ``solid'' $\theta$ $ux$ $uy$ $uz$)

\smallskip

\noindent where, ``target'' is the name of the solid which is obtained
by rotating ``solid'' by an angle of $\theta$ radians relative to the
axis defined by the unit vector with components $ux$, $uy$ and $uz$
respective along the $x$, $y$ and $z$ axes. For instance, the
rotation specification {\tt ("R1" "Solid" 90.0 1.0 0.0 0.0)} means,
rotate the solid associated with the name ``Solid'' by 90.0 degree
radians relative to the $x$-axis in world coordinate frame and
register this intermediate solid as ``R1''.

@<Read rotation parameters@>=
read_count = fscanf(f, "(\"%[^\"]\" \"%[^\"]\" %lf %lf %lf %lf)",
	   op_target, op_solid, &op_theta, &op_x, &op_y, &op_z);
if (EOF == read_count || 6 != read_count)
        @<Exit after cleanup: failed to read from file@>;

@ In order for a rotation to be applicable, the supplied target
solid must already exists within the system, either by prior
definition as a primitive solid, or as an intermediate solid
defined by a CSG subtree.

@<Find the target solid for the operation@>=
target_solid = find_csg_node(op_solid);
if (NULL == target_solid) {
        solid_name = op_solid;
        @<Exit after cleanup: solid does not exists@>;
}

@ We now have a valid rotation, so create the parameter leaf node.

@<Create rotation parameter node@>=
leaf_node = create_csg_node();
if (NULL == leaf_node) {
        @<Exit after cleanup: failed to create leaf node@>;
} else {
        leaf_node->op = PARAMETER; /* this is a parameter leaf node */
	leaf_node->leaf.r.theta = op_theta; /* angle of rotation */
	leaf_node->leaf.r.axis[0] = op_x;
	leaf_node->leaf.r.axis[1] = op_y;
	leaf_node->leaf.r.axis[2] = op_z;
        @<Set up the matrix for rotation@>;
}

@ Finally, create the rotation operator node, and attach the target
solid and the parameter node.

@<Create and register rotation operator node@>=
internal_node = create_csg_node();
if (NULL == internal_node) {
        @<Exit after cleanup: failed to create internal node@>;
} else {
        internal_node->op = ROTATE; /* this is an operator internal node */
	@<Set line number of node and usage bit of transformed node@>;
        @<Set target, parameters and parent pointers@>;
	register_csg_node(internal_node, op_target);
	++(nodes_repo->stat[ROTATE]);
}


@*2 Scaling.

Scaling increases or decreases the volume of a solid. The amount of
scaling applied with respect to each of the axes are specified using a
|scale| vector. If the scaling factors are all zero, scaling is not
applied. A scaling transformation maintains the shape and the form of
the solid if and only if the scaling factor along all of the axes are
equal. This is referred to as {\sl uniform scaling}@^uniform scaling@>.

@<Scaling parameters@>=
struct {
       Vector scale; /* scaling factor */
} s;

@ @<Read scaling operation@>=
@<Read scaling parameters@>;
@<Find the target solid for the operation@>;
@<Create scaling parameter node@>; 
@<Create and register scaling operator node@>;

@ The scaling parameters are specified in the geometry input file
using the following format:

\smallskip

(``target" ``solid'' $sx$ $sy$ $sz$)

\smallskip

\noindent where, ``target'' is the name of the solid which is obtained
after scaling ``solid'' by the scaling factors $sx$, $sy$ and $sz$
along the $x$, $y$ and $z$ axes, respectively. For instance, the
scaling specification {\tt ("S1" "Solid" 1.0 2.0 3.0)} means, increase
the scale of the solid associated with the name ``Solid'' by a factor
of 2.0 and 3.0 respectively along the $y$ and $z$ axes in world
coordinate frame and register this intermediate solid as ``S1''. 

@<Read scaling parameters@>=
read_count = fscanf(f, "(\"%[^\"]\" \"%[^\"]\" %lf %lf %lf)",
	   op_target, op_solid, &op_x, &op_y, &op_z);
if (EOF == read_count || 5 != read_count)
        @<Exit after cleanup: failed to read from file@>;

@ In order for a scaling to be applicable, the supplied target
solid must already exists within the system, either by prior
definition as a primitive solid, or as an intermediate solid
defined by a CSG subtree.

@<Find the target solid for the operation@>=
target_solid = find_csg_node(op_solid);
if (NULL == target_solid) {
        solid_name = op_solid;
        @<Exit after cleanup: solid does not exists@>;
}

@ We now have a valid scaling, so create the parameter leaf node.

@<Create scaling parameter node@>=
leaf_node = create_csg_node();
if (NULL == leaf_node) {
        @<Exit after cleanup: failed to create leaf node@>;
} else {
        leaf_node->op = PARAMETER; /* this is a parameter leaf node */
	leaf_node->leaf.s.scale[0] = op_x;
	leaf_node->leaf.s.scale[1] = op_y;
	leaf_node->leaf.s.scale[2] = op_z;
	@<Set up the matrix for scaling@>;
}

@ Finally, create the scaling operator node, and attach the target
solid and the parameter node.

@<Create and register scaling operator node@>=
internal_node = create_csg_node();
if (NULL == internal_node) {
        @<Exit after cleanup: failed to create internal node@>;
} else {
        internal_node->op = SCALE; /* this is an operator internal
	node */
	@<Set line number of node and usage bit of transformed node@>;
        @<Set target, parameters and parent pointers@>;
	register_csg_node(internal_node, op_target);
	++(nodes_repo->stat[SCALE]);
}

@*2 Registering a solid.
A CSG tree defines a solid by specifying how primitive solids are
combined using transformation, translation, and boolean
operators. The solid thus defined is represented by the root of the
CSG tree, and any test against the solid must begin at this
root. Hence, if we define $n$ solids in the simulation world,
each of these solids must be represented by $n$ CSG trees. We maintain
a list of CSG solids using a table of pointers to CSG root nodes.

@<Read solid registration@>=
@<Read target solid for registration@>;
@<Find the target solid for the operation@>;
@<Register the target solid@>; 

@ To separate solids from intermediate solid components, we must
explicitly register each solid (represented by a node in the single
CSG tree) using the following format:

\smallskip

(``target'')

\smallskip

\noindent where, ``target'' is the name of the solid to register.

@<Read target solid for registration@>=
read_count = fscanf(f, "(\"%[^\"]\")", op_solid);
if (EOF == read_count || 1 != read_count)
        @<Exit after cleanup: failed to read from file@>;

@ @<Register the target solid@>=
process_and_register_solid(target_solid);
set_ubit(target_solid);

@*2 Defining the simulation world.
The simulation world defines the volume of space that we are
interested in. Only particles that are inside this volume are tracked
during the simulation. Logically, the simulation world is represented
by a bounding box, which is axis-aligned with the world coordinate frame.

@<Global variables@>=
BoundingBox sim_world = {ZERO_VECTOR, ZERO_VECTOR};
uint32_t num_subcuboids = 0;
uint32_t div_subcuboids[3] = {1, 1, 1}; /* division along $x$, $y$ and $z$ */

@ In the geometry input file, the specification of the simulation
world uses the following format:

\smallskip

($lx$ $ly$ $lz$ $ux$ $uy$ $uz$ $l$ $m$ $n$)

\smallskip

\noindent where, the coordinates $(lx, ly, lz)$ and
$(ux, uy, uz)$ respectively give the lower and upper bounds. The
values $l$, $m$ and $n$ gives the number of divisions along the
$x$, $y$ and $z$ axes, so that the simulation world is divided into $l
\times m \times n$ subcuboids. If there are multiple specifications of
the simulation world, the last one read will set the effective bounds.

@<Read simulation world specification@>=
read_count = fscanf(f, "(%lf %lf %lf %lf %lf %lf %u %u %u)",
    &sim_world.l[0], &sim_world.l[1], &sim_world.l[2],
    &sim_world.u[0], &sim_world.u[1], &sim_world.u[2],
    &div_subcuboids[0], &div_subcuboids[1], &div_subcuboids[2]);
if (EOF == read_count || 9 != read_count)
    @<Exit after cleanup: failed to read from file@>;
num_subcuboids = div_subcuboids[0] * div_subcuboids[1] * div_subcuboids[2];
if (num_subcuboids > MAX_SUBCUBOIDS)
    @<Exit after cleanup: invalid division of simulation world@>;

@ @<Exit after cleanup: invalid division of simulation world@>=
goto invalid_subcuboid_division_exit_after_cleanup;

@ Function |print_sim_world(f)| prints the lower and upper bounds of
the simulation world to the I/O stream pointed to by |f|.
@<Global functions@>=
void print_sim_world(FILE *f)
{
    fprintf(f, "Simulation world: [%lf, %lf, %lf : %lf, %lf, %lf]\n",
    sim_world.l[0], sim_world.l[1], sim_world.l[2],
    sim_world.u[0], sim_world.u[1], sim_world.u[2]);
}

@*2 Bounding Boxes.
Every node in the CSG tree is enclosed by a bounding box. This
allows us to accelerate the search for solids by efficiently pruning
the search trees.

@<Information common to all CSG nodes@>=
BoundingBox bb;

@ Function |primitive_bb(p, bb)| calculates the bounding box of a
primitive solid |p| and stores the result in the supplied bounding box
variable |bb|. This function returns |true| if the bounding box was
calculated successfully; otherwise, |false| is returned.

@<Global functions@>=
bool primitive_bb(Primitive *p, BoundingBox *bb)
{
    switch(p->type) {
    case BLOCK: @<Calculate bounding box of primitive block@>;@+ break;
    case SPHERE: @<Calculate bounding box of primitive sphere@>;@+ break;
    case CYLINDER: @<Calculate bounding box of primitive cylinder@>;@+ break;
    case TORUS: @<Calculate bounding box of primitive torus@>;@+ break;
    default: return false; /* invalid primitive */
    }
    return true;
}

@ Note here that the |length|, |height| and |width| are half-lengths
of the block's dimension along the $x$, $y$ and $z$ axes. Furthermore,
the block is centered at the origin of the world coordinate frame, so
that it's centroid coincides with the origin.

@<Calculate bounding box of primitive block@>=
bb->l[0] = -p->b.length;
bb->u[0] = p->b.length;
bb->l[1] = -p->b.height;
bb->u[1] = p->b.height;
bb->l[2] = -p->b.width;
bb->u[2] = p->b.width;

@ The sphere is centered at the origin of the world coordinate frame.
@<Calculate bounding box of primitive sphere@>=
bb->l[0] = bb->l[1] = bb->l[2] = -p->s.radius;
bb->u[0] = bb->u[1] = bb->u[2] = p->s.radius;

@ The centroid of the cylinder is centered at the origin of the world
coordinate frame, and the normals at the center of the two circular
faces of the cylinder are parallel to the $y$-axis.
@<Calculate bounding box of primitive cylinder@>=
bb->l[0] = bb->l[2] = -p->c.radius;
bb->u[0] = bb->u[2] = p->c.radius;
bb->l[1] = -p->c.height;
bb->u[1] = p->c.height;

@ The torus is centered at the origin so that its center coincides
with the origin of the world coordinate frame. Furthermore, the radial
surface of the torus lies on the $xz$-plane.
@<Calculate bounding box of primitive torus@>=
bb->l[0] = bb->l[2] = -(p->t.major + p->t.minor);
bb->u[0] = bb->u[2] = p->t.major + p->t.minor;
bb->l[1] = -p->t.minor;
bb->u[1] = p->t.minor;

@ Function |intersection_bb(n, l, r, a)| calculates the bounding
box |n| of the intersection node along the supplied axis |a| using the
bounding boxes of the left and right nodes, |l| and |r|
respectively. The parameter |a| must only take values 0, 1, and 2,
representing respectively the $x$, $y$ and $z$ axes.

@d X_AXIS 0
@d Y_AXIS 1
@d Z_AXIS 2
@<Global functions@>=
bool intersection_bb(BoundingBox *n, BoundingBox *l, BoundingBox *r, int a)
{
    if (l->l[a] < r->l[a] && l->u[a] > r->u[a]) { /* left fully encloses right */
        n->l[a] = l->l[a];
        n->u[a] = l->u[a];
    } else if (r->l[a] < l->l[a] && r->u[a] > l->u[a]) { /* right fully encloses left */
        n->l[a] = r->l[a];
        n->u[a] = r->u[a];
    } else {
        if (l->l[a] < r->l[a]) { /* intersection with left behind
	    right */
            n->l[a] = r->l[a];
            n->u[a] = l->u[a];
	} else { /* intersection with right behind
	    left */
            n->l[a] = l->l[a];
            n->u[a] = r->u[a];
	}
    }
    return true;
}

@ Function |apply_affine_matrix(m,v,r)| applies the affine
transformation matrix |m| to the vector |v| and stores the transformed
vector in |r|. The original vector |v| is left unmodified. The two
macros |affine_normal()| and |affine_inverse()| use this function to
respectively apply the affine matrix for normal transformations, or
the inverse transformations ({\sl matrix inverse} of the normal
affine matrix).

@d affine_normal(n,v,r) apply_affine_matrix((n)->affine, (v), (r))
@d affine_inverse(n,v,r) apply_affine_matrix((n)->inverse, (v), (r))
@<Global functions@>=
void apply_affine_matrix(Matrix m, Vector v, Vector r)
{
    r[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2] + m[0][3];
    r[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2] + m[1][3];
    r[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2] + m[2][3];
    r[3] = 1.0;
}

@ Function |calculate_bounding_box(n)| calculates the bounding box
of the solid represented by the supplied root of a CSG tree. To
calculate the bounding box at a given node, the function uses the
bounding boxes of the left and right subtrees. Only at the leaves,
where we reach a primitive solid, we calculate the bounding box using
only the definition of the primitive solid. This function returns
|true| if the bounding box was calculated successfully; otherwise,
|false| is returned.

@<Global functions@>=
bool calculate_bounding_box(CSG_Node *n)
{
	CSG_Node *l, *r; /* left and right subtrees */
	Vector temp, c[8]; /* for affine transformation */
	int i, j;
	if (is_primitive(n)) {
	    if (primitive_bb(n->leaf.p, &n->bb)) {
	        @<Calculate affine transformed bounding box@>;
	        return true;
            } else return false;
        }
	if (!calculate_bounding_box(n->internal.left)) return false;
	if (!calculate_bounding_box(n->internal.right)) return false;
	@<Find bounding box for the node using left and right subtrees@>;
	return true;
}

@ The bounding box of a primitive solid must be defined using the
actual location and orientation of the primitive in the world
coordinate frame. This means that, after the bounding box for the
primitive has been calculated under the assumption that the center of
the primitive coincides with the origin of the world coordinate frame,
we must apply the affine transformation for that primitive to the
calculated bounding box, to obtain the actual bounds.

@<Calculate affine transformed bounding box@>=
@<Calculate the eight corners of the bounding box@>;
@<Apply affine transformations to the corners@>;
@<Calculate axis aligned bounding box of the transformed bounding box@>;

@ @<Calculate the eight corners of the bounding box@>=
c[0][0] = c[1][0] = c[2][0] = c[3][0] = n->bb.l[0];
c[0][1] = c[1][1] = c[4][1] = c[5][1] = n->bb.l[1];
c[0][2] = c[2][2] = c[4][2] = c[6][2] = n->bb.l[2];
c[4][0] = c[5][0] = c[6][0] = c[7][0] = n->bb.u[0];
c[2][1] = c[3][1] = c[6][1] = c[7][1] = n->bb.u[1];
c[1][2] = c[3][2] = c[5][2] = c[7][2] = n->bb.u[2];

@ @<Apply affine transformations to the corners@>=
for (i = 0; i < 8; ++i) {
        c[i][3] = 1.0; /* homogenise the corner vector */
        affine_normal(n, &c[i][0], temp);
        vector_copy(&c[i][0], temp);
}

@ After applying an affine transformation, a bounding box may no
longer be axis-aligned. Since we require axis-aligned bounding boxes,
we must re-calculate an axis-aligned bounding box of the transformed 
bounding box by using its corners. 
@<Calculate axis aligned bounding box of the transformed bounding box@>=
for (j = 0; j < 3; ++j) {
        n->bb.l[j] = c[0][j];
        n->bb.u[j] = c[0][j];
        for (i = 1; i < 8; ++i) {
                if (c[i][j] < n->bb.l[j]) n->bb.l[j] = c[i][j]; /*
		find minimum */
                if (c[i][j] > n->bb.u[j]) n->bb.u[j] = c[i][j]; /*
		find maximum */
	}
}
for (i = 0; i < 8; ++i) n->bb.l[3] = n->bb.u[3] = 1.0; /* homogenise bound vectors*/

@ @<Find bounding box for the node using left and right subtrees@>=
l = n->internal.left;
r = n->internal.right;
switch(BIT_MASK_NODE & n->op) {
case UNION: @<Calculate bounding box of union@>;@+ break;
case INTERSECTION: @<Calculate bounding box of intersection@>;@+ break;
case DIFFERENCE: @<Calculate bounding box of difference@>;@+ break;
default: return false;
}

@ In each of the axes, the lowest of the values stored in the left and
right subtree nodes becomes the lower bound for the node. Similarly,
the highest value becomes the upper bound. The resulting bounding box
must enclose both solids on the left and right subtrees.

@<Calculate bounding box of union@>=
n->bb.l[0] = (l->bb.l[0] < r->bb.l[0]) ? l->bb.l[0] : r->bb.l[0];
n->bb.l[1] = (l->bb.l[1] < r->bb.l[1]) ? l->bb.l[1] : r->bb.l[1];
n->bb.l[2] = (l->bb.l[2] < r->bb.l[2]) ? l->bb.l[2] : r->bb.l[2];
n->bb.u[0] = (l->bb.u[0] > r->bb.u[0]) ? l->bb.u[0] : r->bb.u[0];
n->bb.u[1] = (l->bb.u[1] > r->bb.u[1]) ? l->bb.u[1] : r->bb.u[1];
n->bb.u[2] = (l->bb.u[2] > r->bb.u[2]) ? l->bb.u[2] : r->bb.u[2];

@ When the left and right subtrees do not intersect, a bounding box
must not be defined for the node. This is represented simply by making
the upper bound smaller than the lower bound---we only need to do this
for only one axis; here, we choose the $x$-axis.

@d no_intersection_bb(left,right) ((left).u[0] < (right).l[0] ||
    (left).u[1] < (right).l[1] ||
    (left).u[2] < (right).l[2] ||
    (right).u[0] < (left).l[0] ||
    (right).u[1] < (left).l[1] ||
    (right).u[2] < (left).l[2])
@<Calculate bounding box of intersection@>=
if (no_intersection_bb(l->bb,r->bb)) {
        n->bb.l[0] = 1;
        n->bb.u[0] = -1;
} else {
    intersection_bb(&n->bb, &l->bb, &r->bb, X_AXIS);
    intersection_bb(&n->bb, &l->bb, &r->bb, Y_AXIS);
    intersection_bb(&n->bb, &l->bb, &r->bb, Z_AXIS);
}

@ For a solid defined by a boolean difference, the bounding box of the
node should be the bounding box of the left subtree before we subtract
the right subtree.

@<Calculate bounding box of difference@>=
n->bb = l->bb;

@*2 Merge affine transformations.
Every node has an affine transformation matrix, which is initialised
to the identity matrix. When an affine transformation $T$ is applied
to a node, the affine matrix for that node is updated. This affine
matrix |affine| gives the accumulated transformations starting at the
node to the root of the CSG tree. We also sometimes required the
affine matrix |inverse| for applying the transformations in the
reverse order, starting at the root until we reach the node.

@<Information common to all CSG nodes@>=
Matrix affine; /* inverse affine transformations matrix */
Matrix inverse; /* inverse of |affine| matrix */

@ @<Initialise affine matrices to the identity matrix@>=
matrix_copy(slot->affine, identity_matrix);
matrix_copy(slot->inverse, identity_matrix);

@ While the solids are built bottom-up, CSG nodes are added to the
nodes repository. When an affine transformation node is added to this
repository, the transformation is only applied to one target node,
which is an existing node in the nodes repository. However, once the
CSG tree has been built, we want each node to store the accumulated
affine transformations starting at the root of the CSG tree. Hence, we
pre-process the CSG tree so that all of the affine transformations are
merged. Furthermore, after the affine transformations have been
merged, we do not require the affine transformation nodes and its
parameters. Hence, these nodes are then remove from the CSG tree.

\bigskip

\centerline{\epsfig{file=figures/affine-merge,scale=1.2}}

\bigskip

The merge begins by first merging transformations starting at the
current node. As long as the left-hand child node is an affine
transformation node, we accumulate the affine transformation into |t|,
and continue until we have reached a node that is not an affine
transformation node. At this node, we store the accumulated
affine transformation matrix from |t|, and detach the list of affine
transformation nodes that was merged. Then, we merge the node's
subtrees, if they are not leaf nodes. Finally, we the return the last
node found so that the parent can update its child pointers after the
merge (some nodes may have been removed).

@ Function |merge_affine(r)| merges the affine transformations on all
the paths starting at the root of the CSG tree |r| and stores the
accumulated affine transformation matrix in each of the nodes.

@<Global functions@>=
CSG_Node *merge_affine(CSG_Node *r)
{
	CSG_Node *t, *p; /* temporary node and parent node */
	Matrix mm, m = IDENTITY_MATRIX; /* accumulated affine matrix */
	if (NULL == r || is_parameter(r)) return r;
	if (r->parent) matrix_copy(m, r->parent->affine); /*
	initialise with parent's affine matrix */
	if (is_translate(r) || is_rotate(r) || is_scale(r)) {
	    @<Merge sequence of affine transformation nodes@>;
	}
        @<Detach accumulated sequence of affine transformations@>;
	@<Merge affine transformation sequences in subtrees@>;
	return r; /* return the first non-affine transformation node */
}

@ @<Merge sequence of affine transformation nodes@>=
p = r->parent; /* parent of the first affine transformation node to
merge */
do {
        t = r;
	matrix_multiply(m, t->internal.right->affine, mm); /* accumulate using parameter on the right */
	matrix_copy(m, mm);
	r = t->internal.left; /* only left points to the next non-parameter node */
} while (is_translate(r) || is_rotate(r) || is_scale(r));
r->parent = p; /* update parent for first non-affine node at the end
of the sequence */

@ Since affine merge is carried out only once for the entire CSG solid
during registration, we do not need to multiply the accumulated affine
transformation matrix |m| with the affine matrix |r->affine| currently
available in the node, before replacing the value of
|r->affine|. This is because, |r->affine| is always initialised to the
|IDENTITY_MATRIX| matrix when CSG nodes are added to the nodes repository.
@<Detach accumulated sequence of affine transformations@>=
matrix_copy(r->affine, m);
matrix_inverse(r->affine, m);
matrix_copy(r->inverse, m);

@ After the affine transformations in the subtrees have been merged,
we must update the left and right daughter nodes, which are returned by
the recursive subtree merges. It is important to do this because, if
the previous daughter nodes were affine transformations, they will now
be unavailable after the merge.
@<Merge affine transformation sequences in subtrees@>=
if (!is_primitive(r)) {
    if ((t = merge_affine(r->internal.left))) r->internal.left = t;
    if ((t = merge_affine(r->internal.right))) r->internal.right = t;
}

@ Function |print_csg_tree(t,l)| prints the CSG tree |t| using inorder
tree traversal. The value of |l| gives the level-of-indentation to use
to print the node in the current recursive call.

@<Global functions@>=
void print_csg_tree(CSG_Node *t, uint32_t l) {
        int i;
        if (NULL == t) return;
        if (is_primitive(t)) {
                @<Print primitive solid information@>;
                return;
        }
        if (is_parameter(t)) {
                @<Print affine transformation parameters@>;
                return;
        }
	print_csg_tree(t->internal.left, l + 1);
	@<Print intermediate node information@>;
        print_csg_tree(t->internal.right, l + 1);
}

@ @<Print primitive solid information@>=
@<Print indentation@>;
p = t->leaf.p;
switch(p->type) {
case BLOCK:@/
     printf("BLOCK (%u): \"%s\" %lf %lf %lf", get_line(t),
        t->name, p->b.length, p->b.width, p->b.height);
        break;
case SPHERE:@/
     printf("SPHERE (%u): \"%s\" %lf", get_line(t),
        t->name, p->s.radius);
        break;
case CYLINDER:@/
     printf("CYLINDER (%u): \"%s\" %lf %lf", get_line(t),
        t->name, p->c.radius, p->c.height);
        break;
case TORUS:@/
     printf("TORUS (%u): \"%s\" %lf %lf %lf %lf %lf %lf",
        get_line(t), t->name, p->t.phi, p->t.phi_start,
        p->t.theta, p->t.theta_start, p->t.major, p->t.minor);
        break;
default: printf("unknown");
}
@<Print bounding box information@>;
@<Print affine transformation matrices@>;

@ @<Print indentation@>=
for (i = 0; i < l; ++i) printf("\t");

@ @<Print bounding box information@>=
printf(" [%lf, %lf, %lf : %lf, %lf, %lf]\n",
        t->bb.l[0], t->bb.l[1], t->bb.l[2],
        t->bb.u[0], t->bb.u[1], t->bb.u[2]);

@ @<Print affine transformation matrices@>=
matrix_print(stdout, t->affine, 4, 4, l + 1);
printf("\n");
matrix_print(stdout, t->inverse, 4, 4, l + 1);

@ @<Print affine transformation parameters@>=
@<Print indentation@>;
switch(BIT_MASK_NODE & t->parent->op) {
case TRANSLATE:@/
        printf("displacement: (%lf, %lf, %lf)",
	t->leaf.t.displacement[0], t->leaf.t.displacement[1],
	t->leaf.t.displacement[2]);
        break;
case ROTATE:@/
        printf("angle: %lf, axis: (%lf, %lf, %lf)", t->leaf.r.theta,
	t->leaf.r.axis[0], t->leaf.r.axis[1], t->leaf.r.axis[2]);
        break;
case SCALE:@/
        printf("scaling factor: (%lf, %lf, %lf)", t->leaf.s.scale[0],
	t->leaf.s.scale[1], t->leaf.s.scale[2]);
        break;
default: printf("unknown");
}
printf("\n");

@ @<Print intermediate node information@>=
@<Print indentation@>;
switch(BIT_MASK_NODE & t->op) {
case UNION: printf("UNION (%u): %s", get_line(t), t->name);@+ break;
case INTERSECTION: printf("INTERSECTION (%u): %s", get_line(t), t->name);@+ break;
case DIFFERENCE: printf("DIFFERENCE (%u): %s", get_line(t), t->name);@+ break;
case TRANSLATE: printf("TRANSLATE (%u): %s", get_line(t), t->name);@+ break;
case ROTATE: printf("ROTATE (%u): %s", get_line(t), t->name);@+ break;
case SCALE: printf("SCALE (%u): %s", get_line(t), t->name);@+ break;
default: printf("unknown");
}
@<Print bounding box information@>;
@<Print affine transformation matrices@>;

@*2 Forest of solids. We record all of the solids by recording the
root of their CSG tree. Recording the forest of trees is useful for
debugging purposes, and also, while printing out verbose
information. Furthermore, generation of the geometry tables are
significantly simplified, as we shall see in later sections. We do not
use a hash table because all of the solids will be processed
sequentially, and there will be no searching.

@d MAX_CSG_SOLIDS 512
@d reset_forest() memset(&forest_of_solids, 0, sizeof(forest_of_solids))
@<Global variables@>=
struct {
    uint32_t n; /* number of solids in forest */
    CSG_Node *s[MAX_CSG_SOLIDS]; /* root of CSG tree */
} forest_of_solids;

@ Function |print_forest()| prints all of the solids in the forest.
@<Global functions@>=
void print_forest()
{
	uint32_t i, j;
	for (i = 0; i < forest_of_solids.n; ++i) {
	        if (NULL == forest_of_solids.s[i]) continue;
		printf("Solid %s:\n\n", forest_of_solids.s[i]->name);
	        print_csg_tree(forest_of_solids.s[i], 0);
		printf("\n");
		for (j = 0; j < 80; j++) printf("-");
		printf("\n");
	}
}

@ Function |process_and_register_solid(r)| registers the solid
pointed to by the CSG tree root |r|. Before registering the solid, the
CSG tree is processed to the required form by merging the affine
transformations, and by finding the bounding box.

@<Global functions@>=
void process_and_register_solid(CSG_Node *r)
{
	CSG_Node *t;
	if (NULL == r) return;
	t = merge_affine(r); /* returns a non-affine node;
	otherwise returns |NULL| */
	if (NULL == t) t = r; /* in case |r| was a primitive solid */
	calculate_bounding_box(t);
	forest_of_solids.s[(forest_of_solids.n)++] = t;
}

@*1 Particles inside solids.
During the simulation, the {\sl MCS} system must determine which
materials a particle is interacting with at each step of its
track. Since material properties are associated with solids inside the
simulation world, we must first determine the solid which contains the
particle at each step of the particle's trajectory. Once we know the
solid, we can retrieve the relevant material properties associated
with the solid, and then carry out the necessary physics processes.

@ To find the solid which contains a given particle, we require two
algorithms. Firstly, we need an algorithm that will decide if a
particle (i.e., point inside the three-dimensional simulation world)
is located inside a given solid. This algorithm will traverse the CSG
tree that defines the said solid. Secondly, we need an algorithm that
will provide us with a list of solids that could potentially contain the
particle. In the most basic form, we could do an exhaustive search
across all of the solids in the simulation world. However, to improve
efficiency, we must use a space-partitioning scheme to prune the
search-space.

@ A containment test must return one of the following:

$$\vcenter{\halign{\hfil # & #\hfil \cr
OUTSIDE & if the point is outside the solid,\cr
INSIDE & if the point is inside the solid,\cr
SURFACE & if point is on the surface, and\cr
INVALID & if either the solid or the vector is undefined.\cr
}}$$

@<Type definitions@>=
typedef enum {OUTSIDE = 0, INSIDE, SURFACE, INVALID} Containment;

@*2 Containment inside a solid block.
The three pairs of opposite faces of a block determines its
containment range on each of the three axes. Hence, we can test if a
vector is inside, outside, or on the surface of the block by testing
its $x$, $y$ and $z$ components against the corresponding range.

@<Information that defines a primitive block@>=
double x0, x1, y0, y1, z0, z1; /* containment range (i.e., bounding box) */

@ Note here that the dimensions |length|, |height| and |width|
correspond respectively to the $x$, $y$ and $z$ axes in world
coordinate frame, and that we are adding or subtracting half-lengths
of the respective dimensions.

@<Calculate containment range for the block@>=
p->b.x0 = - p->b.length;
p->b.x1 = p->b.length;
p->b.y0 = - p->b.height;
p->b.y1 = p->b.height;
p->b.z0 = - p->b.width;
p->b.z1 = p->b.width;

@ Let $(x, y, z)$ represent the components of the vector to be
tested. Also let the three intervals $[x_0, x_1]$, $[y_0,
y_1]$ and $[z_0, z_1]$ respectively define the containment range for
the block along the $x$, $y$ and $z$ axes in the world coordinate
frame. Then the point defined by the vector |v| is:

$$\vcenter{\halign{\hfil # & # \hfil \cr
outside the block if & $(x < x_0 \vee x > x_1) \vee (y < y_0 \vee y >
y_1) \vee (z < z_0 \vee z > z_1)$,\cr
inside the block if & $(x > x_0 \wedge x < x_1) \wedge (y >
y_0 \wedge y < y_1) \wedge (z > z_0 \wedge z < z_1)$, and\cr
on the surface, & otherwise.\cr
}}$$

In the following implementation, we carry out a surface test instead
of the inside test. This avoids boolean conjuctions. Furthermore, we
validate the surface test by doing an outside test first.

@<Global functions@>=
Containment is_inside_block(const Vector v, const Primitive *p)
{
        if (v[0] < p->b.x0 || v[0] > p->b.x1 || v[1] < p->b.y0 || v[1]
	> p->b.y1 || v[2] < p->b.z0 || v[2] > p->b.z1) 
                return OUTSIDE;
        if (v[0] == p->b.x0 || v[0] == p->b.x1 || v[1] == p->b.y0 || v[1]
	== p->b.y1 || v[2] == p->b.z0 || v[2] == p->b.z1)
                return SURFACE;
	return INSIDE;
}

@*2 Containment inside a solid sphere.
We determine the containment of a point inside a sphere by calculating
the distance of the point from the sphere's origin. Let $\delta$ be
the distance of the point from the origin of the sphere. Then the
point defined by the vector |v| is:

$$\vcenter{\halign{\hfil # & # \hfil \cr
outside the sphere if & $(\delta > radius)$,\cr
inside the sphere if & $(\delta < radius)$, and\cr
on the surface, & otherwise.\cr
}}$$

During containment testing, the origin of the sphere always coincides
with the origin of the world coordinate frame. Hence, to determine
$\delta$, we only need to calculate the magnitude of the vector
|v|.

Note here that we carry out the {\sl outside test} first. This is
because, in a real-world setup, a spherical component is less likely
to occupy a large volume of the simulation space. Thus, it is highly
likely that a particle is outside the sphere most of the time.

@<Global functions@>=
Containment is_inside_sphere(const Vector v, const Primitive *p)
{
	double delta = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	if (delta > p->s.radius) return OUTSIDE; /* highly likely */
	if (delta == p->s.radius) return SURFACE; /* least likely */
	return INSIDE;
}

@*2 Containment inside a solid cylinder.
During containment testing, the origin of the cylinder coincides with
the origin of the world coordinate frame; and the normals at the
center of the circular bases of the cylinder are parallel to the
$y$-axis. Hence, to test containment of a point inside a cylinder, we
first check if the $y$-component of vector |v| is within the range
defined by the two parallel circular faces of the cylinder. Secondly,
we check if the projected distance of the vector |v| from the origin
on the $xz$-plane is within the area subscribed on the same plane by
the circular surfaces of the cylinder.

@<Information that defines a primitive cylinder@>=
double y0, y1; /* containment range of cylinder height */

@ Note here that the dimension |height| corresponds to the $y$-axis
in world coordinate frame, and that we are adding or subtracting
half-lengths of the cylinder height.

@<Calculate containment range for the cylinder@>=
p->c.y0 = - p->c.height;
p->c.y1 = p->c.height;

@ Let $(x, y, z)$ represent the components of the vector |v| that we
wish to test. Also let $\delta$ represent the magnitude of the
two-dimensional projection of vector |v| on the $xz$-plane, and let
the interval $[y_0, y_1]$ give the containment range of the cylinder
in the $y$-axis. Then the point defined by vector |v| is:

$$\vcenter{\halign{\hfil # & # \hfil \cr
outside the cylinder if & $(y < y_0 \vee y > y_1) \vee (\delta > |radius|)$,\cr
inside the cylinder if & $(y > y_0 \wedge y < y_1) \wedge (\delta < |radius|)$, and\cr
on the surface, & otherwise.\cr
}}$$

Note here that calculation of $\delta$ is delayed until it is
absolutely necessary.

@<Global functions@>=
Containment is_inside_cylinder(const Vector v, const Primitive *p)
{
        double delta;
	if (v[1] < p->c.y0 || v[1] > p->c.y1) return OUTSIDE;
	@<Calculate distance of the two-dimensional $xz$-projection@>;
	if (delta > p->c.radius) return OUTSIDE;
	if (v[1] == p->c.y0 || v[1] == p->c.y1 || delta ==
	p->c.radius)
                return SURFACE;
	return INSIDE;
}

@ Since the origin of the cylinder coincides with the origin of the
world coordinate frame during containment testing, the magnitude of
the two-dimensional projection gives the required distance. 

@<Calculate distance of the two-dimensional $xz$-projection@>=
delta = sqrt(v[0] * v[0] + v[2] * v[2]);

@*2 Containment inside a solid torus.
During containment testing, the origin of the torus coincides with the
origin of the world coordinate frame and is radially symmetrical about
the $y$-axis. Hence, to test containment inside the torus, we first
check if the projected distance $\delta$ of the vector |v| on the
$xz$-plane is outside the radial containment range defined by the
major and minor radii of the torus. If it is, |v| is outside.

@<Information that defines a primitive torus@>=
double r0, r1; /* radial containment range on $xz$-plane */
double phi_end, theta_end; /* end angles in degrees */ 

@ @<Calculate radial containment range for the torus@>=
p->t.r0 = p->t.major - p->t.minor;
p->t.r1 = p->t.major + p->t.minor;

@ Calculate end angle $\phi_1$ after subtending $\phi$ degrees from
the start angle $\phi_0$. Do the same for $\theta$.

Note here that $0^\circ < \phi \le 360^\circ$ and $0^\circ \le \phi_0
< 360^\circ$. Furthermore, we ensure that $0^\circ <
\phi_1 \le 360^\circ$. This is important while testing if an angle
lies {\sl outside} a given range.

@<Calculate end angles for the torus@>=
p->t.phi_end = p->t.phi_start + p->t.phi; 
p->t.theta_end = p->t.theta_start + p->t.theta;
if (p->t.phi_end > 360.0) p->t.phi_end -= 360.0; /* periodicity
adjustment */
if (p->t.theta_end > 360.0) p->t.theta_end -= 360.0;

@ Because of {\sl periodicity}@^periodicity@>, it is easy to check if
an angle $0^{\circ} \le \gamma < 360^{\circ}$ falls outside the region
subtended from $0^{\circ} \le \phi_0 < 360^\circ$ to $0^{\circ} <
\phi_1 \le 360^{\circ}$ by testing one of the following conditions:

$$\vcenter{\halign{\hfil {\rm if}\ $#$, & {\rm $\gamma$ is outside
if}\ $#$ \hfil\cr
\phi_0 < \phi_1 & \gamma < \phi_0 \vee \gamma > \phi_1;\cr
\phi_0 > \phi_1 & \gamma < \phi_0 \wedge \gamma > \phi_1.\cr
}}$$

The second condition only occurs when the value of $\phi_1$ was
adjusted to lie in the range $(0, 360]$. When $\phi_0 =
\phi_1$, we could either have a full rotation, or a no
rotation. However, since $\phi > 0^\circ$, we must always interpret
this as a full rotation; hence, $\gamma$ is always in range when
$\phi_0 = \phi_1$.

@<Global functions@>=
bool angle_outside_range(double gamma, double phi0, double phi1)
{
        if (phi0 == phi1) return false; /* full rotation */
        if (phi0 < phi1) return (gamma < phi0 || gamma > phi1);
        if (phi0 > phi1) return (gamma < phi0 && gamma > phi1); /*
	periodicity adjusted */
        return false;
}

@ If |v| is within the radial containment range, we calculate the
angle $\gamma$ subtended by |v| on the $xz$-plane. Let,
$\gamma^{\circ}$ represent the value of $\gamma$ after converting
radians to degrees. If the torus is partial, i.e., $\phi <
360^{\circ}$, we check if $\gamma^{\circ}$ is outside the range
defined by $\phi_0$ and $\phi_1$. If $\gamma^{\circ}$ is not outside
the range, we check if |v| is outside the volume of the tube. Finally,
if none of the previous conditions were satisfied, we can be certain
that |v| is either inside the torus, or on the surface. To decide
which is valid, we do a final surface test.

Notice that the function |angle_outside_rangle()| implicitly tests if
the torus is partial (i.e., we do not need to test the condition
|p->t.phi < 360.0| explicitly).

@<Global functions@>=
Containment is_inside_torus(const Vector v, const Primitive *p)
{
	double gamma, gamma_deg, tau, tau_deg, delta, radial;
	Vector tube_center, from_tube_center_to_v, temp;
	@<Calculate the projected distance $\delta$ of |v| on the $xz$-plane@>;
	if (delta < p->t.r0 || delta > p->t.r1) return OUTSIDE; /* check radial
	containment on $xz$-plane */
	@<Calculate $\gamma$ subtended by |v| on the $xz$-plane@>;
	if (angle_outside_range(gamma_deg, p->t.phi_start, p->t.phi_end))
                return OUTSIDE;
	@<Check if |v| is outside the tube@>;
	@<Check if |v| is on the surface of the tube@>;
	return INSIDE;
}

@ @<Calculate the projected distance $\delta$ of |v| on the $xz$-plane@>=
delta = sqrt(v[0] * v[0] + v[2] * v[2]);

@ @<Calculate $\gamma$ subtended by |v| on the $xz$-plane@>=
vector_copy(temp, v);
temp[1] = 0.0; /* angle must be on the $xz$-plane */
gamma = vector_angle_radian(positive_xaxis_unit_vector, temp);
gamma_deg = convert_radian_to_degree(gamma);


@ The vector |tube_center| gives the center of the tube which subtends
$\gamma$ radians on the $xz$-plane. Using this point, we calculate
|radial|, which is the radial distance of |v| from |tube_center|. If
|radial| is greater than the minor radius, |v| is outside the
torus. If the tube is partial, we must check if the angle $\tau$
subtended by the vector |from_tube_center_to_v|, from the |tube_center|
to |v|, is outside $\theta$. If it is, |v| is outside the volume.

@<Check if |v| is outside the tube@>=
@<Calculate vector |tube_center| on the center of the tube at $\gamma$@>;
@<Calculate radial distance of |v| from |tube_center|@>;
if (radial > p->t.minor) return OUTSIDE;
@<Calculate the radial angle of |v| on the cross-section at |tube_center|@>;
if (angle_outside_range(tau, p->t.theta_start, p->t.theta_end)) return OUTSIDE;

@ @<Calculate vector |tube_center| on the center of the tube at $\gamma$@>=
tube_center[0] = p->t.major * cos(gamma);
tube_center[1] = 0.0; /* the circular tube center always lies on the $xz$-plane */
tube_center[2] = p->t.major * sin(gamma);

@ @<Calculate radial distance of |v| from |tube_center|@>=
vector_difference(v, tube_center, from_tube_center_to_v);
radial = vector_magnitude(from_tube_center_to_v);

@ @<Calculate the radial angle of |v| on the cross-section at |tube_center|@>=
tau = vector_angle_radian(tube_center, from_tube_center_to_v);
tau_deg = convert_radian_to_degree(tau);

@ @<Check if |v| is on the surface of the tube@>=
if (radial == p->t.minor) return SURFACE;
if (p->t.phi < 360.0 && (gamma == p->t.phi_start || gamma == p->t.phi_end))
        return SURFACE;
if (p->t.theta < 360.0 && (tau == p->t.theta_start || tau == p->t.theta_end))
        return SURFACE;

@*2 Test containment inside a solid primitive.
Containment testing uses different approaches depending on the type
of the primitive solid. Nonetheless, during these tests, all of the
origins of the primitive solids are assumed to coincide with the
origin of the world coordinate frame.

Solid primitives are created internally so that their origin coincides
with the origin of the world coordinate frame. Any translation or
transformation applied henceforth is then stored as operators in the
CSG tree. Thus, while carrying out
|@<Test containment after affine transformations@>|, we apply 
the inverse of the transformations or translations directly to the
point being tested, instead of testing the point against the
transformed solid. This will become clearer in the following
sections.

For now, just remember that all containment tests are carried
out assuming that no transformation or translation has been applied to
the primitive solid. In other words, the origin of the solid will
always coincide with the origin of world coordinate frame, and that
the initial orientation of the solid at creation is preserved. Of
course, these initial orientations are primitive specific, as
discussed in the following sections.

@<Test containment inside primitive solid@>=
affine_inverse(root, v, r);
p = root->leaf.p;
switch(p->type) {
case BLOCK: return is_inside_block(r, p);
case SPHERE: return is_inside_sphere(r, p);
case CYLINDER: return is_inside_cylinder(r, p);
case TORUS: return is_inside_torus(r, p);
default: return INVALID; /* invalid solid */
}

@*2 Test containment inside an intermediate solid.
For primitive solids, test for containment is pretty
straight-forward. We either use distance validation, or parameteric
solutions with respect to the solid. However, for intermediate solids,
which are defined by a CSG tree, we test containment by traversing the
CSG tree. In this section, we use the algorithm described in page 312
of {\sl An Integrated Introduction to Computer Graphics and Geometric
Modelling} by Ronald Goldman [CRC Press (2009)]. The function
|solid_contains_vector(root, v)| returns |true| if the vector |v| is
inside the solid defined by the CSG tree rooted at |root|; otherwise,
|false| is returned. This function assumes that vector |v| is a
homogeneous vector.

@<Global functions@>=
Containment solid_contains_vector(CSG_Node *root, Vector v)
{
	if (NULL == root || NULL == v) return INVALID;
	return recursively_test_containment(root, v);
}

@ Function to recursively test containment.
@<Global functions@>=
Containment recursively_test_containment(CSG_Node *root, Vector v)
{
        Containment left, right;
        Vector r; /* used during inverse transformation */
        if (is_primitive(root)) {
                @<Test containment inside primitive solid@>;
	} else {
                if (is_union(root) || is_intersection(root) || is_difference(root)) {
                        @<Test containment in subtrees using boolean operators@>;
                } else {
                        @<Test containment after affine transformations@>;
                }
        }
	return INVALID;
}

@*3 Containment inside boolean solids.
We test containment inside intermediate solids (i.e., solids defined
as a combination of two solids) using the appropriate boolean
tests. When a point is on the surface of either, or both, the left and
the right solids, we must check separately if the same point will be
inside, or on the surface of the combined solid.

Notice that, when we are carrying out the difference, some points
inside the |left| solid will be on the surface of the new solid if
they were adjacent to points on the surface of the |right| solid after
the subtraction. Since it is difficult to determine this adjacency, we
shall assume that points are on the surface of the new solid only if
they were on the surface of the left solid.

@<Test containment in subtrees using boolean operators@>=
left = recursively_test_containment(root->internal.left, v);
right = recursively_test_containment(root->internal.right, v);
switch(BIT_MASK_NODE & root->op) {
case UNION:
	if (INVALID == left || INVALID == right)
                return INVALID; /* handle error */
	if (INSIDE == left || INSIDE == right)
                return INSIDE;
	if (SURFACE == left || SURFACE == right)
                return SURFACE;
	return OUTSIDE;
case INTERSECTION:
	if (INVALID == left || INVALID == right)
                return INVALID; /* handle error */
	if (OUTSIDE == left || OUTSIDE == right)
                return OUTSIDE;
	if (INSIDE == left && INSIDE == right)
                return INSIDE;
        return SURFACE;
case DIFFERENCE:
	if (INVALID == left || INVALID == right)
                return INVALID; /* handle error */
        if (OUTSIDE == right) {
	        if (SURFACE == left) return SURFACE; /* definitely on
		the surface */
	        if (INSIDE == left) return INSIDE; /* could be on the
		surface, but assume it is inside */
        }
	return OUTSIDE;
default: return INVALID;
}

@*3 Containment after affine transformations.
Let $T_i$ represent an {\sl affine transformation} matrix, which
expresses one of translation, rotation or scaling, and let $T^{-1}_i$
represent its inverse. Furthermore, let $T = \{T_0, \ldots, T_{n-1}\}$
represent an ordered sequence of $n$ transformations. If $s$
represents a solid, and $s'$ represents the solid after applying $T$
to $s$, in ascending order starting with $T_0$, and if $r$ represents
the result of applying to a homogeneous vector $v$ the inverse of the
transformations in $T$ in descending order starting with
$T^{-1}_{n-1}$, i.e. $s' = T_{n-1} \times T_{n-2} \times \ldots \times
T_1 \times T_0 \times s$ and $r = T_0^{-1} \times T_1^{-1} \times
\ldots \times T_{n-2}^{-1} \times T_{n-1}^{-1} \times v$, then
checking if $v$ lies inside $s'$, is equivalent to checking if $r$
lies inside $s$. This is because of the matrix property:

\medskip

\centerline{$(T_{n-1} \times T_{n-2} \times \ldots \times
T_1 \times T_0)^{-1} = T_0^{-1} \times T_1^{-1} \times
\ldots \times T_{n-2}^{-1} \times T_{n-1}^{-1}$}

\medskip

This significantly simplifies the testing of points using the CSG
tree: 1) it is easier to calculate $r$ from $v$, and 2)
containment testing is easier with $s$ than it is with $s'$.

@<Test containment after affine transformations@>=
switch(BIT_MASK_NODE & root->op) {
case TRANSLATE: affine_inverse(root->internal.right, v, r);@+ break;
case ROTATE: affine_inverse(root->internal.right, v, r);@+ break;
case SCALE: affine_inverse(root->internal.right, v, r);@+ break;
default: return INVALID;
}
return recursively_test_containment(root->internal.left, r);

@ The affine matrix $T$ for a translation with displacement vector
$(x, y, z)$ is given by:

\medskip

$T = \left(\matrix{1.0 & 0.0 & 0.0 & x\cr
0.0 & 1.0 & 0.0 & y\cr
0.0 & 0.0 & 1.0 & z\cr
0.0 & 0.0 & 0.0 & 1.0\cr
}\right)$

\medskip

This affine matrix $T$ is stored in the two dimensional array |affine|
using {\sl row-major}@^row-major@> form.

NOTE: In the following initialisations of the affine matrices, we only
set the nonzero elements. All of the node fields have already been
initialised to zero during |create_csg_node()|, since it uses
|mem_typed_alloc()| which in turn uses the |calloc()| system call.

@<Set up the matrix for translation@>=
leaf_node->affine[0][0] = leaf_node->affine[1][1] = leaf_node->affine[2][2] = leaf_node->affine[3][3] = 1.0;
leaf_node->affine[0][3] = op_x; /* $x$-axis translation */
leaf_node->affine[1][3] = op_y; /* $y$-axis translation */
leaf_node->affine[2][3] = op_z; /* $z$-axis translation */

@  The affine matrix $R$ for a rotation of $\theta$ radians about the
axis specified by a unit vector $u = (x, y, z)$ in the world
coordinate frame is given by:

\medskip

$R = \left(\matrix{tx^2 + c & txy - sz & txz + sy & 0.0\cr
txy + sz & ty^2 + c & tyz - sx & 0.0\cr
txz - sy & tyz + sx & tz^2 + c & 0.0\cr
0.0 & 0.0 & 0.0 & 1.0\cr
}\right)$

\medskip

\noindent where, $c = \cos(\theta)$, $s = \sin(\theta)$, and $t = 1 -
\cos(\theta)$. This representation is taken from the article {\sl The
Mathematics of the 3D Rotation Matrix} by Diana Gruber
(http://www.fastgraph.com/makegames/3drotation/).

@<Local variables: |read_geometry()|@>=
double sine, cosine, t, tx, ty, tz, txy, txz, tyz, sx, sy, sz;

@ The matrix $R$ is stored in the two dimensional
array |affine| using {\sl row-major}@^row-major@> form.

@<Set up the matrix for rotation@>=
op_theta *= DEGREE_TO_RADIAN;
sine = -sin(op_theta);
cosine = cos(op_theta);
t = 1.0 - cosine;
tx = t * op_x;
ty = t * op_y;
tz = t * op_z;
txy = tx * op_y;
txz = tx * op_z;
tyz = ty * op_z;
sx = sine * op_x;
sy = sine * op_y;
sz = sine * op_z;
leaf_node->affine[0][0] = tx * op_x + cosine;
leaf_node->affine[0][1] = txy - sz;
leaf_node->affine[0][2] = txz + sy;
leaf_node->affine[1][0] = txy + sz;
leaf_node->affine[1][1] = ty * op_y + cosine;
leaf_node->affine[1][2] = tyz - sx;
leaf_node->affine[2][0] = txz - sy;
leaf_node->affine[2][1] = tyz + sx;
leaf_node->affine[2][2] = tz * op_z + cosine;
leaf_node->affine[3][3] = 1.0;

@ The affine matrix $S$ for a scaling with
scaling factors $(x, y, z)$ is given by:

\medskip

$S = \left(\matrix{x & 0.0 & 0.0 & 0.0\cr
0.0 & y & 0.0 & 0.0\cr
0.0 & 0.0 & z & 0.0\cr
0.0 & 0.0 & 0.0 & 1.0\cr
}\right)$

\medskip

The matrix $S$ is stored in the two dimensional array |affine|
using {\sl row-major}@^row-major@> form.

@<Set up the matrix for scaling@>=
leaf_node->affine[0][0] = op_x; /* $x$-axis scaling */
leaf_node->affine[1][1] = op_y; /* $y$-axis scaling */
leaf_node->affine[2][2] = op_z; /* $z$-axis scaling */
leaf_node->affine[3][3] = 1.0;

@ @<Test geometry input@>=
{
        FILE *f;
        Vector point = ZERO_VECTOR;
	int n;
	Containment flag;
	char c;
	bool t;

        if (false == read_geometry("test/test_geometry_input.data")) exit(1);
	print_geom_statistics(stdout);
	print_sim_world(stdout);
	print_forest();
	if ((f = fopen("test/test_geometry_input_points.data", "r")) == NULL)
               exit(1);
	input_file_current_line = 1;
	while ((c = fgetc(f)) != EOF) {
	       @<Discard comments, white spaces and empty lines@>;
               @<Process containment test-case@>;
        }
error_invalid_command:
	fclose(f);
}

@ @<Process containment test-case@>=
if (c == 'i' || c == 'o' || c == 's') {
        @<Read parameters for the test-case@>;
        @<Validate the test-case@>;
} else {
        fprintf(stderr, "Invalid test command '%c' at line %u\n",
	c, input_file_current_line);
        goto error_invalid_command;
}

@ @<Read parameters for the test-case@>=
n = fscanf(f, "(\"%[^\"]\" %lf %lf %lf)",
        op_solid, &point[0], &point[1], &point[2]);
if (EOF == n || 4 != n) {
        printf("Invalid %s test at line %d\n",@/
        'i' == c ? "inside" :
            ('o' == c ? "outside" : "surface"), input_file_current_line);
        exit(1);
}

@ @<Validate the test-case@>=
temp_node = find_csg_node(op_solid);
if (NULL == temp_node)
   printf("[%d] Solid '%s' not found\n", input_file_current_line, op_solid);
else {
        vector_homogenise(point);
        flag = solid_contains_vector(temp_node, point);
        t = false;
        switch(flag) {
        case INSIDE:@+ if ('i' == c) t = true;@+ break;
        case SURFACE:@+ if ('s' == c) t = true;@+ break;
        case OUTSIDE:@+ if ('o' == c) t = true;@+ break;
        case INVALID: printf("error: ");@+ break;
        }
	printf("[%4d] %s test: %s\n", input_file_current_line,@/
            'i' == c ? "inside" :
                ('o' == c ? "outside" : "surface"),
            t ? "OK" : "Fail");
}


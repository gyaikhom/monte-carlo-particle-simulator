@q This file makes CWEAVE treat Vectors, Events, Particles etc. as reserved words.  Furthermore, we forward-declare types of several data structures. @>

@s Vector int
@s Matrix int
@s Event int
@s Vertex int
@s Particle int

@s int8_t int
@s uint8_t int
@s int16_t int
@s uint16_t int
@s int32_t int
@s uint32_t int
@s GeometryTable int


@ Forward-declare types of several data structures, so that we can
easily reshuffle the sections for easy exposition, while avoiding
compiler error.

@<Forward declare functions@>=
typedef struct geomtab_struct GeometryTable;
void process_and_register_solid(CSG_Node *root);
Containment recursively_test_containment(CSG_Node *root, Vector v);
bool build_subcuboid_trees(GeometryTable *g);
bool build_neighbour_table(GeometryTable *g);
bool build_subcuboids_table(GeometryTable *g);


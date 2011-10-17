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


@ Forward-declare types of several data structures, so that we can
easily reshuffle the sections for easy exposition, while avoiding
compiler error.

@<Forward declare functions@>=
Primitive *create_primitive_solid();
CSG_Node *create_csg_node();
CSG_Node *find_csg_node(char *name);
CSG_Node *register_solid(CSG_Node *solid, char *name);
void destroy_primitive_solid(Primitive *primitive);
void reset_list_of_solids();
void destroy_csg_tree(CSG_Node *temp);
void print_csg_tree(CSG_Node *temp, uint32_t indent);
void process_and_register_solid(const char *name, CSG_Node *root);
Containment recursively_test_containment(CSG_Node *root, Vector v);


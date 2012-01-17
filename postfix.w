@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@*1 Evaluation of CSG boolean expression.

This section is concerned with checking if a three-dimensional point
lies inside a solid, where the solid is represented by a postfix
boolean expression.

NOTE:
The performance of this evaluation can be improved significantly by
first normalising the CSG tree into sum-of-product form before
deriving the boolean expression. Normalisation to sum-of-product form
can be done using the algorithm described by Jack Goldfeather, Steven
Molnar, Greg Turk, and Henry Fuchs in {\sl Near Real-Time CSG
Rendering Using Tree Normalization and Geometric Pruning} [IEEE
Computer Graphics and Applications, pp. 20--28, May {\bf 1989}].

@d MAX_NUM_CSG_PRIMITIVES 1024
@d MAX_NUM_CSG_NODES 2047 /* $2^{\lceil \lg n\rceil + 1} - 1$ */
@<Type definitions@>=
typedef struct {
	int np, nc;
	int c[MAX_NUM_CSG_NODES];
        Primitive *p[MAX_NUM_CSG_PRIMITIVES];
} Solid;

@ We use a boolean stack to evaluate the CSG boolean expression.

@d MAX_BOOLEAN_STACK_SIZE 1024
@<Type definitions@>=
typedef struct {
	int tos, size;
	bool v[MAX_BOOLEAN_STACK_SIZE];
} boolean_stack;

@ @<Global functions@>=
bool boolean_stack_init(boolean_stack *s)
{
	if (NULL == s) return false;
	s->tos = 0;
	s->size = MAX_BOOLEAN_STACK_SIZE;
	return true;
}

@ @<Global functions@>=
bool boolean_stack_push(boolean_stack *s, bool v)
{
	if (s->tos == s->size) return false;
	s->v[s->tos++] = v;
	return true;
}

@ @<Global functions@>=
bool boolean_stack_pop(boolean_stack *s, bool *v)
{
	if (0 == s->tos) return false;
	*v = s->v[--s->tos];
	return true;
}

@ @<Global functions@>=
bool is_inside_primitive(Vector v, Primitive *p)
{
	Containment c;
	switch(p->type) {
	case BLOCK: c = is_inside_block(v, p); break;
	case SPHERE: c = is_inside_sphere(v, p); break;
	case CYLINDER: c = is_inside_cylinder(v, p); break;
	case TORUS: c = is_inside_torus(v, p); break;
	default: c = INVALID; /* invalid solid */
	}
	if (INSIDE == c || SURFACE == c) return true;
	return false;
}

@ @(mcs.cu@>=
__device__ bool cuda_is_inside_primitive(Vector v, Primitive *p)
{
	Containment c;
	switch(p->type) {
	case BLOCK: c = cuda_is_inside_block(v, p); break;
	case SPHERE: c = cuda_is_inside_sphere(v, p); break;
	case CYLINDER: c = cuda_is_inside_cylinder(v, p); break;
	case TORUS: c = cuda_is_inside_torus(v, p); break;
	default: c = INVALID; /* invalid solid */
	}
	if (INSIDE == c || SURFACE == c) return true;
	return false;
}

@
@d BOOLEAN_DIFFERENCE -1
@d BOOLEAN_INTERSECTION -2
@d BOOLEAN_UNION -3
@<Global functions@>=
bool is_inside(Vector v, Solid *s, bool *result)
{
	boolean_stack stack;
	bool l, r; /* left and right operands */
	bool cache[MAX_NUM_CSG_PRIMITIVES]; /* boolean cache */
	int i;
	@<Initialise stack and boolean cache@>;
	@<Evaluate the postfix boolean expression using the stack@>;
	boolean_stack_pop(&stack, result); /* must be the only item in
	stack */
	@<Check if the CSG expression was valid@>;
	return true;

	@<Handle irrecoverable error: |is_inside(v, s, result)|@>;
	return false;
}

@ The same primitive could appear more than once inside the CSG
expression. Since, checking whether a point lies inside a primitive is
expensive, we do such evaluations only once and cache the value
for use by subsequent appearances of the primitive. Since, all of the
primitives will appear atleast once, we evaluate them beforehand so
that the cache values are either |true| or |false| for each primitive.

@<Initialise stack and boolean cache@>=
boolean_stack_init(&stack);
for (i = 0; i < s->np; ++i) cache[i] = is_inside_primitive(v, s->p[i]);

@ The CSG expression |c| is an array of integers, where all of the
positive values gives an index inside the primitives array, and any
negative value must be either -1, -2, or -3, denoting respectively a
boolean difference, intersection, or union. Any other value in the
expression is invalid. For instance, the integer sequence $c = \{0, 3,
-2, 1, -1, 2, 3, -2, -3\}$ is a postfix representation for the boolean
expression $((A \cap D) - B) \cup (C \cap D)$, where the primitives
array $p = \{A, B, C, D\}$.

@<Evaluate the postfix boolean expression using the stack@>=
for (i = 0; i < s->nc; ++i) {
        if (BOOLEAN_DIFFERENCE < s->c[i]) {
	        @<Push to stack the boolean containment of $v$ inside primitive@>;
	} else {
	        if (BOOLEAN_UNION > s->c[i]) {
		        fprintf(stderr,
			      "Invalid value '%d' in expression\n",
                              s->c[i]);
			return false;
		}
	        @<Evaluate boolean operator and push result into stack@>;
	}
}

@ @<Push to stack the boolean containment of $v$ inside primitive@>=
if (!boolean_stack_push(&stack, cache[s->c[i]]))
        goto stack_full;

@ The first stack pop gives the right operand |r|, and the second
gives the left operand |l|. Since, boolean difference is
noncommutative, we must preserve the order during its evaluation
(i.e., the negation).

@<Evaluate boolean operator and push result into stack@>=
if (!boolean_stack_pop(&stack, &r)) goto stack_empty;
if (!boolean_stack_pop(&stack, &l)) goto stack_empty;
switch(s->c[i]) {
case BOOLEAN_DIFFERENCE:
        if (!boolean_stack_push(&stack, l && !r)) goto stack_full;
	break;
case BOOLEAN_INTERSECTION:
        if (!boolean_stack_push(&stack, l && r)) goto stack_full;
        break;
case BOOLEAN_UNION:
	if (!boolean_stack_push(&stack, l || r)) goto stack_full;
        break;
default:;
}

@ A CSG boolean expression is valid if the stack is empty after the
result has been popped out. Hence, the following pop for |l| must fail
for valid CSG expressions since it is executed after porring the result.

@<Check if the CSG expression was valid@>=
if (boolean_stack_pop(&stack, &l)) {
        fprintf(stderr, "Invalid CSG tree expression\n");
	return false;
}

@ @<Handle irrecoverable error: |is_inside(v, s, result)|@>=
stack_empty:
	fprintf(stderr, "Boolean stack is empty... ");
	goto exit_error;

stack_full:
	fprintf(stderr, "Boolean stack is full... ");
	goto exit_error;

exit_error:
	fprintf(stderr, "while evaluating '%d' at index %d\n",
	s->c[i], i);

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
if (NULL == temp_node) {
        printf("[%d] Solid '%s' not found\n",
	input_file_current_line, op_solid);
} else {
        vector_homogenise(point);
        flag = solid_contains_vector(temp_node, point);
        t = false;
        switch(flag) {
        case INSIDE:
                if ('i' == c) t = true;
                break;
        case SURFACE:
                if ('s' == c) t = true;
                break;
        case OUTSIDE:
                if ('o' == c) t = true;
                break;
        case INVALID:
                printf("error: ");
                break;
        }
	printf("[%4d] %s test: %s\n", input_file_current_line,@/
            'i' == c ? "inside" :
                ('o' == c ? "outside" : "surface"),
            t ? "OK" : "Fail");
}

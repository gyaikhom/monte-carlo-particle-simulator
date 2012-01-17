@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@*1 Evaluation of CSG boolean expression.

This section is concerned with checking if a three-dimensional point
lies inside a solid, where the solid is represented by a postfix
boolean expression.

@ Function |is_inside_primitive(v,p)| checks if the vector |v| lies
inside the primtive |p|. It returns |true| if |v| is inside |p|; or
|false|, otherwise. The value of |v| and |p| are both left
unmodified.

@<Global functions@>=
bool is_inside_primitive(const Vector v, const Primitive *p)
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

@ Function |is_inside(v,s,r)| checks if the vector |v| is inside the
solid which may be retrieved from the {\it solid indices buffer} using
the index |s|. To carry out the containment test, the boolean postfix
expression that corresponds to the CSG tree of the solid is evaluated
against the vector |v|. If this evaluation was successful, the
function returns |true|, and the actual containment result is stored
in |r|. If there was an error, however, e.g., the stack was full, the
function will return a |false|. In this case, value in |r| must be
discarded.

NOTE:
The performance of this evaluator can be improved significantly by
first normalising the CSG tree into sum-of-product form. We have not
pursue this optimisation here due to time constraints. Normalisation
to sum-of-product form can be done using the algorithm described by
Jack Goldfeather, Steven Molnar, Greg Turk, and Henry Fuchs in {\sl
Near Real-Time CSG Rendering Using Tree Normalization and Geometric
Pruning} [IEEE Computer Graphics and Applications, pp. 20--28, May
{\bf 1989}].

@d BOOLEAN_DIFFERENCE -1
@d BOOLEAN_INTERSECTION -2
@d BOOLEAN_UNION -3
@d BOOLEAN_INVALID -4
@<Global functions@>=
bool is_inside(Vector v, uint32_t solid, bool *result)
{
	uint32_t i, c; /* start index and length of the postfix expression */
	int item = BOOLEAN_INVALID; /* current postfix item */
	boolean_stack stack;
	boolean_stack_init(&stack);
	@<Retrieve start index and length of the postfix expression@>;
	@<Evaluate the boolean postfix expression using the stack@>;
	@<Check if the CSG expression was valid and retrieve the result@>;
	return true;

	@<Handle irrecoverable error: |is_inside(v, s, result)|@>;
	return false;
}

@ We use the {\it solids table} component of the geometry table to
retrieve the start index and length of the boolean postfix
expression representing the solid. These values will then be used to
retrieve the boolean expression from the postfix expression buffer.
 
@<Retrieve start index and length of the postfix expression@>=
i = geotab.s[solid].s;
c = geotab.s[solid].c;

@ The CSG expression is an array of integers, where all of the
positive values give an index inside the {\it primitives table}
component of the geomatry table. Any negative value must be either -1,
-2, or -3, denoting respectively a boolean difference, intersection,
or union. Any other value in the expression is invalid. For instance,
the integer sequence $c = \{0, 3, -2, 1, -1, 2, 3, -2, -3\}$ is a
postfix representation for the boolean expression $((A \cap D) - B)
\cup (C \cap D)$, where the primitives array $p = \{A, B, C, D\}$.

@<Evaluate the boolean postfix expression using the stack@>=
while (c--) {
      item = geotab.pb[i++]; /* lookup current item, and move to next item */
      if (BOOLEAN_DIFFERENCE < item) {
        @<Push to stack the boolean containment of $v$ inside primitive@>;
      } else {
	@<Evaluate boolean operator and push result into stack@>;
      }
}

@ We lookup the {\it primitives table} component of the
geometry table to retrieve the primitive |p| that corresponds to the
current postfix item. We then call |is_inside_primitive(v,p)| to check if
the vector |v| is inside the primitive |p|. The returned value is then
pushed to the stack, to be retrieved later as an operand to a subsequent
boolean operator.
@<Push to stack the boolean containment of $v$ inside primitive@>=
bool t = is_inside_primitive(v, &geotab.p[item].p);
if (!boolean_stack_push(&stack, t)) goto stack_full;

@ The first stack pop gives the right operand |r|, and the second
gives the left operand |l|. Since, boolean difference is
noncommutative, we must preserve the order during its evaluation.

@<Evaluate boolean operator and push result into stack@>=
bool l, r; /* left and right operands */
if (!boolean_stack_pop(&stack, &r)) goto stack_empty;
if (!boolean_stack_pop(&stack, &l)) goto stack_empty;
switch (item) {
case BOOLEAN_DIFFERENCE:
        if (!boolean_stack_push(&stack, l && !r)) goto stack_full;
	break;
case BOOLEAN_INTERSECTION:
        if (!boolean_stack_push(&stack, l && r)) goto stack_full;
        break;
case BOOLEAN_UNION:
	if (!boolean_stack_push(&stack, l || r)) goto stack_full;
        break;
default:
	fprintf(stderr, "Invalid value '%d' in expression\n", item);
	return false;
}

@ A CSG boolean expression is valid if the stack is empty after the
result has been popped out. Hence, the following test pop for |t| must fail
for valid CSG expressions.

@<Check if the CSG expression was valid and retrieve the result@>=
{
	bool t;
	boolean_stack_pop(&stack, result); /* must be the only item in stack */
	if (boolean_stack_pop(&stack, &t)) {
           fprintf(stderr, "Invalid CSG tree expression\n");
	   return false;
	}
}

@ @<Handle irrecoverable error: |is_inside(v, s, result)|@>=
stack_empty:
	fprintf(stderr, "Boolean stack is empty... ");
	goto exit_error;

stack_full:
	fprintf(stderr, "Boolean stack is full... ");
	goto exit_error;

exit_error:
	fprintf(stderr, "while evaluating '%d' at index %d\n", item, i);

@ Function |solid_contains_particle(s,p)| checks if the solid |s|
contains the particle |p|. If it does, |true| is returned; 
otherwise, |false| is returned. To carry out the actual check, this
function tests the containment of the particle's position vector
|p->v| against the CSG expression that defines the solid.

@<Global functions@>=
bool solid_contains_particle(uint32_t s, Particle *p)
{
      bool result;
      if (is_inside(p->v, s, &result)) return result;
      else return false;
}

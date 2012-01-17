@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@*1 Boolean stack.

We use a boolean stack to evaluate the CSG boolean expression.

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


@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Particle repository. All of the particles yet to be simulated are
stored in a particles repository. The particles are grouped according
to their containing subcuboid, so that each group with the same solids
can be simulated in one batch.

@<Type definitions@>=
typedef struct particle_struct Particle;

@ The particle repository is implemented as a paged binary
max-heap, which is allowed to grow. The heap is maintained as 
a linked list of heap pages, where each heap page is a fixed array
of particles.

@d HEAP_PAGE_SIZE 41 /* 4096 bytes per page, including next pointer */
@<Type definitions@>=
typedef struct particle_repository_struct {
	uint32_t count; /* current number of particles in heap */
	uint32_t max; /* maximum number of particles allowed in heap */
	uint16_t page_size; /* number of particles per page */
	Particle *head, *tail; /* pointers to first and last pages */
} ParticleRepository;
ParticleRepository particles = {0, 0, HEAP_PAGE_SIZE, NULL, NULL};

@ @<Type definitions@>=
enum {
	HEAP_SUCCESS = 0,
	HEAP_ERROR_UNDEFINED = -1,
	HEAP_ERROR_ALLOC = -2,
	HEAP_EMPTY = 1
};

@ The paging only becomes active when the first allocated heap page is
insufficient to fulfill the required number of particles. Since paging
adds computational overhead, it is expected that the first page is
allocated carefully so as to avoid paging. Keeping this in mind, we
provide two versions of all the heap functions. The functions with
the `fast' suffix does not assume paging, whereas, those with `paged'
suffix assume paging.
@<Global functions@>=
void heap_insert_fast(ParticleRepository *r, const Particle *p)
{
	uint32_t parent, node;
	Particle *page = r->head;
	@<In the first page, bubble up the new node to its rightful place@>;	
	@<In the first page, place node in its rightful place@>;
	r->count++;
}
@ @<In the first page, bubble up the new node to its rightful place@>=
node = r->count + 1;
parent = node >> 1;
while (parent > 0) {
	if (page[parent].subcuboid < p->subcuboid) page[node] = page[parent];
        else break;
	node = parent;
	parent = node >> 1;
}

@ @<In the first page, place node in its rightful place@>=
page[node] = *p;

@ @<Global functions@>=
void heap_remove_fast(ParticleRepository *r, Particle *p)
{
	uint32_t node, left, right;
	Particle temp, *page;
	page = r->head;
	*p = page[1];
	if (r->count > 1) {
	   temp = page[r->count];
	   @<In the first page, fix the heap by bubbling down the last node@>;
	}
	r->count--;
}

@ @<In the first page, fix the heap by bubbling down the last node@>=
node = 1;
left = 2;
while (left < r->count) {
      right = left + 1;		
      if (page[r->count].subcuboid < page[left].subcuboid) {
      	 if (right < r->count &&
	     page[left].subcuboid < page[right].subcuboid) {
	     page[node] = page[right];
	     node = right;
	 } else {
	     page[node] = page[left];
	     node = left;
	 }
      } else {
	  if (right < r->count &&
	     page[r->count].subcuboid < page[right].subcuboid) {
	     page[node] = page[right];
	     node = right;
	  } else break;
      }
      left = node << 1; /* go deeper until we have passed a leaf */
}
page[node] = temp; /* place last node in rightful place */

@ For a given particle index |n| within the heap |r|, function
|heap_find_pidx(r, n,p,i)| finds the page start address |p| and the
index |i| within that page.

The binary heap is implemented using an array representation of a
complete binary tree. Hence, while maintaining this tree we are
required to find the parent of a given tree node. If paging is active,
we are sometimes required to find the page and index within the page
for a given particle.
@<Global functions@>=
void heap_find_pidx(ParticleRepository *r, uint32_t n,
     Particle **p, uint32_t *idx)
{
	int i;
	Particle *t = r->head;
	n--;
	for (i = n / r->page_size; i; --i) t = (Particle *) *(char **)t;
	*idx = (n % r->page_size) + 1;
	*p = t;
}

@ @<Global functions@>=
void heap_insert_paged(ParticleRepository *r, const Particle *p)
{
	uint32_t parent, node;
	uint32_t parent_idx, node_idx;
	Particle *parent_page, *node_page;
	@<Bubble up the new node until we find its rightful place@>;	
	r->count++;
}

@ @<Bubble up the new node until we find its rightful place@>=
node = r->count + 1; /* start bubbling at the last node */
parent = node >> 1;
while (parent > 0) {
	heap_find_pidx(r, parent, &parent_page, &parent_idx);
        heap_find_pidx(r, node, &node_page, &node_idx);
	if (parent_page[parent_idx].subcuboid < p->subcuboid)
	    node_page[node_idx] = parent_page[parent_idx];
	else {
	    node_page[node_idx] = *p;
	    break; /* rightful place found */
	}
	node = parent;
	parent = node >> 1; /* bubble up */
}
if (parent == 0) {
    heap_find_pidx(r, node, &node_page, &node_idx);
    node_page[node_idx] = *p;
}

@ @<Global functions@>=
void heap_remove_paged(ParticleRepository *r, Particle *p)
{
	uint32_t last_idx, left_idx, right_idx, node_idx;
	uint32_t node, left, right;
	Particle temp;
	Particle *last_page, *left_page, *right_page, *node_page;
	*p = r->head[1]; /* particle currently at the top of the
	heap, which will be removed */
	if (r->count > 1) {
	@<Get the page and index for the last node@>;
	@<Fix the heap by bubbling down the last node@>;
	@<Finalise bubbling down by placing the last node in its rightful place@>;
	}
	r->count--;
}

@ @<Get the page and index for the last node@>=
heap_find_pidx(r, r->count, &last_page, &last_idx);
temp = last_page[last_idx];

@ @<Fix the heap by bubbling down the last node@>=
node = 1;
left = 2;
while (left < r->count) {
	right = left + 1;		
	@<Get page and index for left, right and current node@>;
	@<Bubble down the last node relative to the current node@>;
	left = node << 1; /* go deeper until we have passed a leaf */
}

@ @<Get page and index for left, right and current node@>=
heap_find_pidx(r, left, &left_page, &left_idx);
heap_find_pidx(r, right, &right_page, &right_idx);
heap_find_pidx(r, node, &node_page, &node_idx);

@ @<Bubble down the last node relative to the current node@>=
if (last_page[last_idx].subcuboid
    < left_page[left_idx].subcuboid) {
	if (left_page[left_idx].subcuboid
	    < right_page[right_idx].subcuboid) {
	    node_page[node_idx] = right_page[right_idx];
	    node = right;
	} else {
	    node_page[node_idx] = left_page[left_idx];
	    node = left;
	}
} else {
    if (last_page[last_idx].subcuboid
        < right_page[right_idx].subcuboid) {
	node_page[node_idx] = right_page[right_idx];
	node = right;
    } else break;
}

@ @<Finalise bubbling down by placing the last node in its rightful place@>=
heap_find_pidx(r, node, &node_page, &node_idx); /* find rightful place */
node_page[node_idx] = temp; /* place last node in rightful place */

@ @<Global functions@>=
void heap_print(FILE *f, ParticleRepository *r)
{
	int i, j, p;
	Particle *page;
	fprintf(f, "Number of particles: %u", r->count);
	p = 0;
	j = r->count;
	page = r->head;
	while (page) {
		if (j) {
			fprintf(f, "\nPage %u: ", p);
			for (i = 1; j && i <= r->page_size; ++i, --j)
				fprintf(f, "%u ", page[i].subcuboid);
			page = (Particle *) *(char **) page;
			p++;
		} else break;
	}
	fprintf(f, "\n");
}

@ When there is no space available to insert a new particle, the heap
must be expanded. This involves creating a new page and the associated
particles array.
@<Global functions@>=
int heap_expand(ParticleRepository *r)
{
	Particle* t = mem_typed_alloc(HEAP_PAGE_SIZE + 1, Particle, mem_p);
	if (NULL == t) return HEAP_ERROR_ALLOC;
	*(char **) t = NULL;
	r->max += r->page_size;
	*((char **) r->tail) = (char *) t;
	r->tail = t;
	printf("expanding...\n");
	heap_print(stdout, r);
	fprintf(stdout, "\n");
	return HEAP_SUCCESS;
}

@ @<Global functions@>=
int heap_insert(ParticleRepository *r, Particle *p)
{
	int i;
	if (r->count >= r->max) {
	   i = heap_expand(r);
	   if (i) return i;
	   heap_insert_paged(r, p);
	} else {
	   if (r->count > r->page_size) heap_insert_paged(r, p); 
	   else heap_insert_fast(r, p);
	}
	return HEAP_SUCCESS;
}

@ @<Global functions@>=
int heap_remove(ParticleRepository *r, Particle *p)
{
	if (r->count == 0) return HEAP_EMPTY;
	if (r->count > r->page_size) heap_remove_paged(r, p);
	else heap_remove_fast(r, p);
	return HEAP_SUCCESS;
}

@ @<Global functions@>=
int heap_init(ParticleRepository *r)
{
	r->head = mem_typed_alloc(HEAP_PAGE_SIZE + 1, Particle, mem_p);
	if (NULL == r->head) return HEAP_ERROR_ALLOC;
	r->count = 0;
	r->max = r->page_size = HEAP_PAGE_SIZE;
	*(char **) r->head = NULL;
	r->tail = r->head;
	return HEAP_SUCCESS;
}

@ @<Test particle repository@>=
{
	Particle p;
	heap_init(&particles);
	do {
	    printf("Enter subcuboid (enter 0 to remove from heap): ");
	    scanf("%u", &(p.subcuboid));
	    if (p.subcuboid == 0) {
	        if (heap_remove(&particles, &p) != HEAP_SUCCESS) break;
		printf("%u\n", p.subcuboid);
	    } else {
	        if (heap_insert(&particles, &p) != HEAP_SUCCESS) break;
	    }
	    heap_print(stdout, &particles);
	} while (1);
}

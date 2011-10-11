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

@<Type definitions@>=
typedef struct heap_page_struct {
	Particle *start; /* start of the particle array */
	struct heap_page_struct *next; /* pointer to next page */
} HeapPage;

@
@d heap_nodes_per_block 16 /* each memory block will have 1024 bytes */
@<Global functions@>=
HeapPage *bad_heap_page = NULL, *next_heap_page = NULL;
HeapPage* create_heap_page()
{
	HeapPage *t = next_heap_page;
	if (t == bad_heap_page) {
	   t = mem_typed_alloc(heap_nodes_per_block, HeapPage, mem_p);
	   if (t == NULL) return NULL;
	   else {
	   	next_heap_page = t + 1;
		bad_heap_page = t + heap_nodes_per_block;
	   }
	} else next_heap_page++;
        return t;
}

@
@d HEAP_PAGE_SIZE 4 /* 16384 bytes per page */
@<Type definitions@>=
typedef struct particle_repository_struct {
	uint32_t count; /* current number of particles in heap */
	uint32_t max; /* maximum number of particles allowed in heap */
	uint16_t page_size; /* number of particles per page */
	HeapPage *head, *tail; /* pointers to first and last pages */
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
	Particle *page = r->head->start;
	@<In the first page, bubble up the new node to its rightful place@>;	
	@<In the first page, place node in its rightful place@>;
	r->count++;
}
@ @<In the first page, bubble up the new node to its rightful place@>=
node = r->count + 1;
parent = node >> 1;
while (parent > 0) {
	if (page[parent - 1].subcuboid < p->subcuboid)
	    page[node - 1] = page[parent - 1];
        else break;
	node = parent;
	parent = node >> 1;
}

@ @<In the first page, place node in its rightful place@>=
page[node - 1] = *p;

@ @<Global functions@>=
void heap_remove_fast(ParticleRepository *r, Particle *p)
{
	uint32_t node, left, right;
	Particle temp, *page;
	page = r->head->start;
	*p = page[0];
	if (r->count > 1) {
	   temp = page[r->count - 1];
	   @<In the first page, fix the heap by bubbling down the last node@>;
	}
	r->count--;
}

@ @<In the first page, fix the heap by bubbling down the last node@>=
node = 1;
left = 2;
while (left < r->count) {
      right = left + 1;		
      if (page[r->count - 1].subcuboid < page[left - 1].subcuboid) {
      	 if (right < r->count &&
	     page[left - 1].subcuboid < page[right - 1].subcuboid) {
	     page[node - 1] = page[right - 1];
	     node = right;
	 } else {
	     page[node - 1] = page[left - 1];
	     node = left;
	 }
      } else {
	  if (right < r->count &&
	     page[r->count - 1].subcuboid < page[right - 1].subcuboid) {
	     page[node - 1] = page[right - 1];
	     node = right;
	  } else break;
      }
      left = node << 1; /* go deeper until we have passed a leaf */
}
page[node - 1] = temp; /* place last node in rightful place */

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
     HeapPage **p, uint32_t *idx)
{
	int i;
	*p = r->head;
	for (i = n / r->page_size; i; --i) *p = (*p)->next;
	*idx = n % r->page_size;
}

@ @<Global functions@>=
void heap_insert_paged(ParticleRepository *r, const Particle *p)
{
	uint32_t parent, node;
	uint32_t parent_idx, node_idx;
	HeapPage *parent_page, *node_page;
	@<Bubble up the new node until we find its rightful place@>;	
	@<Place node in its rightful place@>;
	r->count++;
}

@ @<Bubble up the new node until we find its rightful place@>=
node = r->count + 1; /* start bubbling at the last node */
parent = node >> 1;
while (parent > 0) {
	heap_find_pidx(r, parent - 1, &parent_page, &parent_idx);
	if (parent_page->start[parent_idx].subcuboid
	    < p->subcuboid) {
		heap_find_pidx(r, node - 1, &node_page, &node_idx);
		node_page->start[node_idx]
			= parent_page->start[parent_idx];
	} else break; /* rightful place found */
	node = parent;
	parent = node >> 1; /* bubble up */
}

@ @<Place node in its rightful place@>=
heap_find_pidx(r, node - 1, &node_page, &node_idx);
node_page->start[node_idx] = *p;

@ @<Global functions@>=
void heap_remove_paged(ParticleRepository *r, Particle *p)
{
	uint32_t last_idx, left_idx, right_idx, node_idx;
	uint32_t node, left, right;
	Particle temp;
	HeapPage *last_page, *left_page, *right_page, *node_page;
	*p = r->head->start[0]; /* particle currently at the top of the
	heap, which will be removed */
	if (r->count > 1) {
	@<Get the page and index for the last node@>;
	@<Fix the heap by bubbling down the last node@>;
	@<Finalise bubbling down by placing the last node in its rightful place@>;
	}
	r->count--;
}

@ @<Get the page and index for the last node@>=
heap_find_pidx(r, r->count - 1, &last_page, &last_idx);
temp = last_page->start[last_idx];

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
heap_find_pidx(r, left - 1, &left_page, &left_idx);
heap_find_pidx(r, right - 1, &right_page, &right_idx);
heap_find_pidx(r, node - 1, &node_page, &node_idx);

@ @<Bubble down the last node relative to the current node@>=
if (last_page->start[last_idx].subcuboid
    < left_page->start[left_idx].subcuboid) {
	if (left_page->start[left_idx].subcuboid
	    < right_page->start[right_idx].subcuboid) {
	    node_page->start[node_idx] = right_page->start[right_idx];
	    node = right;
	} else {
	    node_page->start[node_idx] = left_page->start[left_idx];
	    node = left;
	}
} else {
    if (last_page->start[last_idx].subcuboid
        < right_page->start[right_idx].subcuboid) {
	node_page->start[node_idx] = right_page->start[right_idx];
	node = right;
    } else break;
}

@ @<Finalise bubbling down by placing the last node in its rightful place@>=
heap_find_pidx(r, node - 1, &node_page, &node_idx); /* find rightful place */
node_page->start[node_idx] = temp; /* place last node in rightful
place */

@ When there is no space available to insert a new particle, the heap
must be expanded. This involves creating a new page and the associated
particles array.
@<Global functions@>=
int heap_expand(ParticleRepository *r)
{
	HeapPage* t = create_heap_page();
	if (NULL == t) return HEAP_ERROR_ALLOC;
	t->start = mem_typed_alloc(r->page_size, Particle, mem_p);
	if (NULL == t->start) return HEAP_ERROR_ALLOC;
	t->next = NULL;
	r->max += r->page_size;
	r->tail->next = t;
	r->tail = t;
	return HEAP_SUCCESS;
}

@ @<Global functions@>=
int heap_insert(ParticleRepository *r, Particle *p)
{
	int i;
	if (r->count < r->max) {
	   if (r->count < r->page_size) heap_insert_fast(r, p);
	   else heap_insert_paged(r, p);
	} else {
	   i = heap_expand(r);
	   if (i) return i;
	   heap_insert_paged(r, p);
	}
	return HEAP_SUCCESS;
}

@ @<Global functions@>=
int heap_remove(ParticleRepository *r, Particle *p)
{
	if (r->count == 0) return HEAP_EMPTY;
	if (r->count < r->page_size) heap_remove_fast(r, p);
	else heap_remove_paged(r, p);
	return HEAP_SUCCESS;
}

@ @<Global functions@>=
int heap_init(ParticleRepository *r)
{
	r->count = 0;
	r->max = r->page_size = HEAP_PAGE_SIZE;
	r->head = create_heap_page();
	if (NULL == r->head) return HEAP_ERROR_ALLOC;
	r->head->start = mem_typed_alloc(r->page_size, Particle, mem_p);
	if (NULL == r->head->start) return HEAP_ERROR_ALLOC;
	r->head->next = NULL;
	r->tail = r->head;
	return HEAP_SUCCESS;
}

@ @<Global functions@>=
void heap_print(FILE *f, ParticleRepository *r)
{
	int i, j, p;
	HeapPage *page;
	fprintf(f, "Number of particles: %u", r->count);
	p = 0;
	j = r->count;
	page = r->head;
	while (page) {
		if (j) {
			fprintf(f, "\nPage %u: ", p);
			for (i = 0; j && i < r->page_size; ++i, --j)
				fprintf(f, "%u ", page->start[i].subcuboid);
			page = page->next;
			p++;
		} else break;
	}
	fprintf(f, "\n");
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

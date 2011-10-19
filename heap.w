@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Particle repository. All of the particles yet to be simulated are
stored in a particles repository. The particles are grouped according
to their containing subcuboid, so that each group with the same solids
can be simulated in one batch.

The particle repository is implemented as a paged binary
{\sl max-heap} (we could have use a {\sl min-heap}, but this doesn't
matter), which is allowed to grow. The heap is maintained as a linked
list of heap pages, where each heap page is a fixed array of
particles, with the first array element used as a
pointer-to-a-pointer to maintain the linked list of heap pages.

\bigskip

\centerline{\epsfig{file=figures/particle-repository.mps,scale=1}}

\smallskip

The paging only becomes active when the first allocated heap page is
insufficient to fulfill the required number of particles. Since paging
adds computational overhead, it is expected that the first page is
allocated carefully so as to avoid paging. Keeping this in mind, we
provide two versions of all the heap functions. The functions with
the `fast' suffix does not assume paging, whereas, those with `paged'
suffix assume paging.

@d HEAP_PAGE_SIZE 4 /* 4096 bytes per page, including linked list `next' pointer */
@<Type definitions@>=
typedef struct particle_struct Particle;
typedef struct particle_repository_struct {
	uint32_t pid; /* particle identifier to use next */
	uint32_t count; /* current number of particles in heap */
	uint32_t max; /* maximum number of particles allowed in heap */
	uint16_t page_size; /* number of particles per page */
	Particle *head, *tail; /* pointers to first and last pages */
} ParticleRepository;
ParticleRepository particles = {1, 0, 0, HEAP_PAGE_SIZE, NULL, NULL};

@ Function |heap_insert_fast(pr,t)| inserts a new particle |t| into the
particles repository |pr|. This insertion does not assume heap paging.
In the first page, we {\sl bubble up} the new node to its rightful
place using the containing subcuboid as the comparison key. This
bubbling is done by climbing up the tree starting at the array
position dictated by the array implementation of a complete binary
tree, and moving the node from parent to child as we climb. The
insertion has completed when the node is finally placed at its rightful
place.
@<Global functions@>=
void heap_insert_fast(ParticleRepository *pr, Particle *t)
{
	uint32_t p, n; /* indices of parent and current node */
	Particle *fp;
	fp = pr->head; /* first page */
	n = ++pr->count; /* start bubbling from the last node */
	p = n >> 1;
	while (p > 0) { /* bubble up the tree */
	      if (fp[p].subcuboid < t->subcuboid)
	          fp[n] = fp[p];  /* move from parent to child */
              else break; /* rightful place found */
	      n = p; /* climb up the tree */
	      p = n >> 1;
	}
	if (0 == t->id) t->id = pr->pid++; /* give unique id to particle */
	fp[n] = *t; /* place node */
}

@ Function |heap_remove_fast(pr,t)| removes the particle at the top of
the max-heap representing the particles repository |pr|, and stores the
particle into |t|. This removal does not assume heap paging. From the
first page, we remove the first node, which is always the top of the
heap. We then fill this void by promoting the last node of the
complete binary tree to become the root. Finally, we rebalance the
heap by {\sl bubbling down} the root until we find its rightful place,
or we find that it has becomes a leaf. While bubbling down, we move
either the left or right subtree root to the current node as we climb
down the tree.

@<Global functions@>=
void heap_remove_fast(ParticleRepository *pr, Particle *t)
{
	uint32_t n, l, r, q; /* indices to node, left, right and last */
	Particle temp, *fp;
	fp = pr->head; /* first page */
	*t = fp[1]; /* particle to return */
	q = pr->count--;
	if (q) {
	   temp = fp[q]; /* the last node to be promoted */
	   n = 1; /* start at the root */
	   l = 2;
	   while (l < q) { /* bubble down the tree choosing left or
	   right subtree */
      	      r = l + 1;		
      	      if (fp[q].subcuboid < fp[l].subcuboid) {
      	          if (r < q && fp[l].subcuboid < fp[r].subcuboid) {
	     	      fp[n] = fp[r];
	     	      n = r;
	 	  } else {
	     	      fp[n] = fp[l];
	     	      n = l;
	 	  }
      	      } else {
	          if (r < q && fp[q].subcuboid < fp[r].subcuboid) {
	     	      fp[n] = fp[r];
	     	      n = r;
	  	  } else break; /* rightful place found */
      	      }
      	      l = n << 1; /* go deeper until we have passed a leaf */
	  }
          fp[n] = temp; /* place node */
	}
}


@ For a given particle index |n| within the heap |pr|, function
|heap_find_pidx(pr,n,p,i)| finds the page start address |p| and the
index |i| within that page. The binary heap is implemented using an
array representation of a complete binary tree. Hence, while
maintaining this tree we are required to find the parent, or children,
of a given tree node, and if paging is active, we are required to find
its page and index within the heap.
@<Global functions@>=
void heap_find_pidx(ParticleRepository *pr, uint32_t n,
     Particle **p, uint32_t *i)
{
	Particle *t = pr->head;
	while (n > pr->page_size) {
	    t = (Particle *) *(char **)t;
	    n -= pr->page_size;
	}
	*i = n;
	*p = t;
}

@ Function |heap_insert_paged(pr,t)| inserts a new particle |t| into the
particles repository |pr|. This insertion uses the same method as
|heap_insert_fast(pr,t)|, except that it assumes heap paging; hence,
to reference, place, or move a node, we first find the correct page
and index for each node using |heap_find_pidx(pr,n,p,i)|.

@<Global functions@>=
void heap_insert_paged(ParticleRepository *pr, Particle *t)
{
	uint32_t p, n; /* indices of parent and current node within
	heap */
	Particle *p_pg, *n_pg; /* heap pages of parent and current node */
	uint32_t p_idx, n_idx; /* indices within heap pages of parent and current node */
	n = ++pr->count; /* start bubbling from the last node */
	p = n >> 1;
	while (p > 0) { /* bubble up the tree */
	    heap_find_pidx(pr, p, &p_pg, &p_idx);
            heap_find_pidx(pr, n, &n_pg, &n_idx);
	    if (p_pg[p_idx].subcuboid < t->subcuboid)
	        n_pg[n_idx] = p_pg[p_idx];
	    else break; /* rightful place found */
	    n = p; /* climb up the tree */
	    p = n >> 1;
       	}
	if (0 == t->id) t->id = pr->pid++; /* give unique id to particle */
        heap_find_pidx(pr, n, &n_pg, &n_idx);
	n_pg[n_idx] = *t;
}

@ Function |heap_remove_paged(pr,t)| removes a new particle |t| into the
particles repository |pr|. This removal uses the same method as
|heap_remove_fast(pr,t)|, except that it assumes heap paging; hence,
to reference, or change a node, we first find the correct page
and index for each node using |heap_find_pidx(pr,n,p,i)|.
@<Global functions@>=
void heap_remove_paged(ParticleRepository *pr, Particle *t)
{
	uint32_t n, l, r, q;  /* indices to node, left, right and last */
	Particle *n_pg, *l_pg, *r_pg, *q_pg; /* heap pages */
	uint32_t n_idx, l_idx, r_idx, q_idx; /* indices within heap pages */
	Particle temp;
	*t = pr->head[1]; /* particle to return */
	q = pr->count--;
	if (q) {
	    heap_find_pidx(pr, q, &q_pg, &q_idx);
	    temp = q_pg[q_idx]; /* the last node to be promoted */
	    n = 1; /* start at the root */
	    l = 2;
	    while (l < q) { /* bubble down the tree choosing left or
	   right subtree */
	        r = l + 1;		
		heap_find_pidx(pr, l, &l_pg, &l_idx);
		heap_find_pidx(pr, r, &r_pg, &r_idx);
		heap_find_pidx(pr, n, &n_pg, &n_idx);
		if (q_pg[q_idx].subcuboid < l_pg[l_idx].subcuboid) {
		    if (l_pg[l_idx].subcuboid < r_pg[r_idx].subcuboid) {
	    	        n_pg[n_idx] = r_pg[r_idx];
	    		n = r;
		    } else {
	    	        n_pg[n_idx] = l_pg[l_idx];
	    		n = l;
		    }
		} else {
    		    if (q_pg[q_idx].subcuboid < r_pg[r_idx].subcuboid) {
		    	n_pg[n_idx] = r_pg[r_idx];
			n = r;
		    } else break; /* rightful place found */
		}
		l = n << 1; /* go deeper until we have passed a leaf */
            }
	    heap_find_pidx(pr, n, &n_pg, &n_idx);
            n_pg[n_idx] = temp; /* place node */
	}
}

@ The following are the return codes that must be checked by the caller of the following functions.
@<Type definitions@>=
enum {
	HEAP_SUCCESS = 0,
	HEAP_EMPTY,
	HEAP_FULL,
	HEAP_ERROR_ALLOC,
	HEAP_ERROR_UNDEFINED
};

@ Function |heap_expand(r)| expands the heap |r| by adding a new page
at the end of the linked list.
@<Global functions@>=
int heap_expand(ParticleRepository *r)
{
	Particle* t = mem_typed_alloc(HEAP_PAGE_SIZE + 1, Particle, mem_p);
	if (NULL == t) return HEAP_ERROR_ALLOC;
	*(char **) t = NULL; /* make last page: `next' points to |NULL| */
	r->max += r->page_size;
	*((char **) r->tail) = (char *) t;
	r->tail = t;
	return HEAP_SUCCESS;
}

@ Function |heap_insert(r,p,e)| inserts a new particle |p| into the
particle repository |r|. If |e| is |true| and there is not enough
space in the heap, the heap will be expanded to fit |p|.  

The macro |heap_has_space(r)| may be used to check if there is empty
space in |r| before making an insertion call.
@d heap_has_space(r) ((r).count < (r).max)
@<Global functions@>=
int heap_insert(ParticleRepository *r, Particle *p, bool e)
{
	int i;
	if (r->count < r->max) {
	   if (r->count > r->page_size) heap_insert_paged(r, p); 
	   else heap_insert_fast(r, p);
	} else {
	   if (e) {
	       i = heap_expand(r);
	       if (i) return i;
	       heap_insert_paged(r, p);
	   } else return HEAP_FULL;
	}
	return HEAP_SUCCESS;
}

@ Function |heap_remove(r,t)| removes the particle at the top of
the max-heap representing the particles repository |r|, and stores the
particle into |p|.
@<Global functions@>=
int heap_remove(ParticleRepository *r, Particle *p)
{
	if (r->count == 0) return HEAP_EMPTY;
	if (r->count > r->page_size) heap_remove_paged(r, p);
	else heap_remove_fast(r, p);
	return HEAP_SUCCESS;
}

@ Function |heap_init(r)| initialises the max-heap |r| by creating the
first page.
@<Global functions@>=
int heap_init(ParticleRepository *r)
{
	r->head = mem_typed_alloc(HEAP_PAGE_SIZE + 1, Particle, mem_p);
	if (NULL == r->head) return HEAP_ERROR_ALLOC;
	r->count = 0;
	r->pid = 1;
	r->max = r->page_size = HEAP_PAGE_SIZE;
	*(char **) r->head = NULL;
	r->tail = r->head;
	return HEAP_SUCCESS;
}

@ Function |heap_print(f,r)| prints the heap |r| to the I/O stream
pointed to by |f|.
@<Global functions@>=
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
	    } else if (heap_insert(&particles, &p, true) != HEAP_SUCCESS) break;
	    heap_print(stdout, &particles);
	} while (1);
}

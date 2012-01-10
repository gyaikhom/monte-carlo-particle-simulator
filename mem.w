@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@*1 Memory management. Memory space must be allocated to store the
data that are required by the simulations. This include data
structures for representing the simulation world, the particles and
their trajectories, the physics tables etc. To reduce the number of
system calls required to allocate memory, and to simplify deallocation
of memory during clean-up, we borrow the memory management scheme
described by Donald E. Knuth in the book {\sl The Stanford
GraphBase}~[Addison-Wesley ({\bf 1993})].

With the concept of memory areas, we separate the conceptual
representation of the data from their actual storage in physical
memory. We store all of the raw data inside a linked list of memory
blocks, each block storing a collection of data of the same data
type. Within each data item, we store the information that are
required to process the higher-level conceptual representations.

@ Each memory area is represented by a variable of type |Area|. Before
using a memory area |s|, we must first initialise it by using the
\CEE/ macro |mem_init(s)|.

@d mem_init(s) *s = NULL
@<Type definitions@>=
typedef struct memory_area {
       char *first; /* pointer to the first usable location */
       struct memory_area *next; /* pointer to previously allocated block */
} *Area[1];

@ Function |mem_free(s)| frees the memory blocks currently allocated
under the memory area |s|. We start at the head of the linked list and
free each of the blocks. At the end of this function, the memory area
variable |s| will be in a state similar to the one when it was first
initialised.
@<Global functions@>=
void mem_free(Area s)
{
	Area t;
	while (*s) {
	      *t = (*s)->next;
	      free((*s)->first);
	      *s = *t;
	}
}

@ To allocate storage space to fit |n| elements of a given data type
|t| under the memory area |s|, we use the \CEE/ macro
|mem_typed_alloc(n,t,s)|.
@d mem_typed_alloc(n,t,s) (t*)mem_alloc((n) * sizeof(t), s)

@ Function |mem_alloc(n,s)| allocates a block of storage space under
the memory area |s| that will fit |n| consecutive bytes of data. In
addition to the raw data, a block must also allocate enough space to
store the two pointers defined by a memory area. These are used to
implement the linked list of memory blocks. The memory area pointers 
are stored at the end of each block. To improve efficiency of memory
allocation, we allocate a block as a consecutive array of 256
bytes.
@<Global functions@>=
char *mem_alloc(size_t n, Area s)
{
	size_t m = sizeof(char *); /* size of a pointer variable */
	char *block; /* address of the new block */
	Area t;
	@<Check if the requested size of the block is valid@>;
	@<Allocate space as a consecutive array of 256 bytes each@>;
	@<Add new block to the area by updating the linked list@>;
	return block;
}

@ For a request to be valid, the number of bytes |n > 0| and $|n| \le
(2^{16} - 1) - 2|m|$.

@<Check if the requested size of the block is valid@>=
if (1 > n || (0xffff00 - 2*m) < n) {
   return NULL;
}

@ We first round up, if necessary, the number of bytes |n| so that it
is a multiple of the size of a pointer variable. This ensures that
the area pointers are placed at a location |n| bytes after the
pointer returned by |calloc()|. No matter what that address is, the
placement after |n| bytes is unaffected as long as the size of a
pointer variable is a divisor of 256, since we allocate the space as
an array of 256 bytes each.
@d round_up_mult(x,y) ((((x) + (y) - 1) / (y)) * (y))
@<Allocate space as a consecutive array of 256 bytes each@>=
n = round_up_mult(n,m); /* round up |n| to a multiple of |m| */
block = (char *) calloc((unsigned) ((n + 2*m + 255) / 256), 256);

@ To update the linked list, we first locate the address of the area
pointers. Then, we append the new block at the beginning of the linked
list that is currently pointed by the memory area |s|.
@<Add new block to the area by updating the linked list@>=
if (block) {
   *t = (struct memory_area *) (block + n);
   (*t)->first = block;
   (*t)->next = *s;
   *s = *t;
}

@ To allocate a two-dimensional array which will fit $r \times c$
elements of a given data type |t| under the memory area |s|, we use
the \CEE/ macro |mem_typed_alloc2d(r,c,t,s)|.
@d mem_typed_alloc2d(r,c,t,s) (t**)mem_alloc_2d((r),(c),sizeof(t),s)

@ @<Global functions@>=
char **mem_alloc_2d(uint32_t r, uint32_t c, size_t s, Area a)
{
    char **p = NULL;
    uint32_t i;
    p = (char **) mem_alloc(r * sizeof(char *), a);
    if (NULL == p) return NULL;
    for (i = 0; i < r; ++i) {
	p[i] = mem_alloc(c * s, a);
        if (NULL == p[i])  return NULL;
    }
    return p;
}

@ We shall use only one memory area for the entire application.
@<Global variables@>=
Area mem_phase_one = {NULL}; /* initialisation */
Area mem_phase_two = {NULL}; /* persistent */

@** GPU Memory Management using Nvidia CUDA.
To use CUDA threads, we must make all of the relevant data structures
available on the device memory. Hence, similar to the host CPU, we
must implement a strategy for memory allocation and deallocation on
the GPU devices. We use memory areas, but with CUDA APIs for memory
allocation. Note that we use |cudaMalloc()|, hence, we must not forget
to initialise the memory using |cudaMemset()|, since there is no
interface similar to |calloc()|.
@<Type definitions@>=
typedef struct memory_area *GPUArea[1];

@
@d gpu_memblocks_per_hostblock 32 /* 512 bytes per block */
@<Global functions@>=
static void **next_gpu_memblock = NULL;
static void **bad_gpu_memblock = NULL;
void **create_gpu_memblock(GPUArea s) {
        void **slot = next_gpu_memblock;
	if (slot == bad_gpu_memblock) {
	   slot = mem_typed_alloc(gpu_memblocks_per_hostblock, void *, s);
	   if (NULL == slot) return NULL;
	   else {
	   	next_gpu_memblock = slot + 1;
		bad_gpu_memblock = slot + gpu_memblocks_per_hostblock;
	   }
	} else next_gpu_memblock++;
        return slot;
}

@ Function |cuda_mem_free(s)| frees the memory blocks in the GPU
device currently allocated under the memory area |s|. We start at the
head of the linked list and free each of the blocks. At the end of
this function, the memory area variable |s| will be in a state similar
to the one when it was first initialised.
@<Global functions@>=
void cuda_mem_free(GPUArea s)
{
	Area t;
	uint16_t i;
	while (*s) {
	      *t = (*s)->next;
	      for (i = 0; i < gpu_memblocks_per_hostblock && (*s)->first + i; ++i)
	      cudaFree(*(void **)((*s)->first + i));
	      free((*s)->first);
	      *s = *t;
	}
}

@ To allocate storage space to fit |n| elements of a given data type
|t| inside the GPU device under the memory area |s|, we use the \CEE/
macro |cuda_mem_typed_alloc(n,t,s)|.
@d cuda_mem_typed_alloc(n,t,s) (t*)cuda_mem_alloc((n) * sizeof(t), s)

@ Function |cuda_mem_alloc(n,s)| allocates a block of storage space
inside the GPU device under the memory area |s| that will fit |n|
consecutive bytes of data. In addition to the raw data, a block must
also allocate enough space to store the two pointers defined by a
memory area. These are used to implement the linked list of memory
blocks. The memory area pointers are stored at the end of each
block. To improve efficiency of memory allocation, we allocate a block
as a consecutive array of 256 bytes. The details of the implementation
are similar to that of |mem_alloc(n,s)|, so we will not repeat them
here.
@<Global functions@>=
void *cuda_mem_alloc(size_t n, GPUArea s)
{
	size_t m = sizeof(void *), q; /* size of a pointer variable */
	void *block, **t; /* address of the new block */
	if (1 > n || (0xffff00 - 2*m) < n) return NULL;
	n = round_up_mult(n,m); /* round up |n| to a multiple of |m| */
	q = (size_t) ((n + 255) / 256) * 256;
	if (cudaSuccess == cudaMalloc((void **)&block, q)) {
	    t = create_gpu_memblock(s);
	    if (NULL == t) {
	        cudaFree(block);
		return NULL;
	    }
	    cudaMemset(block, 0, q); /* get the |calloc()| effect */
	    *t = block;
	}
	return block;
}

@ To allocate a two-dimensional array inside the GPU device that will
fit $r \times c$ elements of a given data type |t| under the memory
area |s|, we use the \CEE/ macro
|cuda_mem_typed_alloc2d(p,r,c,t,s)|. The pitch value is returned in
|p|, which is the width in bytes of the two-dimensional allocation. It
is used to calculate the two-dimensional address of the array elements.

@d cuda_mem_typed_alloc2d(p,r,c,t,s)
(t*)cuda_mem_alloc_2d((p),(r),(c),sizeof(t),s)

@ @<Global functions@>=
void *cuda_mem_alloc_2d(size_t *p, uint32_t r, uint32_t c, size_t s, GPUArea a)
{
    void *block, **t;
    if (cudaSuccess == cudaMallocPitch((void**)&block, p, c * s, r)) {
	    t = create_gpu_memblock(a);
	    if (NULL == t) {
	        cudaFree(block);
		return NULL;
	    }
	    cudaMemset2D(block, *p, 0, c * s, r); /* get the |calloc()| effect */
	    *t = block;
	}
    return block;
}

@ We shall use only one GPU memory area for the entire application.
@<Global variables@>=
GPUArea mem_gpu = {NULL}; /* GPU memory area */

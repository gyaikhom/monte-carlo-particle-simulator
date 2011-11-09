@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Nvidia CUDA.
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

@ We shall use only one GPU memory area for the entire application.
@<Global variables@>=
GPUArea mem_gpu = {NULL}; /* GPU memory area */

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

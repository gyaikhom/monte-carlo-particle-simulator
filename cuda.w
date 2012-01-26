@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

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

@*1 Memory organisation.

Inside the GPU memory, we create a table which contains the geometry
information and the necessary physics tables. This has the following
structure:
@<Type definitions@>=
typedef struct gpu_table {
    uint32_t nc; /* number of entries in subcuboids table */
    uint32_t ns; /* number of entries in solids table */
    uint32_t np; /* number of entries in primitives table */
    uint32_t npb; /* number of entries in postfix expression buffer */
    uint32_t nsb; /* number of entries in solid indices buffer */
    uint32_t l, m, n; /* divisions along $x$, $y$ and $z$ axes */
    uint32_t nct; /* number of items in the subcuboids search tree */
    struct subcuboids_table_item *ctab; /* subcuboids lookup table */
    struct solids_table_item *s; /* solids lookup table */
    struct primitives_table_item *p; /* primitives lookup table */
    int32_t *pb; /* pointer to the postfix expression buffer */
    uint32_t *sb; /* pointer to the solid indices buffer */
    double *ctree; /* subcuboid search tree */
    int8_t *iltab; /* index lookup table for two-tiered neighbour table */
    uint32_t *ntab; /* subcuboid neighbour table */
    size_t pitch; /* width in bytes for the neighbour table */
    BoundingBox sw; /* the simulation world cuboid */
} GPUTables;

@ @<Global functions@>=
GPUTables *transfer_tables_to_gpu(GeometryTable *g)
{
    GPUTables *t, k;
    k.nc = g->nc;
    k.ns = g->ns;
    k.np = g->np;
    k.npb = g->npb;
    k.nsb = g->nsb;
    k.l = g->l;
    k.m = g->m;
    k.n = g->n;
    k.sw = g->sw;
    k.ctab = cuda_mem_typed_alloc(k.nc,struct subcuboids_table_item,mem_gpu);
    k.s = cuda_mem_typed_alloc(k.ns,struct solids_table_item,mem_gpu);
    k.p = cuda_mem_typed_alloc(k.np,struct primitives_table_item,mem_gpu);
    k.pb = cuda_mem_typed_alloc(k.npb,int32_t,mem_gpu);
    k.sb = cuda_mem_typed_alloc(k.nsb,uint32_t,mem_gpu);
    k.ctree = cuda_mem_typed_alloc(2 * (k.l + k.m +
        k.n),double,mem_gpu);
    k.iltab = cuda_mem_typed_alloc(MAX_SFIELD + 1,int8_t,mem_gpu);
    k.ntab = cuda_mem_typed_alloc2d(&k.pitch,k.nc,NUM_NEIGHBOURS,uint32_t,mem_gpu);
    t = cuda_mem_typed_alloc(1,GPUTables,mem_gpu);
    cudaMemcpy(k.ctab, g->ctab, k.nc * sizeof(struct
        subcuboids_table_item), cudaMemcpyHostToDevice);
    cudaMemcpy(k.s, g->s, k.ns * sizeof(struct solids_table_item),
        cudaMemcpyHostToDevice);
    cudaMemcpy(k.p, g->p, k.np * sizeof(struct primitives_table_item),
        cudaMemcpyHostToDevice);
    cudaMemcpy(k.pb, g->pb, k.npb * sizeof(int32_t),
        cudaMemcpyHostToDevice);
    cudaMemcpy(k.sb, g->sb, k.nsb * sizeof(uint32_t),
        cudaMemcpyHostToDevice);
    k.nct = 2 * (k.l + k.m + k.n);
    cudaMemcpy(k.ctree, g->ctree, k.nct *
        sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(k.iltab, g->iltab, (MAX_SFIELD + 1) *
        sizeof(int8_t), cudaMemcpyHostToDevice);
    cudaMemcpy2D(k.ntab, k.pitch, g->ntab, NUM_NEIGHBOURS *
        sizeof(uint32_t), NUM_NEIGHBOURS * sizeof(uint32_t),
        k.nc, cudaMemcpyHostToDevice);
    cudaMemcpy(t, &k, sizeof(k), cudaMemcpyHostToDevice);
    return t;
}

@ @<Global functions@>=
void transfer_tables_from_gpu(GPUTables *g, GPUTables *h)
{
    void *t;
    cudaMemcpy(h, g, sizeof(GPUTables), cudaMemcpyDeviceToHost);
    t = h->ctab;
    h->ctab = mem_typed_alloc(h->nc,struct subcuboids_table_item,mem_phase_two);
    cudaMemcpy(h->ctab, t, h->nc * sizeof(struct
        subcuboids_table_item), cudaMemcpyDeviceToHost);
    t = h->s;
    h->s = mem_typed_alloc(h->ns,struct
    solids_table_item,mem_phase_two);
    cudaMemcpy(h->s, t, h->ns * sizeof(struct solids_table_item),
        cudaMemcpyDeviceToHost);

    t = h->p;
    h->p = mem_typed_alloc(h->np,struct
    primitives_table_item,mem_phase_two);
    cudaMemcpy(h->p, t, h->np * sizeof(struct primitives_table_item),
        cudaMemcpyDeviceToHost);

    t = h->pb;
    h->pb = mem_typed_alloc(h->npb,int32_t,mem_phase_two);
    cudaMemcpy(h->pb, t, h->npb * sizeof(int32_t),
        cudaMemcpyDeviceToHost);

    t = h->sb;
    h->sb = mem_typed_alloc(h->nsb,uint32_t,mem_phase_two);
    cudaMemcpy(h->sb, t, h->nsb * sizeof(uint32_t),
        cudaMemcpyDeviceToHost);

    t = h->ctree;
    h->ctree = mem_typed_alloc(h->nct,double,mem_phase_two);
    cudaMemcpy(h->ctree, t, h->nct * sizeof(double),
        cudaMemcpyDeviceToHost);

    t = h->iltab;
    h->iltab = mem_typed_alloc(MAX_SFIELD + 1,int8_t,mem_phase_two);
    cudaMemcpy(h->iltab, t, (MAX_SFIELD + 1) * sizeof(double),
        cudaMemcpyDeviceToHost);

    t = h->ntab;
    h->ntab = mem_typed_alloc(h->nc * NUM_NEIGHBOURS,uint32_t,mem_phase_two);
    cudaMemcpy2D(h->ntab, NUM_NEIGHBOURS *
        sizeof(uint32_t), t, h->pitch, NUM_NEIGHBOURS * sizeof(uint32_t),
        h->nc, cudaMemcpyDeviceToHost);
}

@ @<Global functions@>=
void print_tables(FILE *f, GPUTables *g)
{
	uint32_t i, j;
	fprintf(f, "G E O M E T R Y  T A B L E\nPrimitives table:\n");
	for (i = 0; i < g->np; ++i) {
	    fprintf(f, "Affine:\n");
	    matrix_print(f, g->p[i].a, 4, 4, 0);
	    fprintf(f, "\nInverse:\n");
	    matrix_print(f, g->p[i].i, 4, 4, 0);
	    for (j = 0; j < 80; ++j) fprintf(f, "-");
	    fprintf(f, "\n");
	}
	fprintf(f, "Solids table:\n");
	for (i = 0; i < g->ns; ++i) fprintf(f, "%u %u\n", g->s[i].s, g->s[i].c);
	fprintf(f, "Postfix buffer:\n");
	for (i = 0; i < g->npb; ++i) fprintf(f, "%d ", g->pb[i]);
	fprintf(f, "\nSubcuboids table:\n");
	for (i = 0; i < g->nc; ++i) fprintf(f, "%u %u\n", g->ctab[i].s, g->ctab[i].c);
	fprintf(f, "Solid indices buffer:\n");
	for (i = 0; i < g->nsb; ++i) fprintf(f, "%u ", g->sb[i]);
	fprintf(f, "\n");
}

@ @<Test gpu tables@>=
{
	if (false == read_geometry("test/test_gpu_table.data")) exit(1);
	print_forest();
	create_geotab(&geotab);
	GPUTables *gpu, host;
	gpu = transfer_tables_to_gpu(&geotab);
	transfer_tables_from_gpu(gpu, &host);
	print_tables(stdout, &host);
}


@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

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


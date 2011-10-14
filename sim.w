@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Simulation. The particle simulations are carried out in batches of.

After a simulation has begun, it will only exit under two conditions:

1) All of the particles have been processed. This includes all of the
primary particles generated by the particle gun, and all of the
secondary particles generated by physics processes.

2) There was an error while processing the particles. There are
several causes for error, for instance, there is not enough space to
accommodate the heap expansion required for adding secondary
particles, etc.

@*1 Particle.
A {\sl particle}@^particle@> is the lowest-level data structure that
\.{MCS} manipulates. This data structure represents a fundamental
particle in physics, which interacts with the materials in the
world.

In \.{MCS}, all of the primary particles to be simulated are generated
by the particle gun. The origin of these particles are specified by the
vertex used by the particle gun. A particle gun can generate particles
from various vertices.

@<Type definitions@>=
struct particle_struct {
       @<Physical properties of a particle @>;
       @<Auxilliary data for managing a particle@>;
};

@ These are the properties used by the physics processes.

@<Physical properties of a particle @>=
Vector v; /* position */
Vector mo; /* momentum */
Vector po; /* polarisation */
double m; /* mass */
double c; /* charge */

@ Every particle knows which subcuboid it is currently inside. It also
maintains a unique identifier, which is supplied by the particle
repository. Finally, to determine the particle hierarchy, each
particle also maintains the identifier of its parent and a sibling
index, which ranks the children of a particle.

@<Auxilliary data for managing a particle@>=
uint32_t subcuboid; /* index of subcuboid containing particle */
uint32_t id; /* particle identifier */
uint32_t pi; /* parent identifier */
uint16_t si; /* sibling index (rank among siblings) */
uint16_t nd; /* number of daughter particles */
uint8_t s; /* the $s$-field */
uint8_t active; /* particle type */

@ @<Print particle information to |stdout|@>=
fprintf(stdout, "%u (%lf, %lf, %lf) %u (%lf, %lf, %lf) "
		"(%lf, %lf, %lf) %lf %lf %u %u %u\n",
	p->id, p->v[0], p->v[1], p->v[2], p->s,
	p->mo[0], p->mo[1], p->mo[2],
	p->po[0], p->po[1], p->po[2],
	p->m, p->c, p->pi, p->si, p->nd);

@ @<Print particle information to |gpfile|@>=
fprintf(gpfile, "%u (%lf, %lf, %lf) %u (%lf, %lf, %lf) "
		"(%lf, %lf, %lf) %lf %lf %u %u %u\n",
	gpbuff[i].id, gpbuff[i].v[0], gpbuff[i].v[1], gpbuff[i].v[2], gpbuff[i].s,
	gpbuff[i].mo[0], gpbuff[i].mo[1], gpbuff[i].mo[2],
	gpbuff[i].po[0], gpbuff[i].po[1], gpbuff[i].po[2],
	gpbuff[i].m, gpbuff[i].c, gpbuff[i].pi, gpbuff[i].si, gpbuff[i].nd);

@ Function |save_particle(p)| saves particle |p| in the buffer to be
written to the filesystem, evetually. When a particle is no longer
active, they must be transferred to the filesystem for later
analysis. To reduce the number of system calls, we use a particles
buffer to hold particles until they can be written in a batch write to
the file system. We maintain the buffer, the I/O stream and the number
of particles as global variables. These must be initialised
appropriately before calling |run_simulation()|.

@d MAX_PARTICLE_BUFFER 60 /* number of particles in output buffer */
@<Global functions@>=
Particle gpbuff[MAX_PARTICLE_BUFFER];
uint32_t gpbuffsize = 0; /* number of particles in buffer */
FILE *gpfile = NULL; /* the I/O stream to write particle data */
void save_particle(Particle *p)
{
	if (gpbuffsize == MAX_PARTICLE_BUFFER) @<Flush particles in buffer to filesystem@>;
	gpbuff[gpbuffsize++] = *p;
}

@ @<Flush particles in buffer to filesystem@>=
{
	uint32_t i;
	for (i = 0; i < gpbuffsize; ++i) @<Print particle information to |gpfile|@>;
	gpbuffsize = 0;
}

@*1 Vertex.
To generate particles, a particle gun must first be placed inside the
world by choosing a vertex. When the particle gun is activated,
primary particles are generated from this vertex. A vertex is not
associated with a particle gun, however; they are stored with respect
to an event. Hence, an event can have multiple vertices associated
with it. These are stored as a linked-list.

@f uint32_t int
@<Type definitions@>=
typedef struct vertex_struct {
       Vector v; /* particle gun position vector */
       Vector mo; /* momentum */
       Vector po; /* polarisation */
       double e; /* energy */
       double c; /* charge */
       uint32_t np; /* number of particles required */
       struct vertex_struct *next;
} Vertex;

@ Function |create_vertex(e)| creates a new vertex under the supplied
event |e|. This functions returns a pointer to the vertex, or |NULL|
if a new vertex could not be created.

@d vertices_per_block 10 /* TODO: */
@<Global functions@>=
static Vertex *next_vertex = NULL;
static Vertex *bad_vertex = NULL;
Vertex *create_vertex(Event *e)
{
        Vertex *v = next_vertex;
	if (v == bad_vertex) {
	   v = mem_typed_alloc(vertices_per_block, Vertex, mem_p);
	   if (NULL == v) return NULL;
	   else {
	   	next_vertex = v + 1;
		bad_vertex = v + vertices_per_block;
	   }
	} else ++next_vertex;
	@<Insert vertex to event@>;
        return v;
}

@ @<Insert vertex to event@>=
e->nv++;
v->next = e->v;
e->v = v;

@*1 Events.
An {\sl event}@^event@> is the highest-level simulation object. It
provides a link between the user and the simulator. Users specify the
number of events they wish to simulate. An event is associated with
various vertices which specify the locations where a particle gun
could be placed. To generate primary particles, the particle gun must
first be placed on one of the vertices of the event being simulated.

@<Type definitions@>=
typedef struct event_struct {
       uint32_t id; /* event identifier */
       uint32_t nv; /* number of vertices */
       Vertex *v; /* head of the vertex linked list */
       struct event_struct *next; /* pointer to next event */
} Event;

@ Function |create_event()| creates a new event with a unique
identifier. This functions returns a pointer to the event, or |NULL|
if a new event could not be created.

@d events_per_block 10 /* TODO: */
@<Global functions@>=
static Event *next_event = NULL;
static Event *bad_event = NULL;
static uint32_t next_event_id = 0;
Event *create_event()
{
        Event *e = next_event;
	if (e == bad_event) {
	   e = mem_typed_alloc(events_per_block, Event, mem_p);
	   if (NULL == e) return NULL;
	   else {
	   	next_event = e + 1;
		bad_event = e + events_per_block;
	   }
	} else ++next_event;
	e->id = ++next_event_id;
        return e;
}

@*1 Particle gun.

@<Generate primary particle using particle gun@>=
vector_zero(p.v);
vector_zero(p.mo);
vector_zero(p.po);
p.m = p.c = 0.0;
p.id = p.pi = p.si = 0;
p.active = true; /* only active particles inside heap */
p.s = 0x0; /* renew: $s$-field are changed by physics processes */
p.nd = 0; /* no daughters yet */
p.subcuboid = find_subcuboid(subcuboid_search_tree, p.v);

@*1 Batch simulation.

@ Function |generate_primaries()| fills up the particles repository
with primary partices generated using a particle gun. It returns zero
if all of the primary particles required by the event vertices have
already been processed; otherwise, the number of new primaries are
returned. The batch simulation loop will continue to be valid as long
as the return value of this function is nonzero.

While generating primaries, the aim is to fill up the max-heap without
expanding the heap, so that only primary particles that are needed to
fill up a batch are generated. Of course, since this function will be
called multiple times by |get_particle()| while creating a simulation
batch, all of the particles that are required by all of the vertices
in all of the events will be processed, eventually.

During multiple calls to |generate_primaries()|, we maintain from
one call to the next the current event, the current vertex within that
event, and the number of particles already generated for the current
vertex. Hence, before running a simulation, i.e., calling
|run_simulation()|, these global variables must be initialised
appropriately.

@<Global functions@>=
Event *ge = NULL; /* current simulation event */
Vertex *gv = NULL; /* current vertex within simulation event */
uint32_t gc = 0; /* current particle count for current vertex */
int generate_primaries()
{
    Particle p; /* primary particle to be generated */
    int n = 0; /* number of particles generated in this call */
    while (ge) {
        while (gv) {
             while (gc < gv->np) {
	         if (heap_has_space(particles)) {
                     @<Generate primary particle using particle gun@>;
   		     heap_insert(&particles, &p, false);
   		     ++n;
   		     ++gc;
		 } else goto heap_is_full;
             }
	     gv = gv->next; /* move to the next vertex */
	     gc = 0;
	}
	ge = ge->next; /* move to the next event */
	if (ge) { 
	    gv = ge->v;
	    gc = 0;
	}
    }
    heap_is_full:
    return n;
}

@ Function |get_particle(p)| retrieves a particle from the particle
repository and stores the particle in |p|. It returns 1 if a particle
was retrieved successfully; otherwise, returns 0. If there are no
particles in the particles repository, it will first attempt to fill
up the particles repository, by generating primary particles using the
particle gun, and then retry to retrieve particle from the heap.

@<Global functions@>=
int get_particle(Particle *p)
{
	if (HEAP_EMPTY == heap_remove(&particles, p)) {
	    if (generate_primaries())
	        heap_remove(&particles, p); /* try again */
	    else return 0; /* no more particles left */
	}
	return 1; /* one particle retrieved */
}

@ All of the primary particles that was added from the particles
repository are stored in the particles array |p|, whereas, the
particles generated as secondaries by the compute threads as a result
of physics processes are stored in the particles array |s|. Before the
compute threads have worked on the block, the number of secondaries
|ns| are always set to zero by the function |create_batch()|. After
the block has been processed, the primary particles array |p| could be
empty; however, the secondary particles array |s| will never be empty.

@d MAX_PRIMARIES_BLOCK 4 /* must depend on available memory */
@d MAX_DAUGHTERS 4 /* maximum number of daughters allowed */
@d MAX_SECONDARIES_BLOCK (MAX_DAUGHTERS * MAX_PRIMARIES_BLOCK)
@<Type definitions@>=
typedef struct block_struct {
	uint32_t subcuboid; /* index of the associated subcuboid */
	uint16_t np; /* number of primaries */
	uint16_t ns; /* number of secondaries */
	Particle p[MAX_PRIMARIES_BLOCK];
	Particle s[MAX_SECONDARIES_BLOCK];
} Block;

@ Each simulation batch comprises of several simulation blocks that
contain data about the particles and the subcuboid which contains
those particles. Each of these blocks could be processed exclusively
by a single thread, or a group of compute threads.

We have used fixed-size arrays, instead of dynamically allocated
memory blocks, so that the input data as they are available on the
host application, or the output data as they are available in the
compute memory, could be trasferred using a single memory copy of the
batch. For instance, if we wish to use multiple threads that are
available on graphics processing units, we can easily transfer the
batches back and forth using a single memory copy between the host
application and the graphics hardware. This is also important because,
for cache efficiency, we required the memory footprint on the compute
devices to be quite simple and compact.

@d MAX_BLOCKS_BATCH 4 /* must depend on available compute blocks */
@<Type definitions@>=
typedef struct sim_batch_struct {
	uint16_t n; /* number of blocks with particles in them */
	Block b[MAX_BLOCKS_BATCH]; /* array of blocks */
} Batch;

@ Function |create_batch(b)| creates a simulation batch and stores the
batch in |b|. If successful, it returns a positive number that gives
the number of particles in the batch; otherwise, it returns zero to
mark a successful end of simulation.

A simulation batch is created by starting at the first slot of the
first block. We continue filling each slot in each of the blocks by
taking particles from the particles repository, until all of the
blocks are full. Furthermore, since all of the particles in a block
must be contained by the same subcuboid that is associated with that
block, we move forward to the next block when the particle retrieved
from the particles repository has a smaller subcuboid index than the
subcuboid index associated with the current block. Since the particles 
repository is implemented as a max-heap, where the subcuboid index is
used as a comparison key, the |get_particle()| call will never return
a particle with subcuboid index greater than the subcuboid index of
the current block. The blocks, therefore, effectively group the
particles that belong in the same subcuboid.

@<Global functions@>=
int create_batch(Batch *b)
{
	Particle p; /* particle retrieved from heap */
	uint32_t n; /* number of particles in batch */
	uint32_t i, j; /* current block, and current slot within block */
	i = j = 0; /* start at first slot of first block */
	b->n = n = 0; /* reset block and particle count */
	if (get_particle(&p)) {
	    b->b[i].subcuboid = p.subcuboid; /* subcuboid index for the block */
	    b->b[i].p[j++] = p; /* add particle to block */
	    b->n++; /* increment number of active blocks */
	    ++n;
	    while (get_particle(&p)) {
	        if (b->b[i].subcuboid != p.subcuboid ||
		    MAX_PRIMARIES_BLOCK == j) { /* change block */
		    b->b[i].np = j; /* finalise current block */
		    b->b[i].ns = 0;
		    if (MAX_BLOCKS_BATCH == ++i) { /* move to next block */
		        heap_insert(&particles, &p, false); /* batch is full: retain particle for next batch */
		        goto batch_is_full;
		    }
		    b->b[i].subcuboid = p.subcuboid; /* subcuboid index for the new block */
                    b->n++; /* increment number of active blocks */
		    j = 0; /* first slot of new block */
	        }
	        b->b[i].p[j++] = p; /* add particle to block */
		++n;
	    }
    	    b->b[i].np = j; /* finalise current block */
	    b->b[i].ns = 0;
    	}
batch_is_full:
	return n;
}

@ Function |update_repository(b)| updates the max-heap, following a
successful simulation of the batch |b|, by moving primary and
secondary particles from the active blocks in that batch to the
particles repository.

@<Global functions@>=
int update_repository(Batch *b)
{
    uint16_t i, j;
    Particle p;
    for (i = 0; i < b->n; ++i) {
	@<Move unprocessed primary particles to repository@>;
	@<Process particles in the secondary particles array after a batch simulation@>;
    }
    return 0;
}

@ It is unnecessary and inefficient to move primary particles from the
block back to the particles repository. However, for now, we take the
simpler route. Future changes must attempt to update the blocks
without having to removed unprocessed primary particles that are
already in the block.

@<Move unprocessed primary particles to repository@>=
for (j = 0; j < b->b[i].np; ++j)
    heap_insert(&particles, &b->b[i].p[j], true);

@ After a block has been processed by the compute threads, all of the
particles in the block's secondary particles array |s| are of three
types: 1) primary particles that are no longer active, 2) primary
particles that are still active, and 3) daughter particles generated
by the physics processes. All of these particles must first be saved
to the file system for later analysis. In the second and third cases,
we prepare the particle for further processing, by updating its
containing subcuboid using the $s$-field. After this, we only move the
particle to the particles repository if it continues to stay inside
the simulation world.

@<Process particles in the secondary particles array after a batch simulation@>=
for (j = 0; j < b->b[i].ns; ++j) {
    p = b->b[i].s[j];
    save_particle(&p);
    if (p.active) {
        @<Prepare particle for further processing@>;
	if (OUTSIDE_WORLD == p.subcuboid) continue;
	heap_insert(&particles, &p, true);
    }
}

@ After a primary particle has been processed, or a secondary particle
has been generated by physics processes, the particle's position
vector |v| reflects its new location. Using this new location, and the 
bounding box |bb| of the containing subcuboid, we must first update the
particles $s$-field. The $s$-field is then used in combination with
the index of the current subcuboid to determine the index of the new
subcuboid that will now contain the particle. Once this is done,
we reset the $s$-field for the next simulation batch.

@<Prepare particle for further processing@>=
update_sfield(&(p.s), &(subcuboids[p.subcuboid].bb), p.v);
p.subcuboid = get_neighbour(p.subcuboid, p.s);
p.s = 0x0; /* renew particle */
p.id = 0; /* get a new identifier from particle repository */

@ @<Global functions@>=
int simulate_batch(Batch *b)
{
	uint32_t i = 0, j, k;
	while(i < b->n) {
	    j = 0;
	    k = 0;
	    while(j < b->b[i].np) {
	        b->b[i].s[k] = b->b[i].p[j];
		b->b[i].s[k].active = false;
		++j;
		++k;
	    }
	    b->b[i].ns = k;
	    b->b[i].np = 0;
	    ++i;
	}
	return 0;
}

@ Function |run_simulation()| runs the simulation iteratively by 
creating and simulating batches in each iteration, until all of the
particles, both primary and secondary, have been processed. It returns
0 if the simulation was a success; otherwise, a non-zero value is returned.
@<Global functions@>=
int run_simulation()
{
	int i = 0;
	Batch b;
	while ((i = create_batch(&b))) {
	    if (i < 0) goto handle_error;
	    if (simulate_batch(&b)) goto handle_error;
	    else if (update_repository(&b)) goto handle_error;
	}
	return 0; /* simulation done */

handle_error:
	fprintf(stderr, "Error\n");
	return 1;
}

@ @<Test simulation batch@>=
{
	Event *e, *ee;
	Vertex *v, *vv;

	/* first event */
	ee = create_event();
	vv = create_vertex(ee); vv->np = 5;
	v = create_vertex(ee); v->np = 5;
	v = create_vertex(ee); v->np = 5;
	v = create_vertex(ee); v->np = 5;

	/* second event */
	e = create_event();
	v = create_vertex(e); v->np = 5;
	v = create_vertex(e); v->np = 5;
	v = create_vertex(e); v->np = 5;
	v = create_vertex(e); v->np = 5;	
	e->next = ee;
	ee = e;

	/* third event */
	e = create_event();
	v = create_vertex(e); v->np = 5;
	v = create_vertex(e); v->np = 5;
	v = create_vertex(e); v->np = 5;
	v = create_vertex(e); v->np = 5;	
	e->next = ee;
	ee = e;

	e = ee;
	while (e) {
	    v = e->v;
	    while (v) {
	        fprintf(stderr, "%u %u\n", e->id, v->np);
	        v = v->next;
	    }
	    e = e->next;
	}
	gpfile = stdout;
	BoundingBox bb = {{100.0, 100.0, 1000.0, 1.0},{0.0, 0.0, 0.0, 1.0}};
	build_subcuboid_trees(&bb, 4, 4, 4);
	heap_init(&particles);
	ge = ee;
	if (ge) {
	    gv = ge->v;
	    gc = 0;
        }
	run_simulation();
	@<Flush particles in buffer to filesystem@>;
}

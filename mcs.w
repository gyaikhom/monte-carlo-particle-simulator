\input miniltx 
\input epsfig.sty 
\resetatcatcode

\def\hash{\char35}

@* Introduction. The {\sl Monte Carlo Simulator}, \.{MCS} in short, is
a simulation package for high-energy physics. It uses Monte Carlo
techniques for the randomised simulation of particle interactions. The
system borrows concepts defined in the Geant4 system; however, they
have been reimplemented using ANSI \CEE/ using new data structures and
algorithms with the aim that the system will be ported eventually to
GPGPUs.

The \.{MCS} system is in effect an event processor. The user specifies
the number of events they wish to simulate, and provides the system
with the necessary details to process each of the events. These
details include procedures and data for generating the particles for
each event, and the {\sl physics processes}@^physics processes@>
required by the type of the particle being simulated and the
properties of the materials involved in each of the interactions. We
are only interested in interactions that are happening inside a closed
three-dimensional system, referred to as the {\it
world}@^world@>. Since the world is defined by the components present,
the user must provide the geometry of the components in addition to
the properties of the materials the components are made of. In this
document, these input data will be referred to {\sl simulation
parameters}@^simulation parameters@>.

Once the \.{MCS} system is supplied with a valid set of simulation
parameters, it processes the geometries and creates a data structure
that allows efficient location of an interaction point. This is an
important step because, to process a large number of interactions
rapidly, we must efficiently find the location inside the world where
the interactions are happening, as it determines the materials
required by the processes.

We then define the {\sl particle gun}@^particle gun@>, which generates
the required number of {\sl primary particles}@^primary particles@>
for each of the events. Primary particles are particles that are
generated explicitly by a particle gun, and {\sl secondary
particles}@^secondary particles@> are those that are generated as a
result of particle-matter interaction. Since multiple primary
particles may be generated by an event, a {\sl vertex}@^vertex@> is
defined for each event. For any given event the associated vertex
provides the location of the particle gun, so that all of the
particles generated originates from this location. When the particle
gun is activated, a set of primary particles are generated. These
particles have the same origin, as defined by the vertex location,
however, they carry different properties, such as energy, momentum,
etc. These values are set using a {\sl random number
generator}@^random number generator@>.

Once a particle has been generated, it travels through the materials
until it either comes to a stop, or exits the closed three-dimensional
world. Throughout this journey, the {\sl trajectory}@^trajectory@> of
the particle changes depending on its interaction with the
materials. Hence, the \.{MCS} system uses a {\sl tracker}@^tracker@>
to track each of the particles throughout its journey.
A snapshot of the particle in its trajectory is referred to as a
{\sl track}@^track@>, and is processed independently of any previous
tracks. In other words, a track does not remember its past. Hence, to
keep a record of the trajectory, the tracker records old tracks as new
ones are derived. Each of these points on the trajectory is referred
to as the {\sl step point}@^step point@>, and the process of moving to
the next step point as {\sl stepping}@^stepping@>. Hence, we invoke
multiple stepping commands to chart the entire trajectory of a
particle.

To move a particle from one track to the next, each stepping decides
the length with which the particle must progress. This is referred to
as a {\sl step}@^step@>. Each step has a start point and an end
point. These points store information for retrieving the material
properties that are required by each of the physics processes that are
valid in that stepping. The length of a step is determined by the
process which requires the shortest space-time {\sl interaction
length}@^interaction length@>. In other words, the next track must be
located in space-time so that all of the valid physics processes are
applicable.

We use a {\sl stack}@^stack@> to store all of the unprocessed
particles, both primary and secondary. At the beginning of each event,
each of the primary particles generated by the particle gun are first
converted to a track, so that the tracker can process them. Then, the
tracks are all pushed into the {\sl track stack}@^track stack@>. This
stack is then passed to the tracker, which continuously pops a track
from the track stack and charts the particle's trajectory. If
secondary particles are generated as a result of an interaction, these
are first converted to a track, and then pushed into the track
stack. In other words, if we consider the entire simulation as a tree,
where the root of the tree represents the particle gun, the nodes
represent the particle-matter interactions, and the leaves represent
interactions which resulted in the particle to either stop, or exit
the closed three-dimensional world, then the tracker is simply
undergoing a {\sl depth-first tree traversal}@^depth-first traversal@>.

@* History.
The {\sl Monte Carlo Simulator} project began in June 2011, when
Dr.~Gagarine Yaikhom was a WIMCS Research Fellow under Prof.~David
W. Walker at Cardiff University. The initial aim of the project was to
port the Geant4 system to run on CUDA GPUs. However, after studying
the Geant4 code-base, we concluded that it will be prohibitive to port
the entire Geant4 system given the short duration of the funding (8
months, from June 2011 until January 2012) and that only Dr.~Yaikhom
will be carrying out the design and implementation. We, therefore,
decided to use the Geant4 system only as a guideline system
architecture, and to reimplement the performance-intensive concepts
and their dependencies for a simplified simulator. This
implementation, which began on 5 August 2011, uses data structures
and algorithms that are more appropriate for eventual multithreaded
parallelisation on GPUs.

@i types.w
@i common.w
@i csg.w
@i sim.w
@i main.w

@**Index.


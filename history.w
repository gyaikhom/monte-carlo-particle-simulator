@* History.
The {\it Monte Carlo Simulator} project began in June 2011, when
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

@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@* Error. This section defines data-structures and functions for
handling and reporting errors.

There are three message categories, which are printed using the
following macros: |fatal|, |warn|, and |info|.

@d fatal(...) fprintf(stderr, __VA_ARGS__)
@d warn(...) if (verbose) {
     fprintf(stderr, __VA_ARGS__ );
}
@d info(...)  if (verbose) {
     fprintf(stderr, __VA_ARGS__ );
}
@<Global variables@>=
bool verbose = false; /* verbose output */

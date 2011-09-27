@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@* Error. This section defines data-structures and functions for
handling and reporting errors.

There are three message categories, which are printed using the
following macros: |fatal|, |warn|, and |info|.

@d fatal(X) fprintf(stderr, "%s[%5d] %s\n", __FILE__, __LINE__, X);
@d warn(X) fprintf(stderr, "%s[%5d] %s\n", __FILE__, __LINE__, X);
@d info(X) fprintf(stderr, "%s[%5d] %s\n", __FILE__, __LINE__, X);

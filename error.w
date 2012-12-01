@q This file is part of the Monte Carlo Simulator (c) G. Yaikhom, Cardiff University 2011, 2012 @>

@* Error messaging. This section defines data-structures and functions for
handling and reporting errors. There are three message categories,
which are printed using the following macros:

\begin{itemize}
\item We use the |fatal| macro to advise the user of irrecoverable
errors. They are usally commnicated before the application is about to
exit due to the errors. Fatal messages must provide enough information
to help the user diagnose the issues that caused the exit. Messages
communicated with the |fatal| macro should never is suppressed, even
when |verbose == true|.

\item We use the |warn| macro to advise the user of non-fatal
recoverable issues, such as wrong input data. They are usually
communicated without exiting the application. Messages communicated
with the |warn| macro are suppressed when |verbose == false|.

\item We use the |info| macro to advise the user on the general state
of the application, for instance, the current stage in the
processing. They are usually communicated without exiting the
application. Messages communicated with the |info| macro are
suppressed when |verbose == false|.
\end{itemize}

NOTE:
Some of the messaging facilities are defined as {\sl variadic macros},
using the ISO C99 syntax. Please ensure that your compiler fully
supports variadic macros; otherwise, redefine these macros as functions.

@<Global variables@>=
#define fatal(...)@+ fprintf(stderr, __VA_ARGS__)
#define warn(...)@+ if (verbose) fprintf(stderr, __VA_ARGS__)
#define info(...)@+ if (verbose) fprintf(stderr, __VA_ARGS__)
bool verbose = false; /* verbose output if true */

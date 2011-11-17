@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@** Common. This section defines common types, physical constants,
data-structures and functions that are used in various sections. 

@ We specify Boolean variables using the |bool| data type.
@<Type definitions@>=
typedef enum {false = 0, true} bool;

@ The following are mathematical constants that are use throughout.

@d PI 3.14159265358979323846
@d TWICE_PI 6.283185307
@d RADIAN_TO_DEGREE (180.0 / PI) /* $|RADIAN_TO_DEGREE| = 180 / \pi$ */
@d DEGREE_TO_RADIAN (PI / 180.0)

@ Function |convert_radian_to_degree(angle)| determines the positive
angle in degrees that is equivalent to the specified |angle| in
radians.
@<Global functions@>=
double convert_radian_to_degree(double angle)
{
      angle *= RADIAN_TO_DEGREE;
      if (angle < 0.0) angle += 360.0; /* positive angle required */
      return angle;
}

@(mcs.cu@>=
__device__ double cuda_convert_radian_to_degree(double angle)
{
      angle *= RADIAN_TO_DEGREE;
      if (angle < 0.0) angle += 360.0; /* positive angle required */
      return angle;
}

@ Function |roundup_pow2(v)| rounds up a positive number to the
nearest power of two.

@<Global functions@>=
unsigned roundup_pow2(unsigned v)
{
        if (v == 0) return 1;
        if (v & (v - 1)) { /* check if already power of two */
                unsigned k;
		@<Find non-zero most significant bit@>;
                return (1 << (k + 1));
        }
        return v;
}

@ @<Find non-zero most significant bit@>=
for (k = sizeof(unsigned) * 8 - 1; ((1 << k) & v) == 0; --k);

@i vector.w
@i matrix.w
@i mem.w

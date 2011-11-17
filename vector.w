@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@*1 Vector. This section defines data-structures and functions for
representing and manipulating three-dimensional vectors.

In {\tt MCS}, we use the {\sl homogeneous coordinate
system}@^Homogeneous coordinate system@> introduced by August
Ferdinand M\"obius@^August Ferdinand M\"obius@> in his work {\sl Der
barycentrische Calc\"ul} [Leipzig, {\bf 1827}]. By expressing vectors
using homogeneous coordinates, we simplify {\sl affine
transformations}@^Affine transformations@> of a vector, such as
translation, rotation and scaling, as these can be expressed
efficiently as a left-hand multiplication of the vector by a $4 \times
4$ affine transformation matrix.

The homogeneous coordinates of a finite point $(x', y', z')$ in the
three-dimensional cartesian volume are any four numbers $(x, y,
z, w)$ for which $x' = x/w$, $y' = y/w$, and $z' =
z/w$. Before applying any affine transformation matrix to a vector, we
must ensure that $w = 1$. We shall refer to this as the $w$-component.

@<Type definitions@>=
typedef double Vector[4];

@ Function |vector_print(f,v)| prints the components of vector
|v| to the I/O stream pointed to by |f|. Vector |v| is left
unmodified.
@<Global functions@>=
void vector_print(FILE *f, Vector v)
{
	fprintf(f, "(%lf, %lf, %lf, %lf)", v[0], v[1], v[2], v[3]);
}

@ Function |vector_zero(v)| modifies the vector |v| by setting all of
its components to zero, except for its $w$-component, which is set
to 1.

@d ZERO_VECTOR { 0.0, 0.0, 0.0, 1.0 }
@<Global functions@>=
void vector_zero(Vector v)
{
	v[0] = v[1] = v[2] = 0.0;
	v[3] = 1.0;
}

@ @(mcs.cu@>=
__device__ void cuda_vector_zero(Vector v)
{
	v[0] = v[1] = v[2] = 0.0;
	v[3] = 1.0;
}

@ Function |vector_homogenise(v)| modifies the vector |v|, so that its
$w$-component equals 1.
@<Global functions@>=
void vector_homogenise(Vector v)
{
	if (1.0 == v[3]) return; /* already homogenised */
	if (0.0 != v[3]) {
	        v[0] /= v[3];
	        v[1] /= v[3];
	        v[2] /= v[3];
		v[3] = 1.0;
        }
}

@ @(mcs.cu@>=
__device__ void cuda_vector_homogenise(Vector v)
{
	if (1.0 == v[3]) return; /* already homogenised */
	if (0.0 != v[3]) {
	        v[0] /= v[3];
	        v[1] /= v[3];
	        v[2] /= v[3];
		v[3] = 1.0;
        }
}


@ Function |vector_copy_f(u,v)| modifies vector |u| by copying all of
the components of vector |v| to vector |u|. Vector |v| is left
unmodified. We may prefer to use the macro |vector_copy(X, Y)| instead.
@d vector_copy(X, Y) memcpy((X), (Y), 4*sizeof(double))
@<Global functions@>=
void vector_copy_f(Vector u, Vector v)
{
	u[0] = v[0];
        u[1] = v[1];
        u[2] = v[2];
	u[3] = v[3];
}

@ @(mcs.cu@>=
__device__ void cuda_vector_copy(Vector u, Vector v)
{
	u[0] = v[0];
        u[1] = v[1];
        u[2] = v[2];
	u[3] = v[3];
}


@ Function |vector_magnitude(v)| returns the {\sl vector magnitude} of
the vector |v|. Vector |v| is left unmodified.@^vector magnitude@>
@<Global functions@>=
double vector_magnitude(Vector v)
{
	return sqrt(v[0] * v[0] +
		    v[1] * v[1] +
		    v[2] * v[2]);
}

@ @(mcs.cu@>=
__device__ double cuda_vector_magnitude(Vector v)
{
	return sqrt(v[0] * v[0] +
		    v[1] * v[1] +
		    v[2] * v[2]);
}


@ Function |vector_normalize(v,r)| {\sl normalises} the vector |v| and
stores the result in |r|. Vector |v| is left unmodified. The
normalised vector |r| is an {\sl unit vector} with magnitude 1. This
function assumes that the $w$-component of |v| is already 1.
@^unit vector@>
@^vector normalisation@>
@<Global functions@>=
void vector_normalize(Vector v, Vector r)
{
	double m = vector_magnitude(v);
	/* assumes r[3] = 1.0 */
	if (m > 0.0) {
		r[0] = v[0] / m;
		r[1] = v[1] / m;
		r[2] = v[2] / m;
	} else {
		r[0] = 0.0;
		r[1] = 0.0;
		r[2] = 0.0;
	}
}

@ @(mcs.cu@>=
__device__ void cuda_vector_normalize(Vector v, Vector r)
{
	double m = cuda_vector_magnitude(v);
	/* assumes r[3] = 1.0 */
	if (m > 0.0) {
		r[0] = v[0] / m;
		r[1] = v[1] / m;
		r[2] = v[2] / m;
	} else {
		r[0] = 0.0;
		r[1] = 0.0;
		r[2] = 0.0;
	}
}

@ Function |vector_difference(u,v,r)| calculates the {\sl vector
difference} by subtracting vector |v| from vector |u|, and stores the
result in vector |r|. Both vectors |u| and |v| are left unmodified. This
function assumes that the $w$-components in |v| and |u| are already 1.
@^vector difference@>
@<Global functions@>=
void vector_difference(Vector u, Vector v, Vector r)
{
	r[0] = u[0] - v[0];
	r[1] = u[1] - v[1];
	r[2] = u[2] - v[2];
	r[3] = 1.0;
}

@ @(mcs.cu@>=
__device__ void cuda_vector_difference(Vector u, Vector v, Vector r)
{
	r[0] = u[0] - v[0];
	r[1] = u[1] - v[1];
	r[2] = u[2] - v[2];
	r[3] = 1.0;
}


@ Function |vector_dot(u,v)| returns the {\sl dot} (also known as the
{\sl scalar}) product of the two vectors |u| and |v|. Vectors |u| and |v|
are left unmodified. Vector dot products are {\sl
commutative}@^commutative@>---the order of the parameters are
irrelevant.@^dot product@>@^scalar product@>
@<Global functions@>=
double vector_dot(Vector u, Vector v)
{
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

@ @(mcs.cu@>=
__device__ double cuda_vector_dot(Vector u, Vector v)
{
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

@ Function |vector_cross(u,v,r)| calculates the {\sl cross} product
between vector |u| and vector |v|, and stores the result in vector
|r|. Vectors |u| and |v| are left unmodified.

Vector |r| is perpendicular to both |v| and |u|, with a direction
given by the {\sl right-hand rule}@^right-hand rule@> (when |u| and
|v| respectively point towards the index finger and the middle finger,
|r| points in the direction of the thumb). The {\sl magnitude}
@^vector magnitude@> of |r| is equal to the area of the parallelogram 
that the vectors span from |u| towards |v|. Vector cross products are
{\sl noncommutative}@^noncommutative@>---the order of the parameters are
significant.
@^cross product@>
@<Global functions@>=
void vector_cross(Vector u, Vector v, Vector r)
{
	r[0] = (u[1] * v[2] - u[2] * v[1]);
	r[1] = (u[2] * v[0] - u[0] * v[2]);
	r[2] = (u[0] * v[1] - u[1] * v[0]);
}

@ @(mcs.cu@>=
__device__ void cuda_vector_cross(Vector u, Vector v, Vector r)
{
	r[0] = (u[1] * v[2] - u[2] * v[1]);
	r[1] = (u[2] * v[0] - u[0] * v[2]);
	r[2] = (u[0] * v[1] - u[1] * v[0]);
}


@ Function |vector_angle_radian(u,v)| returns the angle $\theta$ in
radians between the vectors |u| and |v|, where rotation started
at vector |u|. The value of $\theta$ is in the range $[0, 2\pi]$.
Vectors |u| and |v| are left unmodified.

@d TWICE_PI 6.283185307
@<Global functions@>=
double vector_angle_radian(Vector u, Vector v)
{
        Vector a, b, c = ZERO_VECTOR;
	double angle;
        vector_normalize(u, a);
	vector_normalize(v, b);
        angle = acos(vector_dot(a, b));
        vector_cross(u, v, c);
	if (c[3] < 0.0) return (TWICE_PI - angle);
	return angle;
}

@ @(mcs.cu@>=
__device__ double cuda_vector_angle_radian(Vector u, Vector v)
{
        Vector a, b, c = ZERO_VECTOR;
	double angle;
        cuda_vector_normalize(u, a);
	cuda_vector_normalize(v, b);
        angle = acos(vector_dot(a, b));
        cuda_vector_cross(u, v, c);
	if (c[3] < 0.0) return (TWICE_PI - angle);
	return angle;
}


@ Function |vector_angle_degree(u,v)| returns the angle $\theta$ in
degrees between the vectors |u| and |v|, where rotation began
at vector |u|. The value of $\theta$ is in the range $[0^\circ, 
360^\circ]$.  Vectors |u| and |v| are left unmodified.
@<Global functions@>=
double vector_angle_degree(Vector u, Vector v)
{
	return RADIAN_TO_DEGREE * vector_angle_radian(u, v);
}

@ @(mcs.cu@>=
__device__ double cuda_vector_angle_degree(Vector u, Vector v)
{
	return RADIAN_TO_DEGREE * vector_angle_radian(u, v);
}


@ Function |vector_distance(u, v)| returns the distance between the
three-dimensional points represented by the vectors |u| and
|v|. Vectors |u| and |v| are left unmodified.
@<Global functions@>=
double vector_distance(Vector u, Vector v)
{
	double x, y, z;
	x = u[0] - v[0];
	y = u[1] - v[1];
	z = u[2] - v[2];
	return sqrt(x * x + y * y + z * z);
}

@ @(mcs.cu@>=
__device__ double cuda_vector_distance(Vector u, Vector v)
{
	double x, y, z;
	x = u[0] - v[0];
	y = u[1] - v[1];
	z = u[2] - v[2];
	return sqrt(x * x + y * y + z * z);
}

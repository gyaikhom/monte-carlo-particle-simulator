@q This file is part of the Monte Carlo Simulator (c) Cardiff University 2011 @>

@*1 Matrix. This section defines data-structures and functions for
representing and manipulating matrices.

In {\tt MCS}, we use a $4 \times 4$ matrix to represent affine
transformations such as translation, rotation and scaling. These
matrices are applied to vectors in homogeneous coordinates, with
$w$-component equal to 1, to determine the transformed resultant
vector.

@<Type definitions@>=
typedef double Matrix[4][4]; /* a $4 \times 4$ matrix */

@ The {\sl identity matrix} represents a null transformation. Applying
an affine transformation represented by an identify matrix leaves the
original vector unmodified. All new affine transformation matrices
must be initialised to the |IDENTITY_MATRIX|. When we wish to convert
an existing matrix into the identity matrix, we simply copy the values
stored in the immutable global variable |identity_matrix|.

@d IDENTITY_MATRIX {
   {1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0},
   {0.0, 0.0, 0.0, 1.0}
}
@<Global variables@>=
const Matrix identity_matrix = IDENTITY_MATRIX;

@ Function |matrix_print(f,m,nr,nc,t)| prints the $|nr| \times |nc|$
matrix |m| to the I/O stream pointed to by |f| using an indentation of
|t| blank tabs.
@<Global functions@>=
void matrix_print(FILE *f, Matrix m, uint8_t nr, uint8_t nc, uint8_t t)
{
        uint8_t i, j, k;
	for (i = 0; i < nr; ++i) {
                for (k = 0; k < t; ++k) fprintf(f, "\t");
                fprintf(f, "| ");
	        for (j = 0; j < nc; ++j) fprintf(f, "%8.3lf ", m[i][j]);
                fprintf(f, "|\n");
        }
}

@ Function |matrix_copy_f(u,v)| modifies the matrix |u| by copying all
of the 16 elements from matrix |v| to matrix |u|. Matrix |v| is left
unmodified.
@d matrix_copy(X, Y) memcpy((X), (Y), 16*sizeof(double))
@<Global functions@>=
void matrix_copy_f(Matrix u, Matrix v)
{
	u[0][0] = v[0][0];
	u[0][1] = v[0][1];
	u[0][2] = v[0][2];
	u[0][3] = v[0][3];
	u[1][0] = v[1][0];
	u[1][1] = v[1][1];
	u[1][2] = v[1][2];
	u[1][3] = v[1][3];
	u[2][0] = v[2][0];
	u[2][1] = v[2][1];
	u[2][2] = v[2][2];
	u[2][3] = v[2][3];
	u[3][0] = v[3][0];
	u[3][1] = v[3][1];
	u[3][2] = v[3][2];
	u[3][3] = v[3][3];
}

@ Function |matrix_multiply(a,b,c)| multiplies the $4 \times 4$ matrix
|a| to the $4 \times 4$ matrix |b|, and stores the result in matrix
|c|. Matrix multiplication is {\sl noncommutative}---the order of the
parameters are significant.@^matrix multiplication@>@^noncommutative@>
@<Global functions@>=
void matrix_multiply(Matrix a, Matrix b, Matrix c)
{
c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0] + a[0][3] * b[3][0];
c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1] + a[0][3] * b[3][1];
c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2] + a[0][3] * b[3][2];
c[0][3] = a[0][0] * b[0][3] + a[0][1] * b[1][3] + a[0][2] * b[2][3] + a[0][3] * b[3][3];
c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0] + a[1][3] * b[3][0];
c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1] + a[1][3] * b[3][1];
c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2] + a[1][3] * b[3][2];
c[1][3] = a[1][0] * b[0][3] + a[1][1] * b[1][3] + a[1][2] * b[2][3] + a[1][3] * b[3][3];
c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0] + a[2][3] * b[3][0];
c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1] + a[2][3] * b[3][1];
c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2] + a[2][3] * b[3][2];
c[2][3] = a[2][0] * b[0][3] + a[2][1] * b[1][3] + a[2][2] * b[2][3] + a[2][3] * b[3][3];
c[3][0] = a[3][0] * b[0][0] + a[3][1] * b[1][0] + a[3][2] * b[2][0] + a[3][3] * b[3][0];
c[3][1] = a[3][0] * b[0][1] + a[3][1] * b[1][1] + a[3][2] * b[2][1] + a[3][3] * b[3][1];
c[3][2] = a[3][0] * b[0][2] + a[3][1] * b[1][2] + a[3][2] * b[2][2] + a[3][3] * b[3][2];
c[3][3] = a[3][0] * b[0][3] + a[3][1] * b[1][3] + a[3][2] * b[2][3] + a[3][3] * b[3][3];
}

@ Function |matrix_determinant(m)| returns the {\sl determinant}
of a $4 \times 4$ matrix. Matrix |m| is left unmodified.
@^matrix determinant@>
@<Global functions@>=
double matrix_determinant(Matrix m)
{
   return
   m[0][3] * m[1][2] * m[2][1] * m[3][0] -
   m[0][2] * m[1][3] * m[2][1] * m[3][0] -
   m[0][3] * m[1][1] * m[2][2] * m[3][0] +
   m[0][1] * m[1][3] * m[2][2] * m[3][0] +
   m[0][2] * m[1][1] * m[2][3] * m[3][0] -
   m[0][1] * m[1][2] * m[2][3] * m[3][0] -
   m[0][3] * m[1][2] * m[2][0] * m[3][1] +
   m[0][2] * m[1][3] * m[2][0] * m[3][1] +
   m[0][3] * m[1][0] * m[2][2] * m[3][1] -
   m[0][0] * m[1][3] * m[2][2] * m[3][1] -
   m[0][2] * m[1][0] * m[2][3] * m[3][1] +
   m[0][0] * m[1][2] * m[2][3] * m[3][1] +
   m[0][3] * m[1][1] * m[2][0] * m[3][2] -
   m[0][1] * m[1][3] * m[2][0] * m[3][2] -
   m[0][3] * m[1][0] * m[2][1] * m[3][2] +
   m[0][0] * m[1][3] * m[2][1] * m[3][2] +
   m[0][1] * m[1][0] * m[2][3] * m[3][2] -
   m[0][0] * m[1][1] * m[2][3] * m[3][2] -
   m[0][2] * m[1][1] * m[2][0] * m[3][3] +
   m[0][1] * m[1][2] * m[2][0] * m[3][3] +
   m[0][2] * m[1][0] * m[2][1] * m[3][3] -
   m[0][0] * m[1][2] * m[2][1] * m[3][3] -
   m[0][1] * m[1][0] * m[2][2] * m[3][3] +
   m[0][0] * m[1][1] * m[2][2] * m[3][3];
}

@ Function |matrix_inverse(m,i)| calculates the {\sl inverse} of a $4
\times 4$ matrix and stores the result in the $4 \times 4$ matrix
|i|. Matrix |m| is left unmodified.@^matrix inverse@>
@<Global functions@>=
void matrix_inverse(Matrix m, Matrix i) {
   double det = matrix_determinant(m);
   i[0][0] = (m[1][2] * m[2][3] * m[3][1] -
             m[1][3] * m[2][2] * m[3][1] +
	     m[1][3] * m[2][1] * m[3][2] -
	     m[1][1] * m[2][3] * m[3][2] -
	     m[1][2] * m[2][1] * m[3][3] +
	     m[1][1] * m[2][2] * m[3][3]) / det;

   i[0][1] = (m[0][3] * m[2][2] * m[3][1] -
             m[0][2] * m[2][3] * m[3][1] -
	     m[0][3] * m[2][1] * m[3][2] +
	     m[0][1] * m[2][3] * m[3][2] +
	     m[0][2] * m[2][1] * m[3][3] -
	     m[0][1] * m[2][2] * m[3][3]) / det;

   i[0][2] = (m[0][2] * m[1][3] * m[3][1] -
             m[0][3] * m[1][2] * m[3][1] +
	     m[0][3] * m[1][1] * m[3][2] -
	     m[0][1] * m[1][3] * m[3][2] -
	     m[0][2] * m[1][1] * m[3][3] +
	     m[0][1] * m[1][2] * m[3][3]) / det;

   i[0][3] = (m[0][3] * m[1][2] * m[2][1] -
             m[0][2] * m[1][3] * m[2][1] -
	     m[0][3] * m[1][1] * m[2][2] +
	     m[0][1] * m[1][3] * m[2][2] +
	     m[0][2] * m[1][1] * m[2][3] -
	     m[0][1] * m[1][2] * m[2][3]) / det;

   i[1][0] = (m[1][3] * m[2][2] * m[3][0] -
             m[1][2] * m[2][3] * m[3][0] -
	     m[1][3] * m[2][0] * m[3][2] +
	     m[1][0] * m[2][3] * m[3][2] +
	     m[1][2] * m[2][0] * m[3][3] -
	     m[1][0] * m[2][2] * m[3][3]) / det;

   i[1][1] = (m[0][2] * m[2][3] * m[3][0] -
   	     m[0][3] * m[2][2] * m[3][0] +
	     m[0][3] * m[2][0] * m[3][2] -
	     m[0][0] * m[2][3] * m[3][2] -
	     m[0][2] * m[2][0] * m[3][3] +
	     m[0][0] * m[2][2] * m[3][3]) / det;

   i[1][2] = (m[0][3] * m[1][2] * m[3][0] -
             m[0][2] * m[1][3] * m[3][0] -
	     m[0][3] * m[1][0] * m[3][2] +
	     m[0][0] * m[1][3] * m[3][2] +
	     m[0][2] * m[1][0] * m[3][3] -
	     m[0][0] * m[1][2] * m[3][3]) / det;

   i[1][3] = (m[0][2] * m[1][3] * m[2][0] -
             m[0][3] * m[1][2] * m[2][0] +
	     m[0][3] * m[1][0] * m[2][2] -
	     m[0][0] * m[1][3] * m[2][2] -
	     m[0][2] * m[1][0] * m[2][3] +
	     m[0][0] * m[1][2] * m[2][3]) / det;

   i[2][0] = (m[1][1] * m[2][3] * m[3][0] -
             m[1][3] * m[2][1] * m[3][0] +
	     m[1][3] * m[2][0] * m[3][1] -
	     m[1][0] * m[2][3] * m[3][1] -
	     m[1][1] * m[2][0] * m[3][3] +
	     m[1][0] * m[2][1] * m[3][3]) / det;

   i[2][1] = (m[0][3] * m[2][1] * m[3][0] -
             m[0][1] * m[2][3] * m[3][0] -
	     m[0][3] * m[2][0] * m[3][1] +
	     m[0][0] * m[2][3] * m[3][1] +
	     m[0][1] * m[2][0] * m[3][3] -
	     m[0][0] * m[2][1] * m[3][3]) / det;

   i[2][2] = (m[0][1] * m[1][3] * m[3][0] -
             m[0][3] * m[1][1] * m[3][0] +
	     m[0][3] * m[1][0] * m[3][1] -
	     m[0][0] * m[1][3] * m[3][1] -
	     m[0][1] * m[1][0] * m[3][3] +
	     m[0][0] * m[1][1] * m[3][3]) / det;

   i[2][3] = (m[0][3] * m[1][1] * m[2][0] -
             m[0][1] * m[1][3] * m[2][0] -
	     m[0][3] * m[1][0] * m[2][1] +
	     m[0][0] * m[1][3] * m[2][1] +
	     m[0][1] * m[1][0] * m[2][3] -
	     m[0][0] * m[1][1] * m[2][3]) / det;

   i[3][0] = (m[1][2] * m[2][1] * m[3][0] -
             m[1][1] * m[2][2] * m[3][0] -
	     m[1][2] * m[2][0] * m[3][1] +
	     m[1][0] * m[2][2] * m[3][1] +
	     m[1][1] * m[2][0] * m[3][2] -
	     m[1][0] * m[2][1] * m[3][2]) / det;

   i[3][1] = (m[0][1] * m[2][2] * m[3][0] -
             m[0][2] * m[2][1] * m[3][0] +
	     m[0][2] * m[2][0] * m[3][1] -
	     m[0][0] * m[2][2] * m[3][1] -
	     m[0][1] * m[2][0] * m[3][2] +
	     m[0][0] * m[2][1] * m[3][2]) / det;

   i[3][2] = (m[0][2] * m[1][1] * m[3][0] -
             m[0][1] * m[1][2] * m[3][0] -
	     m[0][2] * m[1][0] * m[3][1] +
	     m[0][0] * m[1][2] * m[3][1] +
	     m[0][1] * m[1][0] * m[3][2] -
	     m[0][0] * m[1][1] * m[3][2]) / det;

   i[3][3] = (m[0][1] * m[1][2] * m[2][0] -
             m[0][2] * m[1][1] * m[2][0] +
	     m[0][2] * m[1][0] * m[2][1] -
	     m[0][0] * m[1][2] * m[2][1] -
	     m[0][1] * m[1][0] * m[2][2] +
	     m[0][0] * m[1][1] * m[2][2]) / det;
}

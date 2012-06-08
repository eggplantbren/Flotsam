#ifndef Matrix_h
#define Matrix_h

#include <gsl/gsl_matrix.h>

// A wrapper around GSL's matrices and vectors
class Matrix
{
	private:
		gsl_matrix* mat;

	public:
		Matrix(int i, int j)
		:mat(gsl_matrix_alloc(i, j))
		{
		}

		Matrix(const Matrix& other)
		{
			gsl_matrix_memcpy(mat, other.mat);
		}

		~Matrix()
		{
			gsl_matrix_free(mat);
		}

		double& operator () (int i, int j)
		{
			return *(gsl_matrix_ptr(mat, i, j));
		}
};

#endif


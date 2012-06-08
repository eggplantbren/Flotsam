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
			mat = gsl_matrix_alloc(other.mat->size1, other.mat->size2);
			gsl_matrix_memcpy(mat, other.mat);
		}

		~Matrix()
		{
			gsl_matrix_free(mat);
		}

		Matrix& operator = (const Matrix& other)
		{
			gsl_matrix_free(mat);
			mat = gsl_matrix_alloc(other.mat->size1, other.mat->size2);
			gsl_matrix_memcpy(mat, other.mat);
			return *this;
		}

		double& operator () (int i, int j)
		{
			return *(gsl_matrix_ptr(mat, i, j));
		}

		double operator () (int i, int j) const
		{
			return *(gsl_matrix_ptr(mat, i, j));
		}

		gsl_matrix* get_gsl_matrix()
		{
			return mat;
		}

		const gsl_matrix* get_gsl_matrix() const
		{
			return mat;
		}
};


#include <gsl/gsl_vector.h>

// A wrapper around GSL's matrices and vectors
class Vector
{
	private:
		gsl_vector* vec;

	public:
		Vector(int i)
		:vec(gsl_vector_alloc(i))
		{
		}

		Vector(const Vector& other)
		{
			vec = gsl_vector_alloc(other.vec->size);
			gsl_vector_memcpy(vec, other.vec);
		}

		~Vector()
		{
			gsl_vector_free(vec);
		}

		Vector& operator = (const Vector& other)
		{
			gsl_vector_free(vec);
			vec = gsl_vector_alloc(other.vec->size);
			gsl_vector_memcpy(vec, other.vec);
			return *this;
		}

		double& operator () (int i)
		{
			return *(gsl_vector_ptr(vec, i));
		}

		double operator () (int i) const
		{
			return *(gsl_vector_ptr(vec, i));
		}

		gsl_vector* get_gsl_vector()
		{
			return vec;
		}

		const gsl_vector* get_gsl_vector() const
		{
			return vec;
		}
};



#endif


#ifndef Matrix_h
#define Matrix_h

#include <gsl/gsl_matrix.h>

// A wrapper around GSL's matrices and vectors
class Matrix
{
	private:
		gsl_matrix* mat;

	public:
		Matrix(int i, int j);
		Matrix(const Matrix& other);
		~Matrix();
		Matrix& operator = (const Matrix& other);
		double& operator () (int i, int j);
		double operator () (int i, int j) const;
		gsl_matrix* get_gsl_matrix();
		const gsl_matrix* get_gsl_matrix() const;
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
			if(vec->size != other.vec->size)
			{
				gsl_vector_free(vec);
				vec = gsl_vector_alloc(other.vec->size);
			}
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


#ifndef Flotsam_Matrix_h
#define Flotsam_Matrix_h

#include <gsl/gsl_matrix.h>

// A wrapper around GSL's matrices
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

#endif


#include "Matrix.h"

Matrix::Matrix(int i, int j)
:mat(gsl_matrix_alloc(i, j))
{
}

Matrix::Matrix(const Matrix& other)
{
	mat = gsl_matrix_alloc(other.mat->size1, other.mat->size2);
	gsl_matrix_memcpy(mat, other.mat);
}

Matrix::~Matrix()
{
	gsl_matrix_free(mat);
}

Matrix& Matrix::operator = (const Matrix& other)
{
	if(mat->size1 != other.mat->size1 || mat->size2 != other.mat->size2)
	{
		gsl_matrix_free(mat);
		mat = gsl_matrix_alloc(other.mat->size1, other.mat->size2);
	}
	gsl_matrix_memcpy(mat, other.mat);
	return *this;
}

double& Matrix::operator () (int i, int j)
{
	return *(gsl_matrix_ptr(mat, i, j));
}

double Matrix::operator () (int i, int j) const
{
	return *(gsl_matrix_ptr(mat, i, j));
}

gsl_matrix* Matrix::get_gsl_matrix()
{
	return mat;
}

const gsl_matrix* Matrix::get_gsl_matrix() const
{
	return mat;
}


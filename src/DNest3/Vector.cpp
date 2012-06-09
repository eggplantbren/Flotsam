#include "Vector.h"

Vector::Vector(int i)
:vec(gsl_vector_alloc(i))
{
}

Vector::Vector(const Vector& other)
{
	vec = gsl_vector_alloc(other.vec->size);
	gsl_vector_memcpy(vec, other.vec);
}

Vector::~Vector()
{
	gsl_vector_free(vec);
}

Vector& Vector::operator = (const Vector& other)
{
	if(vec->size != other.vec->size)
	{
		gsl_vector_free(vec);
		vec = gsl_vector_alloc(other.vec->size);
	}
	gsl_vector_memcpy(vec, other.vec);
	return *this;
}

double& Vector::operator () (int i)
{
	return *(gsl_vector_ptr(vec, i));
}

double Vector::operator () (int i) const
{
	return *(gsl_vector_ptr(vec, i));
}

gsl_vector* Vector::get_gsl_vector()
{
	return vec;
}

const gsl_vector* Vector::get_gsl_vector() const
{
	return vec;
}


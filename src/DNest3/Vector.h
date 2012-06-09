#ifndef Flotsam_Vector_h
#define Flotsam_Vector_h

#include <gsl/gsl_vector.h>

// A wrapper around GSL's vectors
class Vector
{
	private:
		gsl_vector* vec;

	public:
		Vector(int i);
		Vector(const Vector& other);
		~Vector();
		Vector& operator = (const Vector& other);
		double& operator () (int i);
		double operator () (int i) const;
		gsl_vector* get_gsl_vector();
		const gsl_vector* get_gsl_vector() const;
};

#endif


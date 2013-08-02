//
//  vectorc.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/17/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "vectorc.h"
#include "complex_util.h"
#include "gsl/gsl_complex_math.h"

//Need to fix ostream operator
std::ostream& operator<<(std::ostream & os, const vectorc::ComplexProxy& cp)
{
//	os<<"("<<GSL_REAL(gsl_vector_complex_get(cp.theVector.v_m, cp.ComplexIndex))<<","<<GSL_IMAG(gsl_vector_complex_get(cp.theVector.v_m, cp.ComplexIndex))<<")";
	return os;	
}

vectorc::ComplexProxy::ComplexProxy(vectorc& vec, int index) : theVector(vec), ComplexIndex(index)
{}

vectorc::ComplexProxy::operator gsl_complex() const
{
	return gsl_vector_complex_get(theVector.v_m, ComplexIndex);
}

vectorc::ComplexProxy& vectorc::ComplexProxy::operator=(const ComplexProxy& rhs)
{
	gsl_vector_complex_set(theVector.v_m, ComplexIndex, gsl_vector_complex_get(rhs.theVector.v_m, ComplexIndex));
	return *this;
}

vectorc::ComplexProxy& vectorc::ComplexProxy::operator=(gsl_complex c)
{

	gsl_vector_complex_set(theVector.v_m, ComplexIndex, c);
	return *this;
}

vectorc::vectorc(int size) : size_m(size), v_m(gsl_vector_complex_calloc(size))
{}

vectorc::~vectorc()
{
	gsl_vector_complex_free(v_m);
};

vectorc::vectorc(const double values[], const int size) : size_m(size), v_m(gsl_vector_complex_alloc(size))

{
	gsl_complex z;
	for (int i=0; i<size_m; i++)
	{
		GSL_SET_COMPLEX(&z, values[2*i], values[2*i+1]);
		gsl_vector_complex_set(v_m, i, z);
	}	
}


vectorc::vectorc(const vectorc& vec) : size_m(vec.size_m),  v_m(gsl_vector_complex_alloc(vec.size_m))
{
	for (int i=0; i<size_m; i++)
		gsl_vector_complex_set(v_m, i, gsl_vector_complex_get(vec.v_m, i));
}

vectorc vectorc::operator=(const vectorc& vec)
{
	if (this==&vec)
		return *this;
	gsl_vector_complex_free (v_m);
	size_m=vec.size_m;
	v_m=gsl_vector_complex_alloc(size_m);
	for (int i=0; i<size_m; i++)
		gsl_vector_complex_set(v_m, i, gsl_vector_complex_get(vec.v_m, i));
	return *this;
}

bool vectorc::operator==(const vectorc& vec) const
{
	bool equal=true;
	for (int i=0; i<size_m; i++)
	{
		if (GSL_COMPLEX_EQ(gsl_vector_complex_get(v_m, i), gsl_vector_complex_get(vec.v_m, i)))
		{
			equal=false;
			break;
		}
	}
	return equal;
}

const vectorc::ComplexProxy vectorc::operator[](int index) const
{
	return ComplexProxy(const_cast<vectorc&>(*this), index);
	
}

vectorc::ComplexProxy vectorc::operator[](int index)
{
	return ComplexProxy(*this, index);
}

int vectorc::get_size() const
{
	return size_m;
}

void vectorc::get_values_copy(double *vector) const
{
	for (int i=0; i<size_m; i++)
		{
			vector[2*i]=GSL_REAL(gsl_vector_complex_get(v_m, i));
			vector[2*i+1]=GSL_IMAG(gsl_vector_complex_get(v_m, i));
		}
}

vectorc vectorc::operator+(const vectorc & vec) const
{
	try
	{
		if (size_m!=vec.size_m)
			throw("\nCan't add vectors of different sizes.");
	}
	catch(char* str)
	{
		std::cout << "Exception raised: " << str << '\n';
	}
	vectorc tot(*this);
	gsl_vector_complex_add(tot.v_m, vec.v_m);
	return tot;
}

void vectorc::set_values(double *vals) const
{
	gsl_complex z;
	for (int i=0; i<size_m; i++)
		{
			GSL_SET_COMPLEX(&z, vals[2*i], vals[2*i+1]);
			gsl_vector_complex_set(v_m, i, z);
		}
}

vectorc vectorc::operator-(const vectorc& vec) const
{
	try
	{
		if (size_m!=vec.size_m)
			throw("\nCan't add vectors of different sizes.");
	}
	catch(char* str)
	{
		std::cout << "Exception raised: " << str << '\n';
	}
	vectorc tot(*this);
	gsl_vector_complex_sub(tot.v_m, vec.v_m);
	return tot;
}

vectorc vectorc::operator-() const
{
	vectorc tot(size_m);
	gsl_vector_complex_sub(tot.v_m, v_m);
	return tot;	
}

vectorc vectorc::operator*(const gsl_complex n) const
{
	vectorc result(size_m);
	for (int i=0; i<size_m; i++)
		gsl_vector_complex_set(result.v_m, i, gsl_complex_mul(n,gsl_vector_complex_get(v_m,i)));
	return result;
}


gsl_complex vectorc::operator*(const vectorc& vec) const
{
	gsl_complex result;
	for (int i=0; i<size_m; i++)
		gsl_complex_add(result, gsl_complex_mul(gsl_vector_complex_get(v_m, i), gsl_vector_complex_get(vec.v_m, i)));
	return result;
}

vectorc operator*(const gsl_complex n, const vectorc & vec)
{
	vectorc result(vec.size_m);
	for (int i=0; i<vec.size_m; i++)
		gsl_vector_complex_set(result.v_m, i, gsl_complex_mul(n,gsl_vector_complex_get(vec.v_m,i)));
	return result;	
}

std::ostream& operator<<(std::ostream &os, const vectorc& vec)
{
	os<<"[";
	for (int i=0; i<vec.size_m; i++)
	{
		os<<vec[i];
		if (i<vec.size_m-1)
			os<<",";
	}
	os<<"]";
	return os;
}


void vectorc::clear()
{
	gsl_vector_complex_set_zero(v_m);
	
}

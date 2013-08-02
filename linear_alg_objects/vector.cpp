//
//  vector.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 6/29/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "vector.h"
#include "Accelerate/Accelerate.h"

namespace linalg
{

vector::DoubProxy::DoubProxy(vector& vec, int index) : theVector(vec), doubIndex(index)
{}
	
vector::DoubProxy::operator double() const
{
	return gsl_vector_get(theVector.v_m, doubIndex);
}
	
vector::DoubProxy& vector::DoubProxy::operator=(const DoubProxy& rhs)
{
	gsl_vector_set(theVector.v_m, doubIndex, gsl_vector_get(rhs.theVector.v_m, doubIndex));
	return *this;
}
	
vector::DoubProxy& vector::DoubProxy::operator=(double d)
{
		gsl_vector_set(theVector.v_m, doubIndex, d);
		return *this;
}
		
vector::vector(int size) : size_m(size), v_m(gsl_vector_calloc(size))
{}

vector::~vector()
{
	gsl_vector_free(v_m);
}
	
vector::vector(const double values[], const int size) : v_m(gsl_vector_alloc(size)), size_m(size)
	
{
	for (int i=0; i<size_m; i++)
		gsl_vector_set(v_m, i, values[i]);
}
	
	
vector::vector(const vector& vec) : size_m(vec.size_m),  v_m(gsl_vector_alloc(vec.size_m))
{
	gsl_vector_memcpy(v_m, vec.v_m);
}
	
vector vector::operator=(const vector & vec)
{
	if (this==&vec)
		return *this;
	gsl_vector_free(v_m);
	size_m=vec.size_m;
	v_m=gsl_vector_alloc(size_m);
	gsl_vector_memcpy(v_m, vec.v_m);
		return *this;
}
	
bool vector::operator==(const vector& vec) const
{
	bool is_equal=true;
	for (int i=0; i<size_m; i++)
	{
		if (gsl_vector_get(v_m, i) != gsl_vector_get(vec.v_m, i))
		{
			is_equal=false;
			break;
		}
		
	}
	return is_equal;		
}

const vector::DoubProxy vector::operator[](int index) const
{
	return DoubProxy(const_cast<vector&>(*this), index);
		
}
	
vector::DoubProxy vector::operator[](int index)
{
	return DoubProxy(*this, index);	    
}
	
int vector::get_size() const
{
	return size_m;
};

const gsl_vector* vector::get_gsl_vector() const
{
	return v_m;
		
}

void vector::get_values_copy(double *vector) const
{
	for (int i=0; i<size_m; i++)
		vector[i]=gsl_vector_get(v_m, i);
}

vector vector::operator+(const vector & vec) const
{
	try
	{
		if (size_m!=vec.size_m)
			throw("\nCan't add vectors of different sizes.");
	}
	catch( char * str ) 
	{
		std::cout << "Exception raised: " << str << '\n';
	}
	
	vector tot(vec);
	gsl_vector_add(tot.v_m, v_m);
	return tot;
}


void vector::set_values(double *vals) const
{
	for (int i=0; i<size_m; i++)
		gsl_vector_set(v_m, i, vals[i]);
}

void vector::set_values(gsl_vector *vec) const
{
	for(int i=0; i<size_m; i++)
		gsl_vector_set(v_m, i, gsl_vector_get(vec,i));
		
}
	
vector vector::operator-(const vector & vec) const
{
	try
	{
		if (size_m!=vec.size_m)
			throw("\nCan't add vectors of different sizes.");
	}
	catch( char * str )
	{
		std::cout << "Exception raised: " << str << '\n';
	}
	vector tot(vec);
	gsl_vector_sub(tot.v_m, v_m);
	return tot;
}
	
vector vector::operator-() const
{
	vector tot(size_m);
	gsl_vector_sub(tot.v_m, v_m);
	return tot;
}
	
vector vector::operator*(const double n) const
{
	vector result(*this);
	for (int i=0; i<size_m; i++)
		gsl_vector_set(result.v_m, n*gsl_vector_get(v_m, i), i);
	return result;
}
	
	
double vector::operator*(const vector& vec) const
{
	double result=0;
	for (int i=0; i<size_m; i++)
		result+=gsl_vector_get(v_m, i)*gsl_vector_get(vec.v_m,i);
	return result;
}

vector operator*(const double n, const vector & vec)
{
	vector result(vec.size_m);
	for (int i=0; i<vec.size_m; i++)
		gsl_vector_set(result.v_m, i, n*gsl_vector_get(vec.v_m, i));
	return result;	
}
	
std::ostream& operator<<(std::ostream &os, const vector& vec)
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
	
vector vector::e(int ei)
{
	gsl_vector_set_basis(v_m, ei);
	return *this;
}
	
void vector::clear()
{
	gsl_vector_set_zero(v_m);
		
}
}
//
//  vector.h
//  linear_alg_objects
//
//  Created by Jon Lederman on 6/29/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef linear_alg_objects_vector_h
#define linear_alg_objects_vector_h
#include <iostream>
#include "gsl/gsl_vector.h"

namespace linalg
{
class vector
{
protected:
	gsl_vector* v_m;
    int size_m;
	
public:
    class DoubProxy
    {
    public:
        DoubProxy(vector& vec, int index);          //creation
        DoubProxy& operator=(const DoubProxy& rhs); //lvalue uses
        DoubProxy& operator=(double c);             
        operator double() const;                    //rvalue use
		
	private:
        vector& theVector;                          //vector this proxy pertains to
        int doubIndex;                              //char within vector
    };
    
    vector(int size);
    vector(const double values[], const int size);
    vector(const vector& vec);  //copy constructor
    ~vector();
    int get_size() const;
	const gsl_vector* get_gsl_vector() const;
	void get_values_copy(double *vector) const;
	void set_values(double *vals) const;
	void set_values(gsl_vector *vec) const;
    vector operator=(const vector & vec);  			//assignment operator
    const DoubProxy operator[](int index) const; 	//for const Vectors
    DoubProxy operator[](int index);             	//for non-const Vectors
    vector operator+(const vector & vec) const;
    vector operator-(const vector & vec) const;
    vector operator-() const;
    vector operator*(const double n) const;
    double operator*(const vector& vec) const;
	bool operator==(const vector& vec) const;
    friend vector operator*(const double n, const vector& vec);
    friend std::ostream & operator<<(std::ostream & os, const vector& vec);
	vector e(int ei);
	void clear();
    friend class DoubProxy;
};
};
#endif

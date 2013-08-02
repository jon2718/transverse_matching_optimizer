//
//  vectorc.h
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/17/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef linear_alg_objects_vectorc_h
#define linear_alg_objects_vectorc_h
#include "gsl/gsl_eigen.h"

class vectorc
{
protected:
	gsl_vector_complex* v_m;
    int size_m;
	
public:
    class ComplexProxy
    {
    public:
        ComplexProxy(vectorc& vec, int index);          //creation
        ComplexProxy& operator=(const ComplexProxy& rhs); //lvalue uses
        ComplexProxy& operator=(gsl_complex c);
        operator gsl_complex() const;                    //rvalue use
		friend std::ostream & operator<<(std::ostream & os, const ComplexProxy& cp);
		
	private:
        vectorc& theVector;                          //vector this proxy pertains to
        int ComplexIndex;                              //char within vector
    };
    
    vectorc(int size);
    vectorc(const double values[], const int size);
    vectorc(const vectorc& vec);  //copy constructor
    ~vectorc();
    int get_size() const;
	void get_values_copy(double *vector) const;
	void set_values(double *vals) const;
    vectorc operator=(const vectorc & vec);  			//assignment operator
    const ComplexProxy operator[](int index) const; 	//for const Vectors
    ComplexProxy operator[](int index);             	//for non-const Vectors
    vectorc operator+(const vectorc& vec) const;
    vectorc operator-(const vectorc& vec) const;
    vectorc operator-() const;
    vectorc operator*(const gsl_complex n) const;
    gsl_complex operator*(const vectorc& vec) const;
	bool operator==(const vectorc& vec) const;
    friend vectorc operator*(const gsl_complex n, const vectorc& vec);
    friend std::ostream & operator<<(std::ostream & os, const vectorc& vec);
	void clear();
    friend class ComplexProxy;
};



#endif

//
//  eigen.h
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef linear_alg_objects_eigen_h
#define linear_alg_objects_eigen_h
#include "vector.h"
#include "matrix.h"
#include "matrixc.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

class eigen
{
public:
	
	eigen(const matrix& mat);
	eigen(const matrixc& mat);
	eigen(int size);
	eigen(const eigen& eig);
	~eigen();
	eigen& operator=(const eigen& eig);
	const matrix& get_eigenvectors() const;
	const vector& get_eigenvalues() const;
	void set_eigenvalues(gsl_vector &vec);
	void set_eigenvectors(gsl_matrix &mat);
	
	
private:
	matrix eigenvectors;
	vector eigenvalues;
	
};


#endif

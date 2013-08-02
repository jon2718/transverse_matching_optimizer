//
//  matrix_exceptions.h
//  linear_alg_objects
//
//  Created by Jon Lederman on 8/23/12.
//
//

#ifndef __linear_alg_objects__matrix_exceptions__
#define __linear_alg_objects__matrix_exceptions__

#include <iostream>
#include <stdexcept>

class NotSymmetricMatrixException : public std::runtime_error
{
public:
	NotSymmetricMatrixException()
	: std::runtime_error("Matrix is not symmetric."){};
	
};


#endif /* defined(__linear_alg_objects__matrix_exceptions__) */

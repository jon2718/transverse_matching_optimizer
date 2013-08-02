//
//  eigen_unittest.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 8/22/12.
//
//


#include "gtest/gtest.h"
#include "matrix.h"
#include "eigen.h"

// Tests eigenfunctions
TEST(LinearAlgebraObjects, EigenTest)
{
	// This test is named "Function", and belongs to the "MatchingTest"
	// test case.
	eigen ab(4);
	double first[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	double second[16]={1,2,3,4,2,5,8,14,3,8,7,9,4,14,9,12};
	double third[16]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
	
	matrix a(4,4,first);
	matrix b(4,4, second);
	matrix c(4,4);
	c=a*b;
	std::cout<<c;
	b.eig_symm(ab);
//	a.eig_symm(ab);
	std::cout<<ab.get_eigenvalues();
	std::cout<<ab.get_eigenvectors();
	
		
	
	
}
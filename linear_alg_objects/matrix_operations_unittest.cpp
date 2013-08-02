//
//  matrix_operations_unittest.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/2/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>

#include "gtest/gtest.h"
#include "matrix.h"

// Tests creation of functions
TEST(LinearAlgebraObjects, MatrixMultiplication) {
	// This test is named "Function", and belongs to the "MatchingTest"
	// test case.
	
	double first[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	double second[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	
	matrix a(4,4,first);
	matrix b(4,4,second);
	matrix c(4,4);
	c=a*b;
	
	
	EXPECT_EQ((double)c[0][0], 90)<<"[0][0] component is incorrect.";
	EXPECT_EQ((double)c[3][3], 600)<<"[3][3] component is incorrect.";
	EXPECT_EQ((double)a[3][3], 16)<<"[3][3] component is incorrect.";
	EXPECT_EQ((double)a[2][2], 11)<<"[2][2] component is incorrect.";
	EXPECT_EQ((double)a[1][2], 7)<<"[1][2] component is incorrect.";
	
	a[1][3]=12;
	EXPECT_EQ((double)a[1][3], 12)<<"Set of [1][3] component is incorrect.";
	EXPECT_EQ((double)a[0][3], 4)<<"Set of [0][3] component is incorrect.";
	
	double left[9]={1,2,3,
					4,5,6,
					7,8,9};
	
	double right[9]={1,0,0,
					0,1,0,
					0,0,1};
	
	matrix d(3,3,left);
	matrix e(3,3,right);
	matrix f(3,3);
	f=d*e;
	
	EXPECT_EQ((double)f[2][2], 9)<<"Multiplication by identity failed.";
	EXPECT_TRUE(f==d)<<"Identity multiplication failed.";
	std::cout<<f;
	
	
	
}

TEST(LinearAlgebraObjects, MatrixVectorMultiplication) {
	// This test is named "Function", and belongs to the "MatchingTest"
	// test case.
	
	double first[4]={1,2,3,4};
	
	double second[2]={1,1};
	
	matrix a(2,2,first);
	vector b(second, 2);
	vector c(2);
	c=a*b;
	
	
	EXPECT_EQ((double)c[0], 3)<<"0th component is incorrect.";
	EXPECT_EQ((double)c[1], 7)<<"first component is incorrect.";
	
	//std::cout<<c;
	
	
	
	
}


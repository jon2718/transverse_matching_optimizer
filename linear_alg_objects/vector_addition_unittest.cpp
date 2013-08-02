//
//  vector_addition_unittest.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 6/29/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "gtest/gtest.h"
#include "vector.h"

// Tests creation of functions
TEST(LinearAlgebraObjects, Vector_Addition) {
	// This test is named "Function", and belongs to the "MatchingTest"
	// test case.
	
	double first[5]={1,2,3,4,5};
	double second[5]={6,7,8,9,10};
	
	linalg::vector a(first, 5);
	linalg::vector b(second, 5);
	linalg::vector c(5);
	
	c=a+b;
	EXPECT_EQ(c[0], 7)<<"First component is incorrect.";
	EXPECT_EQ(c[1], 9)<<"Second component is incorrect.";
	EXPECT_EQ(c[2], 11)<<"Third component is incorrect.";
	EXPECT_EQ(c[3], 13)<<"Fourth component is incorrect.";
	EXPECT_EQ(c[4], 15)<<"Fifth component is incorrect.";




}
//
//  setting_getting_unittest.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/2/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "gtest/gtest.h"
#include "matrix.h"

// Tests creation of functions
TEST(LinearAlgebraObjects, Setting_Getting) {
	// This test is named "Function", and belongs to the "MatchingTest"
	// test case.
	
	double first[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	//double second[5]={6,7,8,9,10};
	
	
	matrix a(4,4,first);
	
	EXPECT_EQ((double)a[0][0], 1)<<"First component is incorrect.";
	EXPECT_EQ((double)a[1][1], 6)<<"[1][1] component is incorrect.";
//	EXPECT_EQ(a[3][3], 16)<<"[3][3] component is incorrect.";
//	EXPECT_EQ(a[2][2], 11)<<"[2][2] component is incorrect.";
//	EXPECT_EQ(a[1][2], 7)<<"[1][2] component is incorrect.";
	
	a[1][3]=12;
	EXPECT_EQ((double)a[1][3], 12)<<"Set of [1][3] component is incorrect.";

	
	
	
	
}
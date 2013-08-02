//
//  objectivefunctional.h
//  Metrics
//
//  Created by Jon Lederman on 11/11/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Metrics__objectivefunctional__
#define __Metrics__objectivefunctional__

#include <iostream>

#include "objective.h"
#include "interpolate_base.h"
#include "functional.h"

//Treat vector as a function (interpolator_base) and evaluate that instead.

class objectivefunctional : public objective
{
	
	
private:
	interpolator_base* fn_m;			//to hold the function
	functional* fnl_m;					//to hold the functional
	
public:
	objectivefunctional(functional *fnl, interpolator_base* fn);
	virtual double operator()(vector &v);
	
};



#endif /* defined(__Metrics__objectivefunctional__) */

//
//  functional.h
//  Metrics
//
//  Created by Jon Lederman on 11/11/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Metrics__functional__
#define __Metrics__functional__

#include <iostream>
#include "interpolate_base.h"
#include "matrix.h"

class functional
{
	
public:
	functional(){};
	virtual ~functional() {};
	virtual double operator()(interpolator_base& fn)=0;
	virtual linalg::vector deriv1(interpolator_base* B0)=0;
    virtual matrix deriv2(interpolator_base* B0)=0;
	double u(const linalg::vector& partials, const interpolator_base* f0, double h);
	
};


#endif /* defined(__Metrics__functional__) */

//
//  TInterpPlot2D.h
//  Root Plotting Extensions
//
//  Created by Jon Lederman on 9/18/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Root_Plotting_Extensions__TInterpPlot2D__
#define __Root_Plotting_Extensions__TInterpPlot2D__

#include <iostream>
#include "TObject.h"
#include "TInterpPlot2D.h"
#include "TF1.h"
#include "interpolate_base.h"

class TInterpPlot2D : public TF1
{
private:
	interpolator_base* intp_t;
	
public:
	TInterpPlot2D(char name[], interpolator_base* intp, double xmin, double xmax);
	void sety(const linalg::vector& y, const int num_points, const int start=0);
	double operator()(double* x, double* p);
	
	
};


#endif /* defined(__Root_Plotting_Extensions__TInterpPlot2D__) */

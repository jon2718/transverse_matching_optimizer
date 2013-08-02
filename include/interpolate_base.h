//
//  interpolate_base.h
//  Interpolators
//
//  Created by Jon Lederman on 8/28/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//
/* Interpolation abstract base type for different interpolation types including:
	linear, rational, cubic spline, basis spline, etc.
vector x stores the interpolation points along x axis.
vector y stores the interpolation points along the y axis.  
In the case of basis functions such as basis splines, vectors x and y
are not used and the breakpoints and coefficients are stored in member variables.
There's probably a more elegant solution to this class structure.
 */

#ifndef Interpolators_interpolate_base_h
#define Interpolators_interpolate_base_h
#include "vector.h"
#include "gsl/gsl_vector.h"
#include "clonable.h"


class interpolator_base : public clonable
{
	
protected:
	gsl_vector* x_m;
	gsl_vector* y_m;
	int num_points_x_m;
	int num_points_y_m;

public:
	interpolator_base(double x[], double y[], int num_points_x, int num_points_y);
	interpolator_base(int num_points_x, int num_points_y);
	virtual ~interpolator_base(){};
	virtual double operator()(double x)=0;
	virtual int get_num_points_x();
	virtual int get_num_points_y();
	virtual void setx(const linalg::vector& x, const int num_points, const int start=0)=0;
	virtual void sety(const linalg::vector& y, const int num_points, const int start=0)=0;
	virtual void set_bc(double left, double right)=0;
	
};


#endif

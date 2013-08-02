//
//  bspline_gsl.h
//  Interpolators
//
//  Created by Jon Lederman on 8/31/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//
/*B-splines are piecewise (k-1)st order polynomials (they have k coefficients). If the breakpoints are at z_1 ... z_m, then the usual structure is to have t_1 = ... = t_k = z_1, t_i = z_{i-k+1} for i=(k+1)...(m+k-2), and t_{m+k-1} = ... = t_{m+2k-2} = z_m. The number of basis functions is m+k-2, which is n in deBoor’s notation. You can see that my counting is right by taking the trivial cases for k=1 (step bases) and k=2 (triangle bases). m doesn’t appear in deBoor but I’ll use it anyhow.

There is a mapping between deBoor’s numbers and Wikipedia’s numbers:

n(Wiki) = k(deBoor)-1
m(Wiki) = n(deBoor)+k(deBoor)

Thus n(deBoor) = m(Wiki)-k(deBoor) = m(Wiki)-n(Wiki)-1, which agrees with the Wikipedia text just under “Uniform B-Spline.”

The counting only applies in the simpistic case with maximal nontrivial continuity at the knots and no continuity at the ends.  However, it is generally true that if you have n+k knots, you have n basis functions (deBoor’s notation). This follows simply from the recursive definition (or from any of the other definitions).
 */


#ifndef __Interpolators__bspline_gsl__
#define __Interpolators__bspline_gsl__

#include <iostream>
#include "interpolate_base.h"
#include "gsl/gsl_bspline.h"
#include "vector.h"

/*In this derived class from interpolator_base:
 	x_m is a vector of breakpoints
 	y_m is a vector of the coefficients
*/

class bspline_gsl : public interpolator_base
{
	
private:
	int num_basis_fns_m;
	int k_m;
	int nbreak_m;
	gsl_vector* breakpoints_m;
	gsl_vector* coefficients_m;
	gsl_bspline_workspace *ws;
	gsl_vector* b;  //For evaluation of bspline

public:
	bspline_gsl(const int k, const int nbreak, const double breakpoints[]);
	bspline_gsl(const int k, const int nbreak, int a, int b);  //assumes uniform breakpoints
	bspline_gsl(const bspline_gsl& bsp);
	virtual bspline_gsl* clone() const
	{return new bspline_gsl( *this ); }
	~bspline_gsl();
	bspline_gsl& operator=(const bspline_gsl& bsp);
	virtual void setx(const linalg::vector& x, const int num_points, const int start=0);
	virtual void sety(const linalg::vector& y, const int num_points, const int start=0);
	virtual void set_bc(double left, double right);
	virtual double operator()(double x);
	int get_num_basis_fns();
	int get_order();
	int get_num_breakpoints();
	void print_knots();
	virtual int get_num_points();
};

#endif /* defined(__Interpolators__bspline_gsl__) */

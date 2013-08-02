//
//  symp1d.h
//  Phase Space Propagation
//
//  Created by Jon Lederman on 10/25/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Phase_Space_Propagation__symp1d__
#define __Phase_Space_Propagation__symp1d__

#include <iostream>
#include "integrator.h"
#include "interpolate_base.h"
#include "vector.h"
#include "sep_hamiltonian.h"


class symp1d : public integrator
{
private:
    int order_m;
    int k_m;
	enum symp_order{first=1, second=2, fourth=4, sixth=6, eighth=8};
	double *c, *d, *dz;
	static double c1[1], c2[2], c4[4], c6[8], c8[16], d1[1], d2[2], d4[4], d6[8], d8[16];
	static double dz1_dkd[1],dz2_dkd[2],dz4_dkd[4], dz6_dkd[8], dz8_dkd[16];
	static int k1, k2, k4, k6, k8;
	
public:
    symp1d(const int dimensions, const int order);
	virtual void integrate(linalg::vector& state, const sep_hamiltonian& h, const double s0, const double s1, const double step);
	static void init_coeffs();
	static const void calc_coeffs(double *c, double *d, double *w, int m, int k);
};

#endif /* defined(__Phase_Space_Propagation__symp1d__) */


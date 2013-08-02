//
//  match_functional_multi.h
//  Metrics
//
//  Created by Jon Lederman on 11/13/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Metrics__match_functional_multi__
#define __Metrics__match_functional_multi__

#include <iostream>
#include <vector>
#include "integrator.h"
#include "map.h"
#include "functional.h"
#include "sep_hamiltonian.h"
#include "interpolate_base.h"



class match_functional_multi : public functional
{
private:
	integrator* intg_m;
	std::vector<sep_hamiltonian> sep_hamiltonians_m;
	std::vector<matrix> A_L_s0s_m;
	std::vector<matrix> A_L_s1s_m;
	map mp_m;
	matrix mismatch_m;
	int cpe_num_m;
	
	
public:
	match_functional_multi(integrator* intg, double s0, double s1, double step, double cpe_start, int cpe_num, double cpe_step, double B0S0, double B0S1, double qx, double pxp0, interpolator_base *B0_init);
	double operator()(interpolator_base& B0);
	virtual linalg::vector deriv1();

	
};


#endif /* defined(__Metrics__match_functional_multi__) */

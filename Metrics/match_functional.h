//
//  match_functional.h
//  Metrics
//
//  Created by Jon Lederman on 11/11/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Metrics__match_functional__
#define __Metrics__match_functional__

#include <iostream>
#include "integrator.h"
#include "map.h"
#include "functional.h"
#include "sep_hamiltonian.h"
#include "matrix.h"



class match_functional : public functional
{
private:
	integrator* intg_m;
	sep_hamiltonian h_m;
	map mp_m;
	matrix mismatch_m;
	matrix A_L_s0_m;
	matrix A_L_s1_m;
	matrix init_m;
	

public:
	match_functional(integrator* intg, sep_hamiltonian h, double s0, double s1, double step, double B0S0, double B0S1);
	virtual ~match_functional(){};
	double operator()(interpolator_base& B0);
	virtual linalg::vector deriv1(interpolator_base* B0);
	virtual matrix deriv2(interpolator_base* B0);
	
};


#endif /* defined(__Metrics__match_functional__) */

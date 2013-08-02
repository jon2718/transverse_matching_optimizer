//
//  match_functional.cpp
//  Metrics
//
//  Created by Jon Lederman on 11/11/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include "match_functional.h"
#include "physical_constants.h"
#include "math.h"

match_functional::match_functional(integrator* intg, sep_hamiltonian h, double s0, double s1, double step, double B0S0, double B0S1) : intg_m(intg), h_m(h), mp_m(map(2, intg_m, s0, s1, step)), A_L_s0_m(matrix(2,2)), A_L_s1_m(matrix(2,2)), mismatch_m(matrix(2,2)), init_m(matrix(2,2))
{
	
	A_L_s0_m[0][0]=pow((2*h_m.get_p0())/(physical_constants::e_k*abs(B0S0)),0.5);
	A_L_s0_m[1][1]=pow((physical_constants::e_k*abs(B0S0))/(2*h_m.get_p0()),0.5);
	A_L_s0_m[0][1]=0;
	A_L_s0_m[1][0]=0;
	
	A_L_s1_m[0][0]=pow((2*h_m.get_p0())/(physical_constants::e_k*abs(B0S1)),0.5);
	A_L_s1_m[1][1]=pow((physical_constants::e_k*abs(B0S1))/(2*h_m.get_p0()),0.5);
	A_L_s1_m[0][1]=0;
	A_L_s1_m[1][0]=0;
	
	
}

//A_L^-1*M*A_L=R
//A_L maps normal to real coordinates - check this!!!
double match_functional::operator()(interpolator_base& B0)
{
	h_m.set_B0(&B0);
	mp_m.gen(h_m);
	
	mismatch_m=A_L_s1_m.invert()*mp_m.get_matrix()*A_L_s0_m;
	
	double a=mismatch_m[0][0]+mismatch_m[1][1];
	double b=mismatch_m[0][1]-mismatch_m[1][0];
	
	double a_prime=mismatch_m[0][0]-mismatch_m[1][1];
	double b_prime=mismatch_m[0][1]+mismatch_m[1][0];
	
	double coshr=0.5*pow(a*a+b*b,0.5);
	double sinhr=0.5*pow(a_prime*a_prime+b_prime*b_prime,0.5);
	return sinhr;

}


linalg::vector match_functional::deriv1(interpolator_base* B0)
{
    //To be completed
}

matrix match_functional::deriv2(interpolator_base* B0)
{
	//To be completed
}

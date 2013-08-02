//
//  match_functional_multi.cpp
//  Metrics
//
//  Created by Jon Lederman on 11/13/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include "math.h"
#include "match_functional_multi.h"
#include "match_functional_multi.h"
#include "sep_hamiltonian.h"
#include "physical_constants.h"


match_functional_multi::match_functional_multi(integrator* intg, double s0, double s1, double step, double cpe_start, int cpe_num, double cpe_step, double B0S0, double B0S1, double qx, double pxp0, interpolator_base *B0_init) : intg_m(intg), mp_m(map(2, intg_m, s0, s1, step)), mismatch_m(matrix(2,2)), cpe_num_m(cpe_num)
{
	for (int i=0; i<cpe_num; i++)
	{
		sep_hamiltonians_m.push_back(*(new sep_hamiltonian(cpe_start+i*cpe_step, qx, pxp0, B0_init)));
		
		A_L_s0s_m.push_back(*(new matrix(2,2)));
		A_L_s1s_m.push_back(*(new matrix(2,2)));
		
		(A_L_s0s_m[i])[0][0]=pow((2*sep_hamiltonians_m[i].get_p0())/(physical_constants::e_k*abs(B0S0)),0.5);
		(A_L_s0s_m[i])[1][1]=pow((physical_constants::e_k*abs(B0S0))/(2*sep_hamiltonians_m[i].get_p0()),0.5);
		(A_L_s0s_m[i])[0][1]=0;
		(A_L_s0s_m[i])[1][0]=0;
		
		(A_L_s1s_m[i])[0][0]=pow((2*sep_hamiltonians_m[i].get_p0())/(physical_constants::e_k*abs(B0S1)),0.5);
		(A_L_s1s_m[i])[1][1]=pow((physical_constants::e_k*abs(B0S1))/(2*sep_hamiltonians_m[i].get_p0()),0.5);
		(A_L_s1s_m)[i][0][1]=0;
		(A_L_s1s_m[i])[1][0]=0;
		
	}
	
}
 
double match_functional_multi::operator()(interpolator_base& B0)
{
	
	double a, b, a_prime, b_prime, coshr, sinhr;
	double p0_sum=0;
	for (int i=0; i<cpe_num_m; i++)
	{
		sep_hamiltonians_m[i].set_B0(&B0);
		mp_m.gen(sep_hamiltonians_m[i]);
		mismatch_m=A_L_s1s_m[i].invert()*mp_m.get_matrix()*A_L_s0s_m[i];
		a=mismatch_m[0][0]+mismatch_m[1][1];
		b=mismatch_m[0][1]-mismatch_m[1][0];
		a_prime=mismatch_m[0][0]-mismatch_m[1][1];
		b_prime=mismatch_m[0][1]+mismatch_m[1][0];
		coshr=0.5*pow(a*a+b*b,0.5);
		sinhr=0.5*pow(a_prime*a_prime+b_prime*b_prime,0.5);
		p0_sum+=sinhr;
	}
	return p0_sum;
	
}

linalg::vector match_functional_multi::deriv1(interpolator_base* B0)
{
	double infinitesimal=0.00001;
    double denominator=1/(2*infinitesimal);
    double minus, plus;
    interpolator_base* deriv_point=(interpolator_base*)B0->clone();
    
    std::cout<<"\nReceived field is:\n";
    B0->print_y();
    
    std::cout<<"\nCloned deriv_point is:\n";
    deriv_point->print_y();
    interpolator_base* perturbed_point_plus=(interpolator_base*)B0->clone();
    interpolator_base* perturbed_point_minus=(interpolator_base*)B0->clone();

    int num_derivatives=B0->get_num_points_y();  
    vector partials(num_derivatives);    

    for (int j=0; j<num_derivatives; j++)
    {
        perturbed_point_minus->set_basis(j, -1, -infinitesimal, 0);
        
        std::cout<<"\nPerturbed Point Minus is:\n";
        perturbed_point_minus->print_y();
        
        perturbed_point_plus->set_basis(j, -1, infinitesimal, 0);
        
        std::cout<<"\nPerturbed Point Plus is:\n";
        perturbed_point_plus->print_y();

        
        perturbed_point_minus->add(perturbed_point_minus, perturbed_point_minus, deriv_point);
        
        std::cout<<"\nPerturbed Point Minus Add is:\n";
        perturbed_point_minus->print_y();
        
        perturbed_point_plus->add(perturbed_point_plus, perturbed_point_plus, deriv_point);
        
        std::cout<<"\nPerturbed Point Plus Add is:\n";
        perturbed_point_plus->print_y();
        
        minus=this->operator()(*perturbed_point_minus);
        plus=this->operator()(*perturbed_point_plus);
        std::cout<<"\n minus is: "<<minus<<" plus is: "<<plus;
        partials[j]=(plus-minus)/denominator;
    }
    delete deriv_point;
    delete perturbed_point_minus;
    delete perturbed_point_plus;
    return partials;
}

matrix match_functional_multi::deriv2(interpolator_base* B0)
{
    double infinitesimal=0.00001;
    double denominator=1/(infinitesimal*infinitesimal);
    double zero, minus, plus;
    interpolator_base* deriv_point=(interpolator_base*)B0->clone();
    interpolator_base* perturbed_point_plus=(interpolator_base*)B0->clone();
    interpolator_base* perturbed_point_minus=(interpolator_base*)B0->clone();
    int num_derivatives=B0->get_num_points_y();
    matrix hessian(num_derivatives, num_derivatives);     //We do not need derivatives at boundary points
    
    for (int i=0; i<num_derivatives; i++)
    {
        for (int j=0; j<num_derivatives; j++)
        {
            perturbed_point_minus->set_basis(i, j, -infinitesimal, -infinitesimal);
            perturbed_point_plus->set_basis(i, j, infinitesimal, infinitesimal);
            perturbed_point_minus->add(perturbed_point_minus, perturbed_point_minus, deriv_point);
            perturbed_point_plus->add(perturbed_point_plus, perturbed_point_plus, deriv_point);
            minus=this->operator()(*perturbed_point_minus);
            plus=this->operator()(*perturbed_point_plus);
            zero=this->operator()(*deriv_point);
            hessian[i][j]=(plus-2*zero+minus)/denominator;
        }
        
    }
    return hessian;
}

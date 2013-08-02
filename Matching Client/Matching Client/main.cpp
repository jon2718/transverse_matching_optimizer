//
//  main.cpp
//  Matching Client
//
//  Created by Jon Lederman on 11/27/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include <iostream>
#include "amoeba.h"
#include "objectivefunctional.h"
#include "match_functional_multi.h"
#include "symp1d.h"
#include "sep_hamiltonian.h"
#include "bspline_gsl.h"
#include "TApplication.h"
#include "particle_collection.h"
#include "TCanvas.h"
#include "symp1d.h"
#include "eigen.h"


int main(int argc, const char * argv[])
{
	TApplication theApp("App", 0, 0);

	//BSpline
	int k=3;
	int nbreak=10;
	
	//Range
	double s0=0;			//Left
	double s1=10;			//Right
	double B0S0=-10;		//B field left
	double B0S1=10;			//B field right
	
	//Hamiltonians
	double cpe_start=200;
	double qx=0.001;
	double pxp0=0.1;
	int cpe_num=20;			//Number of hamiltonians
	double cpe_step=5;		//Step size of CPE
	double step=0.001;		//Integration step
	
	//Symplectic Integrator
	int dimensions=2;
	int order=4;
	
	//Amoeba
	int iterations=100;
	char amoeba_type[]="right_angled";
	std::string init_simplex(amoeba_type);
	double init_cond[]={-9,1,2,4,-3,3,2,1,9};
	int size=9;
	
	//setup
	symp1d symp_integrator(dimensions, order);
	bspline_gsl *B0_best=new bspline_gsl(k, nbreak, s0, s1);
	bspline_gsl *B0_holder=new bspline_gsl(k, nbreak, s0, s1);
	bspline_gsl *B0_draw=new bspline_gsl(k, nbreak, s0, s1);
	bspline_gsl *B0_const=new bspline_gsl(1, 2, s0, s1);
	
	B0_best->set_bc(-10,10);
	B0_holder->set_bc(-10,10);
	B0_draw->set_bc(-10,10);
	
	//Constant field section
	double coeffs_const[]={10};
	linalg::vector const_vec(coeffs_const, 1);
	B0_const->sety(const_vec, 1);

	match_functional_multi mf(&symp_integrator, s0, s1, step, cpe_start, cpe_num, cpe_step, B0S0, B0S1, qx, pxp0, B0_holder);
	objectivefunctional obj_func(&mf, B0_holder);
	
	//Drawing
	Double_t w = 600;
	Double_t h = 600;
	TCanvas *main_canvas = new TCanvas("Main Canvas", "Multipads", w, h);
	TPad *pad1 = new TPad("p1","titlePad",0,0,1,1);
	pad1->Divide(2,2);
	pad1->SetFillColor(4);
	pad1->Draw();
	
	//Amoeba algorithm
	char name[]="Matching Spline";
	TInterpPlot2D draw_spline(name, B0_draw, 0, 10);
	amoeba amb(iterations, init_simplex, init_cond, size, obj_func,  s0,  s1, &draw_spline, pad1);
	amb.run();
	B0_best->sety(amb.optimal_vector(), size, 1);

    std::cout<<"Optimal vector is: "<<amb.optimal_vector();
    std::cout<<"\nFirst Functional Derivative: \n";
    std::cout<<mf.deriv1(B0_best)<<"\n";
    
    std::cout<<"\nSecond Functional Derivative: \n";
   // std::cout<<mf.deriv2(B0_best)<<"\n";
    
    int num_basis_fns=B0_best->get_num_points_y();
    matrix hessian(num_basis_fns, num_basis_fns);
    hessian=mf.deriv2(B0_best);
    eigen ab(B0_best->get_num_points_y());
	
	hessian.eig_symm(ab);

	std::cout<<ab.get_eigenvalues();

    
    
    
	
	//Track particles through optimum matching section
	sep_hamiltonian best_ham(cpe_start,  qx,  pxp0, B0_best);
	particle_collection myparticles(1000);
	double vals[]={1, 0, 0, 1};
	matrix map(2,2, vals);
	myparticles.gen_ellipse(1, map);
	symp1d tracking_integrator(2,4);
	particle_collection tracked_particles(1000);
	tracked_particles=myparticles.track(s0, s1, best_ham, tracking_integrator, step);
	tracked_particles.draw(pad1, 3);
	
	//Track particles through constant field section
	sep_hamiltonian const_ham(cpe_start,  qx,  pxp0, B0_const);
	particle_collection tracked_const_particles(1000);
	tracked_const_particles=tracked_particles.track(s0, s1, const_ham, tracking_integrator, step);
	tracked_const_particles.draw(pad1, 4);
	theApp.Run();
	
	return 0;
}


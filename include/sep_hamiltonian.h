//
//  sep_hamiltonian.h
//  Phase Space Propagation
//
//  Created by Jon Lederman on 10/9/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Phase_Space_Propagation__sep_hamiltonian__
#define __Phase_Space_Propagation__sep_hamiltonian__

#include <iostream>
#include "interpolate_base.h"

/**********************************************************************************************************************
When you say you have a momentum of 200*10^6 eV/c, what you mean is that pc is 200*10^6*1.6*10^{-16} J.
Thus, pc/e is numerically 200*10^6, and has units of J/C.

 
The user will specify:
	cpe_m (units are Joule/Coulomb)
 	pxpo_m (ratio of x momentum to total momentum (between 0 and 1)
 	qx (initial x displacement).  
  
From this input, the following are calculated:
 	p0_m (actual momentum in kg*m/sec)
 	px_m kg*m/sec
 
A simple example: Take a 10 T solenoid, 200 MeV/c total momentum for the muons. 
Then eB/(2p) = Bc/[2(pc/e)] = 10*3*10^8/(2*200*10^6) = 7.5 m^{-1}. 
Thus, d/ds (px/p) = [eB/(2p)]^2 x = (56.25 m^{-2})
which in general gives you a relatively large change in transverse angle per unit length. 
When you are in a constant solenoid field, you will undergo oscillations such that the two terms 
in the Hamiltonian oscillate such that their maximum values are the same but occur out of sync 
with each other (such that the Hamiltonian remains constant). When the field changes it is more
complicated, which is the point of this exercise in the first place.
 
When you have a momentum of 200*10^6 eV/c, what you mean is that pc is 200*10^6*1.6*10^{-16} J.
Thus, pc/e is numerically 200*10^6, and has units of J/C.
 


A separable hamiltonian.  Assume form is H=V(q)+T(p)
H=alpha*q^2+beta*p^2

Divide Hamiltonian by p0.  Work with phase space variables:  qx and px/p0.
 
H=-p0+(px^2+pz^2)/2p0 + (e^2*b0^2/8p0)*(x^2+z^2)
 
T(p)=(px^2+pz^2)/2p0
V(q)=(e^2*b0^2/8p0)*(x^2+z^2)

Working in 1 variable:

dTdp=px/p0
dVdq=(e^2*b0^2/4p0)*x
 
q_i=q_(i-1) + (tau*c_i)*(p_(i-1)/p0)
p_i=p_(i-1) - (tau*d_i)*(e^2*b0^2/4p0)*(x_i)^2
 
 
Choose more convenient coordinates:
 
 H'=H/p0=-1+(px^2)/(2*p0^2) + (e^2*b0^2/8*p0^2)*(x^2)
 p'=px/p0 and e/p0=r and q=x
 
 H'=-1 + (p'^2)/2  + (b0^2*r^2/8)*q^2
 
 T(p')=(p'^2)/2
 V(q)=(b0^2*r^2/8)*q^2
 
 dTdp=p'
 dVdq=(b0^2*r^2/4)*q
 
 q_i=q_(i-1) + (tau*c_i)*p'
 p'_i=p'_(i-1) - (tau*b0^2*r^2/4)*q_(i-1)
 

***********************************************************************************************************************/

class sep_hamiltonian
{
	
private:
	double p0_m;				//Total momentum
	double cpe_m; 				//cp0/e - units are Joule/Coulomb
	double qx_init_m; 			//Initial condition x
	double pxp0_init_m; 		//Initial condition px/p0.  This is the momentum phase space variable.  
	interpolator_base* B0_m;	//B field
	
	
public:
	sep_hamiltonian(double cpe, double qx, double pxp0, interpolator_base* B0);
	sep_hamiltonian(const sep_hamiltonian& ham);
	~sep_hamiltonian();
	sep_hamiltonian& operator=(const sep_hamiltonian& ham);
	double dVdq(double q, double s) const;
	double dTdp(double p, double s) const;
	void set_B0(interpolator_base *B0);
	interpolator_base* get_B0() const;
	interpolator_base* clone_B0();
	double get_p0() const;
	double get_cpe();
	double get_qx_init();
	double get_pxp0_init();
	
};
#endif /* defined(__Phase_Space_Propagation__sep_hamiltonian__) */

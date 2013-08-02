//
//  integrator.h
//  Phase Space Propagation
//
//  Created by Jon Lederman on 10/19/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Phase_Space_Propagation__integrator__
#define __Phase_Space_Propagation__integrator__

#include <iostream>
#include "sep_hamiltonian.h"
#include "interpolate_base.h"
#include "vector.h"

class integrator  
{
    
protected:
    
    int dimensions_m; 
    
public:
    integrator(const int dimensions);
	virtual ~integrator();
	virtual void integrate(linalg::vector &state, const sep_hamiltonian& h, const double s0, const double s1, const double step)=0;
	
	
};



#endif /* defined(__Phase_Space_Propagation__integrator__) */

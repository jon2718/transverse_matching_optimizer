//
//  objectivefunctional.cpp
//  Metrics
//
//  Created by Jon Lederman on 11/11/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include "objectivefunctional.h"


objectivefunctional::objectivefunctional(functional *fnl, interpolator_base* fn) : fnl_m(fnl), fn_m(fn)
{
	
	
}

/*Sets vector, sets boundary conditions and then returns functional*/

double objectivefunctional::operator()(vector &v)
{
	fn_m->sety(v, fn_m->get_num_points_y()-2, 1);
	return (*fnl_m)(*fn_m);
	
}



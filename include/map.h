//
//  map.h
//  Phase Space Propagation
//
//  Created by Jon Lederman on 11/7/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Phase_Space_Propagation__map__
#define __Phase_Space_Propagation__map__

#include <iostream>
#include "matrix.h"
#include "vector.h"
#include "integrator.h"
#include "interpolate_base.h"

class map
{
private:
	matrix m_m;
    int dimension_m;
	integrator* integrator_m;
    vector sv_m;
	double s0_m, s1_m;				//start and end points of map
	double step_m;
	
public:
    map(int dim, integrator* integrator, double s0, double s1, double step);
    map(const map& mp);                   		//Copy constructor
    map& operator=(const map& map);         	//Assignment operator
    ~map();
    void gen(const sep_hamiltonian& h);
	void gen(const vector& init, const sep_hamiltonian& h);
	void gen(const sep_hamiltonian& h, const double s0, const double s1);
	void gen(const sep_hamiltonian& h, const double s0, const double s1, const vector& init1, const vector& init2);
	matrix calc_AL(const sep_hamiltonian& h, const double s) const;
	integrator* get_intg();
	void set_mapped_pts(const double s0, const double s1);
	void set_step(const double step);
	matrix& get_matrix();
	double get_s0();
	double get_s1();
	void identity();
	friend std::ostream & operator<<(std::ostream & os, const map & m);
};



#endif /* defined(__Phase_Space_Propagation__map__) */

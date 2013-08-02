//
//  particle_collection.h
//  Tracking
//
//  Created by Jon Lederman on 12/15/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Tracking__particle_collection__
#define __Tracking__particle_collection__

#include <iostream>
#include "matrix.h"
#include "sep_hamiltonian.h"
#include "integrator.h"
#include "TGraph.h"
#include "TCanvas.h"

class particle_collection
{
private:
	std::vector<vector> particles_m;
	int num_particles_m;
	TGraph* graph_m;
	
public:
	particle_collection(int num_particles);
	particle_collection(particle_collection& pc);
	~particle_collection();
	void gen_ellipse(double radius, const matrix& map);
	particle_collection track(double s0, double s1, const sep_hamiltonian& h,  integrator &intg, double step);
	particle_collection& operator=(const particle_collection& pc);
	void draw(TPad* pad, int pad_num);
	
};

#endif /* defined(__Tracking__particle_collection__) */

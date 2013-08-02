//
//  particle_collection.cpp
//  Tracking
//
//  Created by Jon Lederman on 12/15/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include "particle_collection.h"
#include "physical_constants.h"
#include "math.h"
#include "map.h"

using namespace physical_constants;

particle_collection::particle_collection(int num_particles) : num_particles_m(num_particles), graph_m(new TGraph(num_particles_m))
{
	vector vec(2);
	vec[0]=0;
	vec[1]=0;
	for (int i=0; i<num_particles_m; i++)
		particles_m.push_back(vec);
	
}

particle_collection::particle_collection(particle_collection& pc) : particles_m(pc.particles_m), num_particles_m(pc.num_particles_m), graph_m(pc.graph_m)
{
	
}


particle_collection::~particle_collection()
{
	delete graph_m;
}

particle_collection& particle_collection::operator=(const particle_collection& pc)
{
	if (this==&pc)
		return (*this);
	else
	
	{
		particles_m=pc.particles_m;
		num_particles_m=pc.num_particles_m;
		*graph_m=*pc.graph_m;
		
	}
	return (*this);
}
void particle_collection::gen_ellipse(double radius, const matrix& map)
{
	//generate unit circle
	vector init(2);
	vector final(2);
	double increment= 2*pi_k/num_particles_m;
	for (int i=0; i<num_particles_m; i++)
	{
		init[0]=cos(i*increment);
		init[1]=sin(i*increment);
		final=map*init;
		graph_m->SetPoint(i, final[0], final[1]);
		particles_m[i]=final;
	}
	
}
particle_collection particle_collection::track(double s0, double s1, const sep_hamiltonian& h,  integrator& intg,  double step)
{
	particle_collection mapped_particles(num_particles_m);
	vector final(2);
	map particle_map(2, &intg, s0, s1, step);
	particle_map.gen(h);
	matrix matrix_map(2,2);
	matrix_map=particle_map.get_matrix();
	std::cout<<matrix_map;
	matrix initial_match(2,2);
	initial_match=particle_map.calc_AL(h, s0);
	std::cout<<initial_match;
	for (int i=0; i<num_particles_m; i++)
	{
		final=matrix_map*initial_match*particles_m[i];
		mapped_particles.graph_m->SetPoint(i, final[0], final[1]);
		mapped_particles.particles_m[i]=final;
	//	std::cout<<mapped_particles.particles_m[i];
		
	}
	return mapped_particles;
}

void particle_collection::draw(TPad* pad, int pad_num)
{
	pad->cd(pad_num);
	graph_m->Draw("AP");
	pad->Update();
	
}

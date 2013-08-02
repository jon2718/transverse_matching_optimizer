//
//  amoeba.h
//  optimizers
//
//  Created by Jon Lederman on 11/13/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __optimizers__amoeba__
#define __optimizers__amoeba__

#include <iostream>
#include "vector.h"
#include <vector>
#include <string>
#include <map>
#include "functional.h"
#include "optimizer.h"
#include "interpolate_base.h"
#include "objective.h"
#include "TCanvas.h"
#include "TInterpPlot2D.h"


class amoeba : public optimizer
{
private:
    int dimension_m;
	std::vector<std::pair<linalg::vector, double> > vectors_m;
	linalg::vector cent, ex, rf, ct;
    int best, second_best, worst, second_worst;
    double alpha, beta, gamma, delta, f_rf, f_xe, f_ct;
    enum termination {term_x, term_f, fail};
 	void rank();
	void init_simplex_right_angled(linalg::vector& x0);
	void init_simplex_regular(linalg::vector& x0);
 	virtual void alg();
	bool reflect();
    bool expand();
 	bool contract();
 	void shrink();
	linalg::vector centroid();
	enum initial{right_angled, regular};
	static std::map<std::string,initial> initial_enum;
	double step_m;
	double best_value_m;
	
public:
	amoeba(int iterations, std::string& init_simplex, double init_cond[], int size, objective& objfn, double s0, double s1, TInterpPlot2D* df, TPad *pad, double a=1.0, double b=0.5, double c=2.0, double d=0.5, double step=1);
	virtual ~amoeba();
	void run();
	void reboot(const linalg::vector& x0);
	friend std::ostream & operator<<(std::ostream & os, const amoeba &amb);
	linalg::vector optimal_vector();
	void Draw(Option_t* option = "");
	void draw();
};

#endif /* defined(__optimizers__amoeba__) */

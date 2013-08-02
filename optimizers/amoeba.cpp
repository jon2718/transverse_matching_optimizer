//
//  amoeba.cpp
//  optimizers
//
//  Created by Jon Lederman on 11/13/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include "amoeba.h"
#include "assign/list_of.hpp"  //Note:: Had to modify boost::assignment_exception.hpp (throw() specifier)
#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"

using namespace boost::assign;
using namespace std;

std::map<std::string,amoeba::initial> amoeba::initial_enum=map_list_of("regular", regular)("right_angled",right_angled);


struct sort_pred {
    bool operator()(const std::pair<linalg::vector,double> &left, const std::pair<linalg::vector,double> &right) {
        return left.second < right.second;
    }
};


amoeba::amoeba(int iterations, std::string& init_simplex, double init_cond[], int size, objective& objfn, double s0, double s1,TInterpPlot2D* df, TPad *pad, double a, double b, double c, double d, double step): optimizer(objfn, init_cond, size, iterations, df, pad), alpha(a), beta(b), gamma(c), delta(d), dimension_m(size+1), cent(size), ex(size), rf(size), ct(size), step_m(step)
{
	best=0;
	second_best=1;
	worst=dimension_m-1;
	second_worst=dimension_m-2;
	amoeba::initial i=(*(amoeba::initial_enum.find(init_simplex))).second;
	
	switch (i)
	{
		case right_angled 		: 	init_simplex_right_angled(init_cond_m);
			break;
		case regular 			: 	init_simplex_regular(init_cond_m);
			break;
	}
}

amoeba::~amoeba()
{}

/*Initial simplex
 The initial simplex S is usually constructed by generating n+1 vertices x0,…,xn around a given input point xin∈Rn . In practice, the most frequent choice is x0=xin to allow proper restarts of the algorithm. The remaining n vertices are then generated to obtain one of two standard shapes of S :
 S is right-angled at x0 , based on coordinate axes, or
 xj:=x0+hjej,j=1,…,n,
 where hj is a stepsize in the direction of unit vector ej in Rn .
 S is a regular simplex, where all edges have the same specified length.
 */

void amoeba::init_simplex_right_angled(linalg::vector& x0)
{
	for (int i=0; i<dimension_m; i++)
		vectors_m.push_back(std::pair<linalg::vector, double>(*(new linalg::vector(dimension_m-1)),0));
	vectors_m[0].first=x0;
	vectors_m[0].second=objfn_m(x0);
    for (int i=1; i<dimension_m; i++)
	{
		vectors_m[i].first=x0+vectors_m[i].first.e(i-1)*step_m;
		vectors_m[i].second=objfn_m(vectors_m[i].first);
	}
	std::sort(vectors_m.begin(), vectors_m.end(), sort_pred());
	best_value_m=vectors_m[0].second;
	pad_m->cd(1);
	this->draw();
	pad_m->cd(2);
}

void amoeba::init_simplex_regular(linalg::vector& x0)
{}


/*Ordering: Determine the indices h,s,l of the worst, second worst and the best vertex, respectively, in the current working simplex S
 fh=maxjfj,fs=maxj≠hfj,fl=minj≠hfj.
 In some implementations, the vertices of S are ordered with respect to the function values, to satisfyf0≤f1≤⋯≤fn−1≤fn . Then l=0 , s=n−1 , and h=n . Consistent tie-breaking rules for this ordering were given by Lagarias et al. (1998).
 */
void amoeba::rank()
{
	
	std::sort(vectors_m.begin(), vectors_m.end(), sort_pred());
	
}

/*Centroid: Calculate the centroid c of the best side—this is the one opposite the worst vertex xh
 c:=1n∑j≠hxj.
 */
linalg::vector amoeba::centroid()
{
//	std::cout<<"\nIn Centroid\n";
	cent.clear();
    for (int i=0; i<dimension_m; i++)
        cent=cent+vectors_m[i].first;
    cent=(1/dimension_m)*cent;
	return cent;
}

void amoeba::run()
{
	alg();
	
}

void amoeba::reboot(const linalg::vector& x0)
{
	for (int i=1; i<dimension_m; i++)
	{
		vectors_m[i].first=x0+vectors_m[i].first.e(i-1)*step_m;
		vectors_m[i].second=objfn_m(vectors_m[i].first);
	}
	std::sort(vectors_m.begin(), vectors_m.end(), sort_pred());
	best_value_m=vectors_m[0].second;
	pad_m->cd(1);
	this->draw();
	pad_m->cd(2);
	std::cout<<"Reboot!";
	
}

void amoeba::alg()
{
	int iterations_since_best_improved=0;
   	for(int i=0; i<iterations_m; i++)
    {
        rank();
		if (iterations_since_best_improved>10)
			{
				reboot(vectors_m[0].first);
				iterations_since_best_improved=0;
			}
		if (vectors_m[0].second<best_value_m)
			{
				std::cout<<"\n"<<vectors_m[0].first<<vectors_m[0].second;
				best_value_m=vectors_m[0].second;
				this->draw();
				iterations_since_best_improved=0;
			}
		iterations_since_best_improved++;
		centroid();
        if (reflect())
            continue;
        else if (expand())
            continue;
        else if (contract())
            continue;
        else
            shrink();
    }
}


//Reflect: Compute the reflection point xr:=c+α(c−xh) and fr:=f(xr) . If fl≤fr<fs , accept xr and terminate the iteration.
bool amoeba::reflect()
{
//	std::cout<<"\nIn Reflect\n";
    rf=cent+alpha*(cent-vectors_m[worst].first);
	f_rf=objfn_m(rf);
    if (vectors_m[best].second<=f_rf && f_rf<vectors_m[second_worst].second)
    {
        vectors_m[worst].first=rf;
		vectors_m[worst].second=f_rf;
        return true;
    }
    else
        return false;
}

/*Expand: If fr<fl , compute the expansion point xe:=c+γ(xr−c) and fe:=f(xe) . If fe<fr , accept xe and terminate the iteration. Otherwise (if fe≥fr), accept xr and terminate the iteration.
 
 This “greedy minimization” approach includes the better of the two points xr , xe in the new simplex, and the simplex is expanded only if fe<fr<fl . It is used in most implementations, and in theory (Lagarias et al., 1998).
 
 The original Nelder-Mead paper uses “greedy expansion”, where xe is accepted if fe<fl and fr<fl , regardless of the relationship between fr and fe . It may happen that fr<fe , so xr would be a better new point than xe , and xe is still accepted for the new simplex. The working simplex is kept as large as possible, to avoid premature termination of iterations, which is sometimes useful for non-smooth functions (see, for example, Rowan, 1990).
 */

bool amoeba::expand()
{
//	std::cout<<"\nIn Expand\n";
    if (f_rf<vectors_m[best].second)
    {
        ex=cent+gamma*(rf-cent);
		f_xe=objfn_m(ex);
		if (f_xe<f_rf)
		{
            vectors_m[worst].first=ex;
			vectors_m[worst].second=f_xe;
        }
		else
		{
            vectors_m[worst].first=rf;
			vectors_m[worst].second=f_rf;
        }
		return true;
    }
    else
        return false;
    
}

/*Contract: If fr≥fs , compute the contraction point xc by using the better of the two points xh and xr .
 Outside: If fs≤fr<fh , compute xc:=c+β(xr−c) and fc:=f(xc) . If fc≤fr , accept xc and terminate the iteration.
 Otherwise, perform a shrink transformation.
 
 Inside: If fr≥fh , compute xc:=c+β(xh−c) and fc:=f(xc) . If fc<fh , accept xc and terminate the iteration.
 Othwise, perform a shrink transformation.
 */

bool amoeba::contract()
{
//	std::cout<<"\nIn Contract\n";
	if (f_rf>=vectors_m[second_worst].second)
	{
		if (vectors_m[second_worst].second<=f_rf && f_rf<vectors_m[worst].second)
		{
			ct=cent+beta*(rf-cent);
			f_ct=objfn_m(ct);
			if (f_ct<=f_rf)
			{
				vectors_m[worst].first=ct;
				vectors_m[worst].second=f_ct;
				return true;
			}
			else
				return false;
		}
		else if (f_rf>=vectors_m[worst].second)
		{
			ct=cent+beta*(vectors_m[worst].first-cent);
			f_ct=objfn_m(ct);
			if (f_ct<vectors_m[worst].second)
			{
				vectors_m[worst].first=ct;
				vectors_m[worst].second=f_ct;
				return true;
			}
			else
				return false;
		}
	}
	else
		return false;
}


/*Shrink: Compute n new vertices xj:=xl+δ(xj−xl) and fj:=f(xj) , for j=0,…,n , with j≠l .
 The shrink transformation was introduced to prevent the algorithm from failing in the following case, described by the quote from the original paper:
 A failed contraction is much rarer, but can occur when a valley is curved and one point of the simplex is much farther from the valley bottom than the others; contraction may then cause the reflected point to move away from the valley bottom instead of towards it. Further contractions are then useless. The action proposed contracts the simplex towards the lowest point, and will eventually bring all points into the valley.
 */

void amoeba::shrink()
{
//	std::cout<<"\nIn Shrink\n";
	for (int i=0; i<dimension_m; i++)
	{
		if (i != best)
		{
			vectors_m[i].first=vectors_m[best].first+delta*(vectors_m[i].first-vectors_m[best].first);
			vectors_m[i].second=objfn_m(vectors_m[i].first);
		}
	}
}

std::ostream & operator<<(std::ostream & os, const amoeba & amb)
{
	std::cout<<"\n";
	for (int i=0;i<amb.dimension_m; i++)
	{
		os<<amb.vectors_m[i].first<<" "<<amb.vectors_m[i].second<<"\n";
	}
	return os;
	
}

linalg::vector amoeba::optimal_vector()
{
	return vectors_m[best].first;
	
}

void amoeba::draw()
{
	draw_function->sety(vectors_m[0].first, dimension_m-1, 1);
	draw_function->Draw();
	pad_m->Update();
	
}

void amoeba::Draw(Option_t* option)
{
	draw_function->sety(vectors_m[0].first, dimension_m-1, 1);
	draw_function->Draw();
//	canvas_m->Update();
}


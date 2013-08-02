//
//  optimizer.h
//  optimizers
//
//  Created by Jon Lederman on 11/13/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __optimizers__optimizer__
#define __optimizers__optimizer__

#include <iostream>
#include "objective.h"
#include "TObject.h"
#include "TCanvas.h"
#include "interpolate_base.h"
#include "TInterpPlot2d.h"
#include "vector.h"


class optimizer 
{
    
protected:
    int iterations_m;
	objective& objfn_m;								//objective function
	linalg::vector init_cond_m;
	TInterpPlot2D* draw_function;
	TPad *pad_m;
	
	
public:
	optimizer(objective& objfn, const double init[], const int size, const int iterations, TInterpPlot2D* df, TPad* pad);
    virtual ~optimizer() {};
    optimizer & operator=(const optimizer &opt);    //assignment operator
    virtual void alg()=0;
    
};



#endif /* defined(__optimizers__optimizer__) */

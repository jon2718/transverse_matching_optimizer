//
//  optimizer.cpp
//  optimizers
//
//  Created by Jon Lederman on 11/13/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include "optimizer.h"




optimizer::optimizer(objective& objfn, const double init[], const int size, const int iterations, TInterpPlot2D* df, TPad* pad) : iterations_m(iterations), init_cond_m(init, size), objfn_m(objfn), draw_function(df), pad_m(pad)
{	
	
}
//
//  objective.h
//  Metrics
//
//  Created by Jon Lederman on 11/11/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Metrics__objective__
#define __Metrics__objective__

#include <iostream>
#include "vector.h"
using namespace linalg;

class objective
{
    
protected:
    
public:
    objective () {};
	//objectivef (const objective& obf);
	
    virtual ~objective() {};
 	//objectivef & operator=(const objective &opt);
    virtual double operator()(vector &v)=0;
    
};


#endif /* defined(__Metrics__objective__) */

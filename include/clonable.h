//
//  clonable.h
//  Design Patterns
//
//  Created by Jon Lederman on 12/29/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#ifndef __Design_Patterns__clonable__
#define __Design_Patterns__clonable__

#include <iostream>
struct clonable
{
    virtual ~clonable() {}
    virtual clonable* clone() const = 0;
};
#endif /* defined(__Design_Patterns__clonable__) */

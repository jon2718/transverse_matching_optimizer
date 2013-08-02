//
//  tracking-test.cpp
//  Tracking
//
//  Created by Jon Lederman on 12/18/12.
//  Copyright (c) 2012 Jon Lederman. All rights reserved.
//

#include "gtest/gtest.h"
#include "particle_collection.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "matrix.h"

TEST(tracking, evaluate)
{
		TApplication theApp("Tracking Test App", 0, 0);
	
		particle_collection myparticles(1000);
		double vals[]={1, 2, 3, 4};
		matrix map(2,2, vals);
		myparticles.gen_ellipse(1, map);
	
		//TCanvas canvas1("Canvas 1", "Graph Draw Options", 200, 10, 600, 400);
		
		Double_t w = 600;
		Double_t h = 600;
		TCanvas *canvas1 = new TCanvas("c", "c", w, h);
		myparticles.draw(canvas1, 1);
		theApp.Run();
	
	
}

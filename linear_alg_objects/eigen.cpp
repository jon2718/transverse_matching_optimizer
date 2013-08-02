//
//  eigen.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "eigen.h"


eigen::eigen(const matrix& mat) : eigenvalues(mat.get_rows()), eigenvectors(mat.get_rows(), mat.get_columns())
{}

eigen::eigen(int size) : eigenvalues(size), eigenvectors(size, size)
{}

eigen::eigen(const eigen& eig) : eigenvalues(eig.get_eigenvalues()), eigenvectors(eig.get_eigenvectors())
{}

eigen::~eigen()
{}

eigen& eigen::operator=(const eigen& eig)
{
	eigenvectors=eig.eigenvectors;
	eigenvalues=eig.eigenvalues;
	return (*this);
}

const matrix& eigen::get_eigenvectors() const
{
	return eigenvectors;	
}

const vector& eigen::get_eigenvalues() const
{
	return eigenvalues;
}

void eigen::set_eigenvalues(gsl_vector &vec)
{
	for (int i=0; i<eigenvalues.get_size(); i++)
		eigenvalues[i]=gsl_vector_get(&vec, i);
	
}

void eigen::set_eigenvectors(gsl_matrix &mat)
{
	for (int i=0; i<eigenvalues.get_size(); i++)
	{
		for (int j=0; j<eigenvalues.get_size(); j++)
			eigenvectors[i][j]=gsl_matrix_get(&mat, i, j);
	}

}

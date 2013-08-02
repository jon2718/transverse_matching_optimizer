//
//  matrixc.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/14/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "matrixc.h"
//#include "Accelerate/Accelerate.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

matrixc::row_proxy::row_proxy(matrixc& mat, int row) : mat_m(mat), row_m(row)
{}

matrixc::row_proxy::column_proxy matrixc::row_proxy::operator [](int row)
{
	return column_proxy(*this, row);
	
}

matrixc::row_proxy::column_proxy::column_proxy(row_proxy& rp, int column) : rp_m(rp), column_m(column)
{}

matrixc::row_proxy::column_proxy& matrixc::row_proxy::column_proxy::operator=(gsl_complex& val)
{

	gsl_matrix_complex_set(rp_m.mat_m.g_m, rp_m.row_m, this->column_m , val);
	return (*this);
	
}

matrixc::row_proxy::row_proxy(const row_proxy& rhs) : mat_m(rhs.mat_m), row_m(rhs.row_m)
{}

matrixc::row_proxy& matrixc::row_proxy::operator=(const row_proxy& rhs)
{
	if (this==&rhs)
		return *this;
	mat_m=rhs.mat_m;
	row_m=rhs.row_m;
	return (*this);
}

matrixc::row_proxy::column_proxy& matrixc::row_proxy::column_proxy::operator=(const column_proxy& rhs)
{
	column_m=rhs.column_m;
	rp_m=rhs.rp_m;
	return (*this);
}


matrixc::row_proxy::column_proxy::operator gsl_complex() 
{
	return gsl_matrix_complex_get (rp_m.mat_m.g_m, rp_m.row_m, this->column_m);
}


matrixc::matrixc(int rows, int columns, const double values[]) : rows_m(rows), columns_m(columns), g_m(gsl_matrix_complex_alloc(rows, columns))
{
	try {
		if (rows<=0)
			throw("Invalid number of rows");
		if (columns<=0)
			throw("Invalid number of columns");
	} catch (char* str) 
	{
		std::cout<<"Invalid argument to matrix constructor.  "<<str;
	}
	gsl_complex temp;
	for (int i=0; i<rows_m; i++)
	{
		for (int j=0; j<columns_m; j++)
		{
			GSL_SET_COMPLEX(&temp, values[2*i*columns+j], values[2*i*columns+j+1]);
			gsl_matrix_complex_set(g_m, i, j, temp); 
		}
	}
}

matrixc::matrixc(int rows, int columns) : rows_m(rows), columns_m(columns), g_m(gsl_matrix_complex_alloc(rows, columns))
{
	gsl_matrix_complex_set_zero(g_m);
}

matrixc::matrixc(const gsl_matrix_complex &gsl_mat) : g_m(gsl_matrix_complex_alloc(gsl_mat.size1, gsl_mat.size2)), rows_m(gsl_mat.size1), columns_m(gsl_mat.size2)
{
	gsl_matrix_complex_memcpy (g_m, &gsl_mat);	
}

matrixc::~matrixc()
{
	gsl_matrix_complex_free(g_m);
}

matrixc::matrixc(const matrixc& mat) : rows_m(mat.rows_m), columns_m(mat.columns_m), g_m(gsl_matrix_complex_alloc(mat.rows_m, mat.columns_m))
{
	for (int i=0; i<rows_m; i++)
	{
		for (int j=0; j<columns_m; j++)
			gsl_matrix_complex_set(g_m, i, j, gsl_matrix_complex_get(g_m, i, j)); 
	}	

}


matrixc& matrixc::operator=(const matrixc& mat) 
{
	if (this==&mat)
		return *this;
	gsl_matrix_complex_free(g_m);
	rows_m=mat.rows_m;
	columns_m=mat.columns_m;
	
	g_m=gsl_matrix_complex_alloc(rows_m, columns_m);
	
	for (int i=0; i<rows_m; i++)
	{
		for (int j=0; j<columns_m; j++)
			gsl_matrix_complex_set(g_m, i, j, gsl_matrix_complex_get(g_m, i, j)); 
	}	
	return *this;	
	
}

matrixc::row_proxy matrixc::operator[](double row)
{
	return row_proxy(*this, row);
	
}


matrixc matrixc::operator*(const matrixc& mat)
{
	gsl_complex alpha, beta;
	gsl_matrix_complex* C=gsl_matrix_complex_alloc(this->rows_m, mat.columns_m);
	GSL_SET_COMPLEX(&alpha, 1, 0);
	GSL_SET_COMPLEX(&beta, 1, 0);

	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, g_m, mat.g_m, beta, C);
	return matrixc(*C);
}

bool matrixc::operator==(matrixc& mat)
{
	return  gsl_matrix_complex_equal(g_m, mat.g_m);
}

int matrixc::get_rows() const
{
	return rows_m;
}

int matrixc::get_columns() const
{
	return columns_m;
}


void matrixc::eig(eigenc& eig)
{	
	
		
}


std::ostream & operator<<(std::ostream & os,  matrixc & mat)
{
	os<<"\n{";
	for(int i=0; i<mat.rows_m; i++)
	{
		for(int j=0; j<mat.columns_m; j++)
		{
//			os<<mat[i][j];
			if (j!=mat.columns_m-1)
				std::cout<<" ";
		}
		if (i!=mat.rows_m-1)
			os<<"\n";
	}
	os<<"}\n";
	return os;
};



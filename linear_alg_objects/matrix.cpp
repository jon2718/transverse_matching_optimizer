//
//  matrix.cpp
//  linear_alg_objects
//
//  Created by Jon Lederman on 6/29/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "matrix.h"
#include "gsl/gsl_blas.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include "eigen.h"
#include "matrix_exceptions.h"

matrix::row_proxy::row_proxy(matrix& mat, int row) : mat_m(mat), row_m(row)
{}


matrix::row_proxy::column_proxy matrix::row_proxy::operator [](int row)
{
	return column_proxy(*this, row);
		
}

const matrix::row_proxy::column_proxy matrix::row_proxy::operator[](int row) const
{
	return column_proxy(const_cast<row_proxy&>(*this), row);
	
}

matrix::row_proxy::column_proxy::column_proxy(row_proxy& rp, int col) : rp_m(rp), column_m(col)
{}


matrix::row_proxy::column_proxy& matrix::row_proxy::column_proxy::operator=(double val)
{
	gsl_matrix_set (rp_m.mat_m.g_m, rp_m.row_m, this->column_m , val);
	return (*this);
}

matrix::row_proxy::row_proxy(const row_proxy& rhs) : mat_m(rhs.mat_m), row_m(rhs.row_m)
{
}

matrix::row_proxy& matrix::row_proxy::operator=(const row_proxy& rhs)
{
	if (this==&rhs)
		return *this;
	mat_m=rhs.mat_m;
	row_m=rhs.row_m;
	return (*this);
}

matrix::row_proxy::column_proxy& matrix::row_proxy::column_proxy::operator=(const column_proxy& rhs)
{
	column_m=rhs.column_m;
	rp_m=rhs.rp_m;
	return (*this);
}


matrix::row_proxy::column_proxy::operator double() const
{
	return gsl_matrix_get(rp_m.mat_m.g_m, rp_m.row_m, this->column_m);
}


matrix::matrix(int rows, int columns, const double values[]) : rows_m(rows), columns_m(columns), g_m(gsl_matrix_alloc(rows, columns))	
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
	
	for (int i=0; i<rows_m; i++)
	{
		for (int j=0; j<columns_m; j++)
			{
				gsl_matrix_set(g_m, i, j, values[i*columns_m+j]);
			}
		
	}
}

matrix::matrix(int rows, int columns) : rows_m(rows), columns_m(columns), g_m(gsl_matrix_alloc(rows, columns))
{
	gsl_matrix_set_zero(g_m);
}

matrix::matrix(gsl_matrix& mat, int rows, int columns) : rows_m(rows), columns_m(columns), g_m(gsl_matrix_alloc(rows, columns))
{
		gsl_matrix_memcpy(g_m, &mat);
}

matrix::~matrix()
{
	gsl_matrix_free(g_m);
}

matrix::matrix(const matrix& mat) : rows_m(mat.rows_m), columns_m(mat.columns_m), g_m(gsl_matrix_alloc(mat.rows_m, mat.columns_m))
{
	gsl_matrix_memcpy(g_m, mat.g_m);
}
		

matrix& matrix::operator=(const matrix& mat) 
{
	if (this==&mat)
		return *this;
	gsl_matrix_free(g_m);
	rows_m=mat.rows_m;
	columns_m=mat.columns_m;
	g_m=gsl_matrix_alloc(rows_m, columns_m);
	gsl_matrix_memcpy(g_m, mat.g_m);
	return *this;	
}

matrix::row_proxy matrix::operator[](int row)
{
	return row_proxy(*this, row);		
}

const matrix::row_proxy matrix::operator[](int row) const
{
	return row_proxy(const_cast<matrix&>(*this), row);
	
}

matrix matrix::operator*(matrix &mat)
{
	
	gsl_matrix *C=gsl_matrix_alloc(rows_m, mat.columns_m);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, g_m, mat.g_m, 0, C);
	
	matrix output(*C,rows_m, mat.columns_m);
	return output;
}

vector matrix::operator*(const vector& vec) const
{
	vector result(vec.get_size());
	gsl_vector* y=gsl_vector_calloc(vec.get_size());
	gsl_blas_dgemv(CblasNoTrans, 1, g_m, vec.get_gsl_vector(), 1, y);
	result.set_values(y);
	return result;
}
	
	

bool matrix::operator==(matrix &mat)
{
 if (gsl_matrix_equal(g_m, mat.g_m))
	return true;
	else 
		return false;
}

void matrix::identity()
{
	gsl_matrix_set_identity(g_m);
}

/*Rewrite inversion as per GSL note:
These functions compute the inverse of a matrix A from its LU decomposition (LU,p), storing the result in the matrix inverse. The inverse is computed by solving the system A x = b for each column of the identity matrix. It is preferable to avoid direct use of the inverse whenever possible, as the linear solver functions can obtain the same result more efficiently and reliably (consult any introductory textbook on numerical linear algebra for details).
 */
 
matrix matrix::invert()
{
	matrix LU(*this);
	matrix inversion=matrix(rows_m, columns_m);
	gsl_permutation * p = gsl_permutation_alloc(rows_m);
	int s;
	gsl_linalg_LU_decomp(LU.g_m, p, &s);
	gsl_linalg_LU_invert (LU.g_m, p, inversion.g_m);
	return inversion;
}

void matrix::set_column(int index, const vector& col)
{
	for (int i=0; i<rows_m; i++)
		gsl_matrix_set(g_m, i, index, col[i]);
}

void matrix::set_row(int index, const vector &row)
{
	for (int i=0; i<columns_m; i++)
		gsl_matrix_set(g_m, index, i, row[i]);
}

int matrix::get_rows() const
{
	return rows_m;
}

int matrix::get_columns() const
{
	return columns_m;
}


void matrix::eig_symm(eigen& eig)
{
	if (!is_symmetric())
		throw NotSymmetricMatrixException();
	gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc(rows_m);
	gsl_matrix *a=gsl_matrix_alloc(rows_m, columns_m);
	gsl_matrix *evec=gsl_matrix_alloc(rows_m, rows_m);
	gsl_vector *eval=gsl_vector_alloc(rows_m);
	gsl_matrix_memcpy(a, g_m);
	gsl_eigen_symmv(a, eval, evec, w);
	eig.set_eigenvalues(*eval);
	eig.set_eigenvectors(*evec);
	
	
	
}
void matrix::eig_non_symm(eigenc& eig)
{
	gsl_eigen_nonsymm_workspace *workspace= gsl_eigen_nonsymm_alloc (rows_m);
	gsl_vector_complex* eigenvalues=gsl_vector_complex_alloc(rows_m);
	gsl_matrix_complex* eigenvectors=gsl_matrix_complex_alloc(rows_m, columns_m);
		
}


bool matrix::is_symmetric()
{
	gsl_matrix* trans=gsl_matrix_alloc(rows_m, rows_m);
	gsl_matrix_transpose_memcpy(trans,g_m);
	if (gsl_matrix_equal(g_m, trans))
		return true;
	else
		return false;
	
}

std::ostream & operator<<(std::ostream & os, const matrix & mat)
{
	os<<"\n{";
	for(int i=0; i<mat.rows_m; i++)
	{
		for(int j=0; j<mat.columns_m; j++)
		{
			os<<mat[i][j];
			if (j!=mat.columns_m-1)
				std::cout<<" ";
		}
		if (i!=mat.rows_m-1)
		os<<"\n";
	}
	os<<"}\n";
	return os;
}



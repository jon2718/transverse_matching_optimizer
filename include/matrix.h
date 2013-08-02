//
//  matrix.h
//  linear_alg_objects
//
//  Created by Jon Lederman on 6/29/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef linear_alg_objects_matrix_h
#define linear_alg_objects_matrix_h
#include "vector.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
using namespace linalg;

class eigen;
class eigenc;

class matrix
{
private:
	int rows_m, columns_m;
	gsl_matrix* g_m;
	
	
public:
	class row_proxy
	{
		private: 
			matrix& mat_m;
			int row_m;
		public:
			class column_proxy
				{
					private:
						row_proxy& rp_m;
						int column_m;
					public:
					column_proxy(row_proxy& rp, int col);
					operator double() const;									//r-value uses
					column_proxy& operator=(const column_proxy& rhs); 	//lvalue uses
					column_proxy& operator=(double val); 
				};
		
			row_proxy(matrix& mat, int row);
			row_proxy(const row_proxy& rhs);
			row_proxy& operator=(const row_proxy& rhs);
			column_proxy operator[](int row);
			const column_proxy operator[](int row) const;
		friend class column_proxy;
		
	};
	matrix(int rows, int columns);
	matrix(int rows, int columns, const double values[]);
	matrix(gsl_matrix& mat, int rows, int columns);
	matrix(const matrix& mat);
	~matrix();
	matrix& operator=(const matrix& mat);
	row_proxy operator[](int row);
	const row_proxy operator[](int row) const; 	//for const Matrices
	matrix operator*(matrix &mat);
	vector operator*(const vector& vec) const;
	bool operator==(matrix &mat);
	void identity();
	matrix invert();
	void set_column(int index, const vector& col);
	void set_row(int index, const vector &row);
	int get_rows() const;
	int get_columns() const;
	void eig_symm(eigen& eig);
	void eig_non_symm(eigenc& eig);
	bool is_symmetric();
	friend std::ostream & operator<<(std::ostream & os, const matrix & mat);
	friend class row_proxy;
};


#endif

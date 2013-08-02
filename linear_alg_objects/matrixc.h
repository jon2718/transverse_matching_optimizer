//
//  matrixc.h
//  linear_alg_objects
//
//  Created by Jon Lederman on 7/14/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef linear_alg_objects_matrixc_h
#define linear_alg_objects_matrixc_h
#include "vector.h"
#include "gsl/gsl_vector_complex.h"
#include "gsl/gsl_matrix.h"
class eigenc;

class matrixc
{
private:
	int rows_m, columns_m;
	gsl_matrix_complex* g_m;
	
public:
	class row_proxy
	{
	private: 
		matrixc& mat_m;
		int row_m;
	public:
		class column_proxy
		{
		private:
			row_proxy& rp_m;
			int column_m;
		public:
			column_proxy(row_proxy& rp, int col);
			operator gsl_complex();									//r-value uses
			column_proxy& operator=(const column_proxy& rhs); 	//lvalue uses
			column_proxy& operator=(gsl_complex& val); 
		};
		row_proxy(matrixc& mat, int row);
		row_proxy(const row_proxy& rhs);
		row_proxy& operator=(const row_proxy& rhs);
		column_proxy operator[](int row);
	};
	matrixc(int rows, int columns);
	matrixc(int rows, int columns, const double values[]);
	matrixc(const matrixc& mat);
	matrixc(const gsl_matrix_complex &gsl_mat);
	~matrixc();
	matrixc& operator=(const matrixc& mat);
	row_proxy operator[](double row);
	matrixc operator*(const matrixc &mat);
	bool operator==(matrixc &mat);
	int get_rows() const;
	int get_columns() const;
	void eig(eigenc& eig);
	friend std::ostream & operator<<(std::ostream & os,  matrixc& mat);	
};



#endif

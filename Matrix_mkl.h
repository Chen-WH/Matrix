#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <fstream>
#include <map>
#include "mkl.h"

using namespace std;

class Matrix
{
private:
	int row, col;
	double* index;
public:
	Matrix(int x = 1, int y = 1) {
		row = x;
		col = y;
		index = (double *)mkl_malloc( row*col*sizeof( double ), 64 );
	}

	Matrix(const Matrix& obj) {
		row = obj.row;
		col = obj.col;
		index = (double *)mkl_malloc( row*col*sizeof( double ), 64 );
		cblas_dcopy(row*col, obj.index, 1, index, 1);
	}

	Matrix(double **array, int x, int y) {
		row = x;
		col = y;
		index = (double *)mkl_malloc( row*col*sizeof( double ), 64 );
		for (int i = 0; i < row; ++i) {
			cblas_dcopy(row*col, array[i], 1, index + i*col, 1);
		}
	}

	~Matrix() {
		mkl_free(index);
	}

	void E() {
		for (int i = 0; i < (row*col); ++i) {
			index[i] = 0;
		}
		for (int i = 0; i < row; ++i) {
			index[i*col+i] = 1;
		}
	}

	void zero() {
		for (int i = 0; i < (row*col); ++i) {
			index[i] = 0;.
	}

	friend istream& operator>>(istream& is, Matrix& obj) {
		for (int i = 0; i < (obj.row*obj.col); ++i) {
			is >> obj.index[i];
		}
		return is;
	}

	friend ostream& operator<<(ostream& os, const Matrix& obj) {
		if (obj.row == 0 || obj.col == 0) {
			os << "ERROR!";          //不符合矩阵运算时自动报错
		}
		else {
			for (int i = 0; i < obj.row; ++i) {
				for (int j = 0; j < obj.col; ++j) {
					os << obj.index[i*obj.col + j] << ' ';
				}
				os << endl;
			}
		}
		return os;
	}

	int Row() {
		return row;
	}

	int Col() {
		return col;
	}

	double& operator()(int x, int y) {
		return index[x*col + y];
	}

	// *this = *this + obj
	Matrix add(const Matrix obj) {
		if (obj.row != row || obj.col != col) {
			printf("相加矩阵维数不同!\n");
			return *this;
		}
		else {
			cblas_daxpy(row*col, 1, obj.index, 1, index, 1);
			return *this;
		}
	}

	// *this = *this - obj
	Matrix minus(const Matrix& obj) {
		if (obj.row != row || obj.col != col) {
			printf("相加矩阵维数不同!\n");
			return *this;
		}
		else {
			cblas_daxpy(row*col, -1, obj.index, 1, index, 1);
			return *this;
		}
	}

	// *this = *this * obj
	Matrix multp(const Matrix& obj) {
		if (col != obj.row) {
			cout << "维数不符合矩阵相乘原则！" << endl;
		}
		else {
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, obj.col, col, 1, index, col, obj.index, obj.col, 0, index, obj.col);
			return *this;
		}
	}

	Matrix multp(const double x) {
		cblas_dscal(row*col, x, index, 1);
		return *this;
	}

	Matrix div(const double x) {
		cblas_dscal(row*col, 1/x, index, 1);
		return *this;
	}

	Matrix operator+(const Matrix& obj) {
		Matrix tmp(row, col);
		if (obj.row != row || obj.col != col) {
			printf("相加矩阵维数不同!\n");
			return *this;
		}
		else {
			tmp.zero();
			cblas_daxpy(row*col, 1, obj.index, 1, tmp.index, 1);
			cblas_daxpy(row*col, 1, index, 1, tmp.index, 1);
			return tmp;
		}
	}

	// *this = *this - obj
	Matrix operator-(const Matrix& obj) {
		Matrix tmp(row, col);
		if (obj.row != row || obj.col != col) {
			printf("相加矩阵维数不同!\n");
			return *this;
		}
		else {
			tmp.zero();
			cblas_daxpy(row*col, -1, obj.index, 1, tmp.index, 1);
			cblas_daxpy(row*col, 1, index, 1, tmp.index, 1);
			return tmp;
		}
	}

	friend Matrix operator*(const Matrix& obj1, const Matrix& obj2) {
		if (obj1.col != obj2.row) {
			cout << "维数不符合矩阵相乘原则！" << endl;
			Matrix tmp;
			return tmp;
		}
		else {
			Matrix tmp(obj1.row, obj2.col);
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, obj1.row, obj2.col, obj1.col, 1, obj1.index, obj1.col, obj2.index, obj2.col, 0, tmp.index, obj2.col);
			return tmp;
		}
	}

	friend Matrix operator*(const Matrix& obj, const double x) {
		Matrix tmp(obj);
		cblas_dscal(obj.row*obj.col, x, tmp.index, 1);
		return tmp;
	}

	friend Matrix operator/(const Matrix& obj, const double x) {
		Matrix tmp(obj);
		cblas_dscal(obj.row*obj.col, 1/x, tmp.index, 1);
		return tmp;
	}

	Matrix& operator=(const Matrix& obj) {
		mkl_free(index);
		Matrix *this(obj);
		return *this;
	}

	//矩阵转置
	Matrix trans() {
		Matrix tmp(col, row);
		for (int i = 0; i < col; ++i) {
			cblas_dcopy(row, index + i, col, tmp.index + i*row, 1);
		}
		return tmp;
	}

	//分块矩阵
	Matrix block(int x, int y, int new_row, int new_col) {
		Matrix tmp(new_row, new_col);
		for (int i = 0; i < new_row; ++i) {
			for (int j = 0; j < new_col; ++j) {
				tmp.index[i][j] = index[x + i][y + j];
			}
			cblas_dcopy(row, index + i, col, tmp.index + i*row, 1);
		}
		return tmp;
	}

	//合并分块矩阵
	

	//行交换
	void swaprow(const int i, const int j) {
		double* tmp = index[i];
		index[i] = index[j];
		index[j] = tmp;
	}

	//行变换
	void rowmultik(int k, double x) {
		for (int i = 0; i < col; ++i) {
			index[k][i] *= x;
		}
	}

	//求矩阵任意元素的余子式
	double cofactor(const int x, const int y) {
		double tmp;
		Matrix M1(row - 1, row - 1);
		for (int i = 0; i < row; ++i) {
			if (i == x) {
				continue;
			}
			for (int j = 0; j < col; ++j) {
				if (j == y) {
					continue;
				}
				if (i < x) {
					if (j < y) {
						M1.index[i][j] = index[i][j];
					}
					else {
						M1.index[i][j - 1] = index[i][j];
					}
				}
				else {
					if (j < y) {
						M1.index[i - 1][j] = index[i][j];
					}
					else {
						M1.index[i - 1][j - 1] = index[i][j];
					}
				}
			}
		}
		tmp = M1.value();
		return tmp;
	}

	//求行列式的值
	double value() {
		Matrix tmp(*this);
		int count = 0;
		int i, j, k;
		double det = 1, tp;
		for (i = 0; i < row; ++i) {
			if (tmp.index[i][i] == 0) {
				for (j = i + 1; j < col; ++j) {
					if (tmp.index[j][i] != 0) {
						tmp.swaprow(i, j);
						++count;
						break;
					}
				}
				if (j == col) {
					return 0;
				}
			}
			for (j = i + 1; j < row; ++j) {
				if (tmp.index[j][i] != 0) {
					tp = tmp.index[j][i] / tmp.index[i][i];
					for (k = i; k < col; ++k) {
						tmp.index[j][k] -= tp * tmp.index[i][k];
					}
				}
			}
		}
		for (i = 0; i < row; ++i) {
			det *= tmp.index[i][i];
		}
		if (count % 2 == 1)
			det = -det;
		return det;
	}

	//判断是否正定
	bool ifposdef() {
		for (int i = 1; i <= row; ++i) {
			if (block(0, 0, i, i).value() <= 0) {
				return false;
			}
		}
		return true;
	}

	//矩阵修正
	void modify() {
		double b1 = 0, tp;
		for (int i = 0; i < row; ++i) {
			tp = 2 * index[i][i];
			for (int j = 0; j < col; ++j) {
				tp -= index[i][j];
			}
			b1 = (tp < b1) ? tp : b1;
		}
		for (int i = 0; i < row; ++i) {
			index[i][i] -= b1;
		}
	}

	//求伴随矩阵
	Matrix adjoint() {
		Matrix tmp(row, col);
		for (int i = 0; i < tmp.row; ++i) {
			for (int j = 0; j < tmp.col; ++j) {
				if ((i + j) % 2 == 0) {
					tmp.index[i][j] = cofactor(j, i);
				}
				else {
					tmp.index[i][j] = cofactor(j, i) * (-1);
				}
			}
		}
		return tmp;
	}

	//求解线性方程组
	Matrix solve(const Matrix b) {
		Matrix A(*this);
		Matrix b2(b);
		int i, j;
		for (i = 0; i < A.row; ++i) {
			j = i;
			while (j < A.row && A.index[j][i] == 0) {
				++j;
			}
			if (j == A.row) {
				cout << "wrong" << endl;
				A.modify();
				return A.solve(b);
			}
			else {
				if (i != j) {
					A.swaprow(i, j);
					b2.swaprow(i, j);
				}
			}
			//将对角位置转换为1
			double tmp = 1.0 / A.index[i][i];
			b2.rowmultik(i, tmp);
			A.rowmultik(i, tmp);
			//将该列非对角位置转换为0
			for (int k = 0; k < A.row; k++) {
				if (k == i || A.index[k][i] == 0) {
					continue;
				}
				for (j = i + 1; j < A.col; ++j) {
					A(k, j) -= A(i, j) * A(k, i);
				}
				b2(k, 0) -= b2(i, 0) * A(k, i);
				A(k, i) = 0;
			}
		}
		return b2;
	}

	//求逆矩阵
	Matrix inv() {
		Matrix tmp(*this);
		Matrix eye(row, col);
		eye.E();
		int i, j, k;
		for (i = 0; i < row; ++i) {
			//判断对角方向的元素是否为0
			j = i;
			while (j < row && tmp.index[j][i] == 0) {
				++j;
			}
			if (j == row) {
				cout << "该矩阵不可逆！" << endl;
				return eye;
			}
			else {
				if (i != j) {
					tmp.swaprow(i, j);
					eye.swaprow(i, j);
				}
			}
			//将对角位置转换为1
			double tp = 1.0 / tmp(i, i);
			eye.rowmultik(i, tp);
			tmp.rowmultik(i, tp);
			//将该列非对角位置转换为0
			for (k = 0; k < row; k++) {
				if (k == i || tmp.index[k][i] == 0) {
					continue;
				}
				for (j = i + 1; j < col; j++) {
					tmp.index[k][j] -= tmp(i, j) * tmp(k, i);
				}
				for (j = 0; j < col; j++) {
					eye.index[k][j] -= eye(i, j) * tmp(k, i);
				}
				tmp.index[k][i] = 0;
			}
		}
		return eye;
	}

	//一维范数
	double norm1() {
		double tmp = 0;
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				tmp += abs(index[i][j]);
			}
		}
		return tmp;
	}

	//二维范数
	double norm2() {
		double tmp = 0;
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				tmp += index[i][j] * index[i][j];
			}
		}
		return sqrt(tmp);
	}
	
	//求取两向量夹角
	friend double theta(Matrix obj1, Matrix obj2)
	{
		if (obj1.row != obj2.row || obj1.col != 1 || obj2.col != 1) {
			cout << "向量维数不相等！" << endl;
			return 0;
		}
		double tmp = 0;
		for (int i = 0; i < obj1.row; ++i) {
			tmp += obj1.index[i][0] * obj2.index[i][0];
		}
		return tmp / obj1.norm2() / obj2.norm2();
	}
};
#ifndef CHOLESKY_LDLT
#define CHOLESKY_LDLT
#define Mat Matrix<_T, Dynamic, Dynamic>

#include <vector>//import std::vector
#include <Eigen/Core>

using std::vector;
using namespace Eigen;

//此函数返回L, D
template<typename _T>
vector<Mat> CholeskyLDLT(const Mat& M/*被分解矩阵*/)
{
	Mat L = M;
	Mat V(L.rows(), 1);
	for (int j = 0; j < L.rows(); ++j)
	{
		for (int i = 0; i < j; ++i)
			V(i) = L(j, i)*L(i, i);
		L(j, j) -= (L.block(j, 0, 1, j)*V.block(0, 0, j, 1))(0, 0);
		L.block(j+1, j, L.rows()-j-1, 1) -= L.block(j+1, 0, L.rows()-j-1, j)*V.block(0, 0, j, 1);
		L.block(j+1, j, L.rows()-j-1, 1) /= L(j, j);
	}
	Mat D = Mat::Zero(L.rows(), L.cols());
	for (int i = 0; i < L.rows(); ++i)
	{
		D(i, i) = L(i, i);
		L(i, i) = 1.0;
	}
	for (int i = 0; i < L.rows(); ++i)
		for (int j = 0; j < L.cols(); ++j)
			if (j > i) {L(i, j) = 0.0;}
	return vector<Mat>({L, D});
}

#endif


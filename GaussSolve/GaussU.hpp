#ifndef GAUSSU
#define GAUSSU
#define Mat Matrix<_T, Dynamic, Dynamic>

#include <Eigen/Core>

using namespace Eigen;

//此函数用于解上三角方程组
template<typename _T>
Mat GaussU(const Mat& U/*被解矩阵*/, const Mat& B/*被解向量*/)
{
	Mat b = B;
	for (int j = B.rows()-1; j > 0; --j)
	{
		b.row(j) /= U(j, j);
		Mat Ub = U.block(0, j, j, 1);
		Mat M(j, b.cols());
		for (int i = 0; i < M.cols(); ++i)
			M.col(i) = b(j, i)*Ub;
		b.topRows(j) -= M;
	}
	b.row(0) /= U(0, 0);
	return b;
}

#endif


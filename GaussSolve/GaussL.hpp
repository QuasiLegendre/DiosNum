#ifndef GAUSSL
#define GAUSSL
#define Mat Matrix<_T, Dynamic, Dynamic>

#include <Eigen/Core>

using namespace Eigen;

//此函数用于解下三角方程组
template<typename _T>
Mat GaussL(const Mat& L/*被解矩阵*/, const Mat& B/*被解向量*/)
{
	Mat b = B;
	for (int j = 0; j < B.rows(); ++j)
	{
		b.row(j) /= L(j, j);
		Mat Lb = L.block(j+1, j, B.rows()-j-1, 1);
		Mat M(B.rows()-j-1, b.cols());
		for (int i = 0; i < M.cols(); ++i)
			M.col(i) = b(j, i)*Lb;
		b.bottomRows(B.rows()-j-1) -= M;
	}
	b.row(B.rows()-1) /= L(B.rows()-1, B.rows()-1);
	return b;
}

#endif


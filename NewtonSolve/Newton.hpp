#ifndef TYPE_DEDUCE
#define TYPE_DEDUCE

template<typename _T>//帮助函数识别匿名函数类型
struct TypeDeduce
{
	using type = _T;
};

#endif
//Newton法解方程组
#ifndef NEWTON
#define NEWTON
#define Mat Matrix<_T, Dynamic, Dynamic>//定义简写方法
#include <functional>//import std::function
#include <Eigen/Core>
#include <Eigen/LU>//MatrixBase::inverse
using std::function;
using std::abs;
using namespace Eigen;

//此函数返回解
template<typename _T>
const Mat Newton(typename TypeDeduce<function<Mat(Mat)> >::type f/*函数*/, typename TypeDeduce<function<Mat(Mat)> >::type J/*f的Jacobi矩阵*/, const Mat& X_0/*初始值*/, const _T E/*容差*/, const unsigned long N/*最大迭代次数*/)
{
	Mat x = X_0;
	unsigned long n = 0;
	_T Max;
	do
	{
		Mat xx = J(x).inverse()*f(x);//求被减去向量
		x -= xx;//一步
		Max = abs(xx(0, 0));//求最大误差
		for (int i = 1; i < xx.rows(); ++i)
		{
			if (abs(xx(i, 0)) > Max) Max = abs(xx(i, 0));
		}
		++n;
	}while(Max >= E && n < N);
	return x;
}

#endif


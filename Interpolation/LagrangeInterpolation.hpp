#ifndef TYPE_DEDUCE
#define TYPE_DEDUCE

template<typename _T>//帮助函数识别匿名函数类型
struct TypeDeduce
{
	using type = _T;
};

#endif
//拉格朗日插值法
#ifndef LAGRANGE_INTERPOLATION
#define LAGRANGE_INTERPOLATION

#include <functional>//import std::function
#include <vector>//import std::vector
#include <numeric>//import std::accumulate
#include <algorithm>//import std::generate, std::transform

using std::function;
using std::vector;
using std::accumulate;
using std::generate;
using std::transform;

//此函数返回一个内含已经计算完毕系数的vector的匿名函数，每次在main中调用时仅需计算return的匿名函数内部的部分
template<typename _T>
function<_T(_T)> LagrangeInterpolation(typename TypeDeduce<function<_T(_T)> >::type f/*被插值函数*/, int n/*插值多项式项数*/, const _T* range/*插值范围*/)
{
	const _T MINVAL = range[0];//最小值
	const _T MAXVAL = range[1];//最大值
	const _T DX = (MAXVAL-MINVAL)/(_T(n));//间隔
	vector<_T> xVec(n+1), fVec(n+1), omegaVec(n+1);
	//分别存储插值点，被插函数值，与联乘值omega（参考教科书）的vector
	_T tmp = MINVAL - DX;
	
	generate(xVec.begin(), xVec.end(), [&]{return tmp += DX;});
	//生成等间距插值点，具体请看http://www.cplusplus.com/reference/algorithm/generate/?kw=generate
	xVec[0] = MINVAL;
	xVec[n] = MAXVAL;
	transform(xVec.begin(), xVec.end(), fVec.begin(), [f](_T x){return f(x);});
	//生成插值函数值，具体请看http://www.cplusplus.com/reference/algorithm/transform/?kw=transform
	vector<vector <_T> > omegaMat(n+1, vector<_T> (n+1));//临时变量，存储x[i] - x[j]
	for (int i = 0; i <= n; ++i)
		for (int j = 0; j <= n; ++j)
			omegaMat[i][j] = xVec[i] - xVec[j];
	for (int i = 0; i <= n; ++i)
		omegaMat[i][i] = 1.0;
	for (int i = 0; i <= n; ++i)
		omegaVec[i] = accumulate(omegaMat[i].begin(), omegaMat[i].end(), 1.0, [](_T x, _T y){return x*y;});
		//行内乘法计算omega，具体请看http://www.cplusplus.com/reference/numeric/accumulate/?kw=accumulate
	return [xVec, fVec, omegaVec, n](_T x)//返回插值匿名函数
	{
		_T numerator = 0.0;
		_T denominator = 0.0;
		for (int i = 0; i <= n; ++i)
		{
			numerator += fVec[i]/((x-xVec[i])*omegaVec[i]);
			denominator += 1.0/((x-xVec[i])*omegaVec[i]);
		}
		return numerator/denominator;
	};
}

#endif


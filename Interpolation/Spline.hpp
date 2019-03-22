#ifndef TYPE_DEDUCE
#define TYPE_DEDUCE
template<typename _T>//帮助识别匿名函数
struct TypeDeduce
{
	using type = _T;
};
#endif

#ifndef SPLINE
#define SPLINE

#include <functional>
#include <vector>
#include <algorithm>
//#include <exception>

using std::function;
using std::vector;
using std::generate;
using std::transform;

template<typename _T>
function<_T(_T)> Spline(typename TypeDeduce<function<_T(_T)> >::type f/*被插函数*/, int n/*插值区间数*/, const _T* range/*插值范围*/, const _T* EndpointDerivative/*边缘导数*/)
{
	const _T MINVAL = range[0];
	const _T MAXVAL = range[1];
	const _T DX = (MAXVAL - MINVAL)/(_T(n));
	vector<_T> xVec(n+1), fVec(n+1), alphaVec(n+1), betaVec(n+1), mVec(n+1);
	_T tmp = MINVAL - DX;
	mVec[0] = EndpointDerivative[0];//f'(x_0)
	mVec[n] = EndpointDerivative[1];//f'(x_n)

	generate(xVec.begin(), xVec.end(), [&]{return tmp += DX;});
	xVec[0] = MINVAL;
	xVec[n] = MAXVAL;
	transform(xVec.begin(), xVec.end(), fVec.begin(), [f](_T x){return f(x);});
	alphaVec[0] = 0.5;//alpha_0
	alphaVec[n] = 0.5;//alpha_n
	//以上部分类似LagrangeInterpolation.hpp
	for (int i = 1; i < n; ++i)
		alphaVec[i] = (xVec[i] - xVec[i-1])/(xVec[i+1] - xVec[i-1]);
	for (int i = 1; i < n; ++i)
		betaVec[i] = 3*((1.0 - alphaVec[i])*(fVec[i] - fVec[i-1])/(xVec[i] - xVec[i-1]) + alphaVec[i]*(fVec[i+1] - fVec[i])/(xVec[i+1] - xVec[i]));
	//三对角矩阵Crout分解开始
	vector<_T> L(n+1), L1(n+1), U1(n+1), z(n+1);//L[i] = l_{ii}, L1[i] = l_{i,i+1}, U1[i] = u_{i,i+1}, z[i] = z_i, from Numerical Analysis 10th edition
	L[1] = 2.0;
	U1[1] = alphaVec[1]/L[1];
	z[1] = (betaVec[1] - (1.0 - alphaVec[1])*mVec[0])/L[1];
	for (int i = 2; i <= n-2; ++i)
	{
		L1[i] = 1.0 - alphaVec[i];
		L[i] = 2.0 - L1[i]*U1[i-1];
		U1[i] = alphaVec[i]/L[i];
		z[i] = (betaVec[i] - L1[i]*z[i-1])/L[i];
	}
	L1[n-1] = 1.0 - alphaVec[n-1];
	L[n-1] = 2.0 - L1[n-1]*U1[n-2];
	z[n-1] = (betaVec[n-1] - alphaVec[n-1]*mVec[n] - L1[n-1]*z[n-2])/L[n-1];
	//以上LU分解
	mVec[n-1] = z[n-1];
	for (int i = n-2; i >= 1; --i)
		mVec[i] = z[i] - U1[i]*mVec[i+1];
	//Crout分解结束
	
	return [xVec, fVec, mVec](_T x)
	{
		//if (x < xVec.front() || x > xVec.back())
		//	throw domain_error();
		int it = find_if(xVec.begin(), xVec.end(), [x](_T y){return x < y;}) - xVec.begin() - 1;
		return (1.0 + 2*(x - xVec[it])/(xVec[it+1] - xVec[it]))*((x - xVec[it+1])/(xVec[it] - xVec[it+1]))*((x - xVec[it+1])/(xVec[it] - xVec[it+1]))*fVec[it]
		+  (1.0 + 2*(x - xVec[it+1])/(xVec[it] - xVec[it+1]))*((x - xVec[it])/(xVec[it+1] - xVec[it]))*((x - xVec[it])/(xVec[it+1] - xVec[it]))*fVec[it+1]
		+  (x - xVec[it])*((x - xVec[it+1])/(xVec[it] - xVec[it+1]))*((x - xVec[it+1])/(xVec[it] - xVec[it+1]))*mVec[it]
		+  (x - xVec[it+1])*((x - xVec[it])/(xVec[it+1] - xVec[it]))*((x - xVec[it])/(xVec[it+1] - xVec[it]))*mVec[it+1];//参考教科书公式
	};
}
#endif 

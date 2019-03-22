#ifndef TYPE_DEDUCE
#define TYPE_DEDUCE
template<typename _T>//帮助识别匿名函数
struct TypeDeduce
{
	using type = _T;
};
#endif

#ifndef PIECEWISE_LINEAR_INTERPOLATION
#define PIECEWISE_LINEAR_INTERPOLATION

#include <functional>//import std::function
#include <vector>//import std::vector
#include <algorithm>//import std::generate, std::tranform, std::find_if
//#include <exception>

using std::function;
using std::vector;
using std::generate;
using std::transform;
using std::find_if;
//using std::domain_error;

template<typename _T>
function<_T(_T)> PiecewiseLinearInterpolation(typename TypeDeduce<function<_T(_T)> >::type f, int n, const _T* range)
{
	const _T MINVAL = range[0];
	const _T MAXVAL = range[1];
	const _T DX = (MAXVAL - MINVAL)/(_T(n));
	vector<_T> xVec(n+1), fVec(n+1);
	_T tmp = MINVAL - DX;
	
	generate(xVec.begin(), xVec.end(), [&]{return tmp += DX;});
	xVec[0] = MINVAL;
	xVec[n] = MAXVAL;
	transform(xVec.begin(), xVec.end(), fVec.begin(), [f](_T x){return f(x);});
	//以上同LagrangeInterpolation.hpp
	return [xVec, fVec](_T x)
	{
		//if (x < xVec.front() || x > xVec.back())
		//	throw domain_error();
		int it = find_if(xVec.begin(), xVec.end(), [x](_T y){return x < y;}) - xVec.begin() - 1;//寻找插值函数所需域, find_if请看http://www.cplusplus.com/reference/algorithm/find_if/?kw=find_if
		return fVec[it]*(x - xVec[it+1])/(xVec[it] - xVec[it+1]) + fVec[it+1]*(x - xVec[it])/(xVec[it+1] - xVec[it]);//计算结果
	};
}
#endif 

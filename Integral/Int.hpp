#ifndef TYPE_DEDUCE
#define TYPE_DEDUCE

template<typename _T>//帮助函数识别匿名函数类型
struct TypeDeduce
{
	using type = _T;
};

#endif
//Romberg积分
#ifndef ROMBERG_INTEGRAL
#define ROMBERG_INTEGRAL

#include <functional>//import std::function
#include <vector>//import std::vector

using std::function;
using std::vector;


//此函数返回积分值常量，每次在main中调用时仅需计算return的匿名函数内部的部分
template<typename _T>
const _T Int(typename TypeDeduce<function<_T(_T)> >::type f/*被积函数*/,const _T E/*容差*/, const _T* RANGE/*积分范围*/)
{
	const _T MINVAL = RANGE[0];//积分下限
	const _T MAXVAL = RANGE[1];//积分上限
	
	_T h = (MAXVAL-MINVAL)/(_T(2.0));//间隔	
	int k=1, n=1;
	vector <_T>T(k), Fours(k);//T:存放T;Fours:存放4的幂;
	T[0] = h*(f(MINVAL)+f(MAXVAL));//初始化
	Fours[0] = 1.0;//初始化
	_T T_0, F;
	
	do
	{
		T_0 = T[0];//存放上一个T[0]用于比较
		F = 0;
		for (int i = 1; i <= n; ++i)
			F += f(MINVAL + (2*i - 1)*h);//计算F
		T.push_back(T[0]/2.0 + h*F);//计算并存放T^(k)
		Fours.push_back(Fours[k-1]*4.0);//计算并存放4^k
		for (int m = 1; m <= k; ++m)
			T[k-m] = (Fours[m]*T[k-m+1] - T[k-m])/(Fours[m] - 1.0);//此处计算别的T，其中无用的T被替换
		h /= 2; n *= 2; ++k;
	}
	while((T[0]>T_0?T[0]-T_0:T_0-T[0]) >= E);//绝对值若大于等于E则继续，否则停止

	return T[0];
}

#endif


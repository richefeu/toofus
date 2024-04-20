#ifndef NR3_THINGS_HPP
#define NR3_THINGS_HPP

template<class T>
inline const T &MAX(const T &a, const T &b)
{
	return b > a ? (b) : (a);
}
// may be replaced by std::max() of cmath

template<class T>
inline void SWAP(T &a, T &b)
{
	T dum = a;
	a = b;
	b = dum;
}
// may be replaced by std::swap(a,b) of cmath

template<class T>
inline T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
// may be replaced by std::copysign(a,b) of cmath

template<class T> inline T SQR(const T a)
{
	return a * a;
}
// may be replaced by Mth::sqr(a) in Mth.hpp

#endif /* end of include guard: NR3_THINGS_HPP */

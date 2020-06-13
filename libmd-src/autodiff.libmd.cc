#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

dual::dual()
{
    //!
    //! Default constructor for dual, does nothing.
    //!
}

dual::dual(ldf f,ldf fx)
{
    //!
    //! Constructor for dual, sets value (<tt>x</tt>) to <tt>f</tt> and derivative (<tt>dx</tt>) to <tt>fx</tt> (default: <tt>fx</tt>=0).
    //!
    x=f;
    dx=fx;
}


// Assignment

template<class X> dual dual::operator=(X a)
{
    //!
    //! Sets value (<tt>x</tt>) to <tt>a</tt> and derivative (<tt>dx</tt>) to 0 and returns <tt>*this</tt>.
    //!
    return *this=dual(a);
}

template<class X> dual::operator X()
{
    //!
    //! Cast operator, returns cast of value (<tt>x</tt>).
    //!
    return x;
}


// Comparisons

bool operator==(dual F,dual G)
{
    return F.x==G.x;
}

bool operator!=(dual F,dual G)
{
    return F.x!=G.x;
}

bool operator<=(dual F,dual G)
{
    return F.x<=G.x;
}

bool operator>=(dual F,dual G)
{
    return F.x>=G.x;
}

bool operator<(dual F,dual G)
{
    return F.x<G.x;
}

bool operator>(dual F,dual G)
{
    return F.x>G.x;
}

template<class X> bool operator==(dual &F,X a)
{
    return F.x==a;
}

template<class X> bool operator!=(dual &F,X a)
{
    return F.x!=a;
}

template<class X> bool operator<=(dual &F,X a)
{
    return F.x<=a;
}

template<class X> bool operator>=(dual &F, X a)
{
    return F.x>=a;
}

template<class X> bool operator<(dual F, X a)
{
    return F.x<a;
}

template<class X> bool operator>(dual F, X a)
{
    return F.x>a;
}

template<class X> bool operator==(X a,dual &F)
{
    return a==F.x;
}

template<class X> bool operator!=(X a,dual &F)
{
    return a!=F.x;
}

template<class X> bool operator<=(X a,dual &F)
{
    return a<=F.x;
}

template<class X> bool operator>=(X a,dual &F)
{
    return a>=F.x;
}

template<class X> bool operator<(X a,dual F)
{
    return a<F.x;
}

template<class X> bool operator>(X a,dual F)
{
    return a>F.x;
}


// Standard operations with other duals

dual operator-(dual F)
{
    return dual(-F.x,-F.dx);
}

dual operator+(dual F,dual G)
{
    return dual(F.x+G.x, F.dx+G.dx);
}

dual operator-(dual F,dual G)
{
    return dual(F.x-G.x, F.dx-G.dx);
}

dual operator*(dual F,dual G)
{
    return dual(F.x*G.x,F.dx*G.x+F.x*G.dx);
}

dual operator/(dual F,dual G)
{
    return dual(F.x/G.x,(F.dx*G.x-F.x*G.dx)/std::pow(G.x,2));
}

dual operator+=(dual& F,dual G)
{
    F.x+=G.x;
    F.dx+=G.dx;
    return F;
}

dual operator-=(dual& F,dual G)
{
    F.x-=G.x;
    F.dx-=G.dx;
    return F;
}

dual operator*=(dual& F,dual G)
{
    return F=F*G;
}

dual operator/=(dual& F,dual G)
{
    return F=F/G;
}


// Standard operations with constants

template<class X> dual operator+(dual F,X a)
{
    return dual(F.x+a,F.dx);
}

template<class X> dual operator-(dual F,X a)
{
    return dual(F.x-a,F.dx);
}

template<class X> dual operator*(dual F,X a)
{
    return dual(F.x*a,F.dx*a);
}

template<class X> dual operator/(dual F,X a)
{
    return dual(F.x/a,F.dx/a);
}

template<class X> dual operator+(X a,dual F)
{
    return F+a;
}

template<class X> dual operator-(X a,dual F)
{
    return (-F)+a;
}

template<class X> dual operator*(X a,dual F)
{
    return F*a;
}

template<class X> dual operator/(X a,dual F)
{
    return dual(a/F.x,-a*F.dx/std::pow(F.x,2));
}

template<class X> dual operator+=(dual& F,X a)
{
    F.x+=a;
    return F;
}

template<class X> dual operator-=(dual& F,X a)
{
    F.x-=a;
    return F;
}

template<class X> dual operator*=(dual& F,X a)
{
    F.x*=a;
    F.dx*=a;
    return F;
}

template<class X> dual operator/=(dual& F,X a)
{
    F.x/=a;
    F.dx/=a;
    return F;
}


// Standard functions

dual sqrt(dual F)
{
    const ldf sqr=std::sqrt(F.x);
    return dual(sqr,F.dx/2.0/sqr);
}

template<class X> dual pow(dual F,X n)
{
    return dual(std::pow(F.x,n),F.dx*n*std::pow(F.x,n-1));
}

dual exp(dual F)
{
    const ldf ex=std::exp(F.x);
    return dual(ex,F.dx*ex);
}

template<class X> dual pow(X a, dual G)
{
    const ldf po=std::pow(a,G.x);
    return dual(po,G.dx*log(a)*po);
}

dual pow(dual F,dual G)
{
    const ldf po=std::pow(F.x,G.x);
    return dual(po,(G.dx*log(F.x)+F.dx*G.x/F.x)*po);
}

dual log(dual F)
{
    return dual(std::log(F.x),F.dx/F.x);
}

dual sin(dual F)
{
    return dual(std::sin(F.x),F.dx*std::cos(F.x));
}

dual cos(dual F)
{
    return dual(std::cos(F.x),-F.dx*std::sin(F.x));
}

dual tan(dual F)
{
    const ldf ta=std::tan(F.x);
    return dual(ta,F.dx*(1.0+std::pow(ta,2)));
}

dual sinh(dual F)
{
    return dual(std::sinh(F.x),F.dx*std::cosh(F.x));
}

dual cosh(dual F)
{
    return dual(std::cosh(F.x),F.dx*std::sinh(F.x));
}

dual tanh(dual F)
{
    const ldf ta=std::tanh(F.x);
    return dual(ta,F.dx*(1.0-std::pow(ta,2)));
}

dual asin(dual F)
{
    return dual(std::asin(F.x),F.dx/std::sqrt(1.0-std::pow(F.x,2)));
}

dual acos(dual F)
{
    return dual(std::acos(F.x),-F.dx/std::sqrt(1.0-std::pow(F.x,2)));
}

dual atan(dual F)
{
    return dual(std::atan(F.x),F.dx/(std::pow(F.x,2)+1.0));
}

dual asinh(dual F)
{
    return dual(std::asinh(F.x),F.dx/std::sqrt(1.0+std::pow(F.x,2)));
}

dual acosh(dual F)
{
    return dual(std::acosh(F.x),-F.dx/std::sqrt(std::pow(F.x,2)-1.0));
}

dual atanh(dual F)
{
    return dual(std::atanh(F.x),F.dx/(1.0-std::pow(F.x,2)));
}

dual abs(dual F)
{
    return F.x<0.0?-F:F;
}

dual heaviside(dual F)
{
    return dual(F.x<0.0?0.0:1.0, F.x==0.0?std::numeric_limits<ldf>::infinity():0.0);
}

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

dual dual::operator=(dual G)
{
    //!
    //! Copies <tt>G</tt> to <tt>this</tt> and returns <tt>*this</tt>.
    //!
    x=G.x;
    dx=G.dx;
    return *this;
}

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
    return dual(F.x/G.x,(F.dx*G.x-F.x*G.dx)/pow(G.x,2));
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
    return dual(a/F.x,-a*F.dx/pow(F.x,2));
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
    return dual(sqrt(F.x),F.dx/2.0/sqrt(F.x));
}

template<class X> dual pow(dual F,X n)
{
    return dual(pow(F.x, n),F.dx *n*pow(F.x,n-1));
}

dual exp(dual F)
{
    return dual(exp(F.x),F.dx*exp(F.x));
}

template<class X> dual pow(X a, dual G)
{
    return dual(pow(a, G.x),G.dx*log(a)*pow(a, G.x));
}

dual pow(dual F,dual G)
{
    return dual(pow(F.x, G.x),(G.dx*log(F.x)+F.dx*G.x/F.x)*pow(F.x, G.x));
}

dual log(dual F)
{
    return dual(log(F.x),F.dx/F.x);
}

dual sin(dual F)
{
    return dual(sin(F.x),F.dx*cos(F.x));
}

dual cos(dual F)
{
    return dual(cos(F.x),-F.dx*sin(F.x));
}

dual tan(dual F)
{
    return dual(tan(F.x),F.dx/pow(cos(F.x),2));
}

dual sinh(dual F)
{
    return dual(sinh(F.x),F.dx*cosh(F.x));
}

dual cosh(dual F)
{
    return dual(cosh(F.x),F.dx*sinh(F.x));
}

dual tanh(dual F)
{
    return dual(tanh(F.x),F.dx/pow(cosh(F.x),2));
}

dual asin(dual F)
{
    return dual(asin(F.x),F.dx/sqrt(1.0-pow(F.x,2)));
}

dual acos(dual F)
{
    return dual(acos(F.x),-F.dx/sqrt(1.0-pow(F.x,2)));
}

dual atan(dual F)
{
    return dual(atan(F.x),F.dx/(pow(F.x,2)+1.0));
}

dual asinh(dual F)
{
    return dual(asinh(F.x),F.dx/sqrt(1.0+pow(F.x,2)));
}

dual acosh(dual F)
{
    return dual(acosh(F.x),-F.dx/sqrt(pow(F.x,2)-1.0));
}

dual atanh(dual F)
{
    return dual(atanh(F.x),F.dx/(1.0-pow(F.x,2)));
}

dual fabs(dual F)
{   
    return F.x<0.0?-F:F;
}

dual heaviside(dual F)
{
    return dual(F.x<0.0?0.0:1.0, F.x==0.0?numeric_limits<ldf>::infinity():0.0);
}



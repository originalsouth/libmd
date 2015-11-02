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
    using namespace std;
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
    const ldf sqr=sqrt(F.x);
    return dual(sqr,F.dx/2.0/sqr);
}

template<class X> dual pow(dual F,X n)
{
    using namespace std;
    return dual(pow(F.x,n),F.dx*n*pow(F.x,n-1));
}

dual exp(dual F)
{
    using namespace std;
    const ldf ex=exp(F.x);
    return dual(ex,F.dx*ex);
}

template<class X> dual pow(X a, dual G)
{
    using namespace std;
    const ldf po=pow(a,G.x);
    return dual(po,G.dx*log(a)*po);
}

dual pow(dual F,dual G)
{
    using namespace std;
    const ldf po=pow(F.x,G.x);
    return dual(po,(G.dx*log(F.x)+F.dx*G.x/F.x)*po);
}

dual log(dual F)
{
    using namespace std;
    return dual(log(F.x),F.dx/F.x);
}

dual sin(dual F)
{
    using namespace std;
    return dual(sin(F.x),F.dx*cos(F.x));
}

dual cos(dual F)
{
    using namespace std;
    return dual(cos(F.x),-F.dx*sin(F.x));
}

dual tan(dual F)
{
    using namespace std;
    const ldf ta=tan(F.x);
    return dual(ta,F.dx*(1.0+pow(ta,2)));
}

dual sinh(dual F)
{
    using namespace std;
    return dual(sinh(F.x),F.dx*cosh(F.x));
}

dual cosh(dual F)
{
    using namespace std;
    return dual(cosh(F.x),F.dx*sinh(F.x));
}

dual tanh(dual F)
{
    using namespace std;
    const ldf ta=tanh(F.x);
    return dual(ta,F.dx*(1-pow(ta,2)));
}

dual asin(dual F)
{
    using namespace std;
    return dual(asin(F.x),F.dx/sqrt(1.0-pow(F.x,2)));
}

dual acos(dual F)
{
    using namespace std;
    return dual(acos(F.x),-F.dx/sqrt(1.0-pow(F.x,2)));
}

dual atan(dual F)
{
    using namespace std;
    return dual(atan(F.x),F.dx/(pow(F.x,2)+1.0));
}

dual asinh(dual F)
{
    using namespace std;
    return dual(asinh(F.x),F.dx/sqrt(1.0+pow(F.x,2)));
}

dual acosh(dual F)
{
    using namespace std;
    return dual(acosh(F.x),-F.dx/sqrt(pow(F.x,2)-1.0));
}

dual atanh(dual F)
{
    using namespace std;
    return dual(atanh(F.x),F.dx/(1.0-pow(F.x,2)));
}

dual abs(dual F)
{
    return F.x<0.0?-F:F;
}

dual heaviside(dual F)
{
    return dual(F.x<0.0?0.0:1.0, F.x==0.0?std::numeric_limits<ldf>::infinity():0.0);
}

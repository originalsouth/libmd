#ifndef libmd_h
#include "../libmd.h"
#endif

// Single differentiation
dual::dual()
{
    dx=1.0;
}

dual::dual(ldf f,ldf fx)
{
    x=f;
    dx=fx;
}

dual dual::operator=(dual y)
{
    x=y.x;
    dx=y.dx;
    return y;
}

void dual::operator+=(dual y)
{
    x+=y.x;
    dx+=y.dx;
}

void dual::operator-=(dual y)
{
    x-=y.x;
    dx-=y.dx;
}

template<class X> X dual::operator=(X y)
{
    x=y;
    return x;
}

template<class X> void dual::operator+=(X y)
{
    x+=y;
}

template<class X> void dual::operator-=(X y)
{
    x-=y;
}

template<class X> void dual::operator*=(X y)
{
    x*=y;
}

template<class X> void dual::operator/=(X y)
{
    x/=y;
}

template<class X> bool dual::operator==(X y)
{
    return x==y;
}

template<class X> bool dual::operator<=(X y)
{
    return x<=y;
}

template<class X> bool dual::operator>=(X y)
{
    return x>=y;
}

template<class X> bool dual::operator<(X y)
{
    return x<y;
}

template<class X> bool dual::operator>(X y)
{
    return x>y;
}

dual operator+(dual F,dual G)
{
    return dual(F.x+G.x, F.dx+G.dx);
}

dual operator-(dual F,dual G)
{
    return dual(F.x-G.x, F.dx-G.dx);
}

dual operator-(dual F)
{
    return dual(-F.x,-F.dx);
}

dual operator*(dual F,dual G)
{
    return dual(F.x*G.x,F.dx*G.x+F.x*G.dx);
}

dual operator/(dual F,dual G)
{
    return dual(F.x/G.x,(F.dx*G.x-F.x*G.dx)/pow(G.x,2));
}

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
    return sin(F)/cos(F);
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
    return sinh(F)/cosh(F);
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

dual fabs(dual F)
{   
    return dual(fabs(F.x), F.dx*((F.x < 0) ? -1 : 1)); // strictly, the differential of |f(x)| is undefined at f(x)=0. this function however returns f'(x).
}

#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> duals<dim>::duals()
{
    //!
    //! Default constructor for duals, does nothing.
    //!
}

template<ui dim> duals<dim>::duals(ldf a)
{
    //!
    //! Constructor for duals, sets value (<tt>x</tt>) to <tt>a</tt> and derivatives to 0.
    //!
    x=a;
    memset(dx,0,sizeof(dx));
    memset(dxdy,0,sizeof(dxdy));
}

template<ui dim> duals<dim>::duals(ldf a,ui i)
{
    //!
    //! Constructor for duals, sets value (<tt>x</tt>) to <tt>a</tt> and derivatives to 0,
    //! with the exception of the derivative with respect to <tt>x[i]</tt>, which is set to 1.
    //!
    x=a;
    memset(dx,0,sizeof(dx));
    memset(dxdy,0,sizeof(dxdy));
    dx[i]=1.0;
}


// Assignment

template<ui dim> duals<dim> duals<dim>::operator=(duals<dim> G)
{
    //!
    //! Copies <tt>G</tt> to <tt>this</tt> and returns <tt>*this</tt>.
    //!
    x=G.x;
    memcpy(dx,G.dx,sizeof(dx));
    memcpy(dxdy,G.dxdy,sizeof(dxdy));
    return *this;
}

template<ui dim> template<class X> duals<dim> duals<dim>::operator=(X a)
{
    //!
    //! Sets value (<tt>x</tt>) to <tt>a</tt> and derivatives to 0 and returns <tt>*this</tt>.
    //!
    return *this=duals<dim>(a);
}

template<ui dim> template<class X> duals<dim>::operator X()
{
    //!
    //! Cast operator, returns cast of value (<tt>x</tt>).
    //!
    return x;
}


// Comparisons

template<ui dim> bool operator==(duals<dim> F,duals<dim> G)
{
    return F.x==G.x;
}

template<ui dim> bool operator!=(duals<dim> F,duals<dim> G)
{
    return F.x!=G.x;
}

template<ui dim> bool operator<=(duals<dim> F,duals<dim> G)
{
    return F.x<=G.x;
}

template<ui dim> bool operator>=(duals<dim> F,duals<dim> G)
{
    return F.x>=G.x;
}

template<ui dim> bool operator<(duals<dim> F,duals<dim> G)
{
    return F.x<G.x;
}

template<ui dim> bool operator>(duals<dim> F,duals<dim> G)
{
    return F.x>G.x;
}

template<ui dim, class X> bool operator==(duals<dim> F,X a)
{
    return F.x==a;
}

template<ui dim, class X> bool operator!=(duals<dim> F,X a)
{
    return F.x!=a;
}

template<ui dim, class X> bool operator<=(duals<dim> F,X a)
{
    return F.x<=a;
}

template<ui dim, class X> bool operator>=(duals<dim> F, X a)
{
    return F.x>=a;
}

template<ui dim, class X> bool operator<(duals<dim> F, X a)
{
    return F.x<a;
}

template<ui dim, class X> bool operator>(duals<dim> F, X a)
{
    return F.x>a;
}

template<ui dim, class X> bool operator==(X a,duals<dim> F)
{
    return a==F.x;
}

template<ui dim, class X> bool operator!=(X a,duals<dim> F)
{
    return a!=F.x;
}

template<ui dim, class X> bool operator<=(X a,duals<dim> F)
{
    return a<=F.x;
}

template<ui dim, class X> bool operator>=(X a,duals<dim> F)
{
    return a>=F.x;
}

template<ui dim, class X> bool operator<(X a,duals<dim> F)
{
    return a<F.x;
}

template<ui dim, class X> bool operator>(X a,duals<dim> F)
{
    return a>F.x;
}


// Standard operations with other duals's

template<ui dim> duals<dim> operator-(duals<dim> F)
{
    return F*(-1.0);
}

template<ui dim> duals<dim> operator+(duals<dim> F,duals<dim> G)
{
    duals<dim> H(F.x+G.x);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]+G.dx[i];
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]+G.dxdy[i][j];
    }
    return H;
}

template<ui dim> duals<dim> operator-(duals<dim> F,duals<dim> G)
{
    duals<dim> H(F.x-G.x);
    for (ui i=0;i<dim; i++)
    {
        H.dx[i]=F.dx[i]-G.dx[i];
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]-G.dxdy[i][j];
    }
    return H;
}

template<ui dim> duals<dim> operator* (duals<dim> F,duals<dim> G)
{
    duals<dim> H(F.x*G.x);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*G.x+F.x*G.dx[i];
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.x*G.dxdy[i][j]+F.dx[i]*G.dx[j]+F.dx[j]*G.dx[i]+F.dxdy[i][j]*G.x;
    }
    return H;
}

template<ui dim> duals<dim> operator/(duals<dim> F,duals<dim> G)
{
    duals<dim>H(F.x/G.x);
    const ldf g2=pow(G.x,2),g3=pow(G.x,3);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=(F.dx[i]*G.x-F.x*G.dx[i])/g2;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(2.0*F.x*G.dx[i]*G.dx[j]-G.x*(F.x*G.dxdy[i][j]+F.dx[i]*G.dx[j]+F.dx[j]*G.dx[i]-G.x*F.dxdy[i][j]))/g3;
    }
    return H;
}

template<ui dim> duals<dim> operator+=(duals<dim>& F,duals<dim> G)
{
    return F=F+G;
}

template<ui dim> duals<dim> operator-=(duals<dim>& F,duals<dim> G)
{
    return F=F-G;
}

template<ui dim> duals<dim> operator*=(duals<dim>& F,duals<dim> G)
{
    return F=F*G;
}

template<ui dim> duals<dim> operator/=(duals<dim>& F,duals<dim> G)
{
    return F=F/G;
}


// Standard operations with constants

template<ui dim, class X> duals<dim> operator+(duals<dim> F,X a)
{
    F.x+=a;
    return F;
}

template<ui dim, class X> duals<dim> operator-(duals<dim> F,X a)
{
    F.x-=a;
    return F;
}

template<ui dim, class X> duals<dim> operator*(duals<dim> F,X a)
{
    F.x*=a;
    for(ui i=0;i<dim;i++)
    {
        F.dx[i]*=a;
        for(ui j=0;j<dim;j++) F.dxdy[i][j]*=a;
    }
    return F;
}

template<ui dim, class X> duals<dim> operator/(duals<dim> F,X a)
{
    F.x/=a;
    for(ui i=0;i<dim;i++)
    {
        F.dx[i]/=a;
        for(ui j=0;j<dim;j++) F.dxdy[i][j]/=a;
    }
    return F;
}

template<ui dim, class X> duals<dim> operator+(X a,duals<dim> F)
{
    return F+a;
}

template<ui dim, class X> duals<dim> operator-(X a,duals<dim> F)
{
    return (-F)+a;
}

template<ui dim, class X> duals<dim> operator*(X a,duals<dim> F)
{
    return F*a;
}

template<ui dim, class X> duals<dim> operator/(X a,duals<dim> F)
{
    duals<dim> H(a/F.x);
    const ldf f2=pow(F.x,2),f3=pow(F.x,3);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=-a*F.dx[i]/f2;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=a*(2.0*F.dx[i]*F.dx[j]-F.x*F.dxdy[i][j])/f3;
    }
    return H;
}

template<ui dim, class X> duals<dim> operator+=(duals<dim>& F,X a)
{
    F.x+=a;
    return F;
}

template<ui dim, class X> duals<dim> operator-=(duals<dim>& F,X a)
{
    F.x-=a;
    return F;
}

template<ui dim, class X> duals<dim> operator*=(duals<dim>& F,X a)
{
    return F=F*a;
}

template<ui dim, class X> duals<dim> operator/=(duals<dim>& F,X a)
{
    return F=F/a;
}


// Standard functions

template<ui dim> duals<dim> sqrt(duals<dim> F)
{
    const ldf r=sqrt(F.x);
    duals<dim> H(r);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]/2.0/r;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]-F.dx[i]*F.dx[j]/2.0/F.x)/2.0/r;
    }
    return H;
}

template<ui dim, class X> duals<dim> pow(duals<dim> F,X n)
{
    duals<dim> H(pow(F.x,n));
    const ldf z1=pow(F.x,n-1),z2=pow(F.x,n-2);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*n*z1;
        for (ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=n*(F.dxdy[i][j]*F.x+(n-1)*F.dx[i]*F.dx[j])*z2;
    }
    return H;
}

template<ui dim> duals<dim> exp(duals<dim> F)
{
    ldf z=exp(F.x);
    duals<dim> H(z);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*z;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]+F.dx[i]*F.dx[j])*z;
    }
    return H;
}

template<ui dim, class X> duals<dim> pow(X a,duals<dim> F)
{
    const ldf z=pow(a,F.x),la=log(a);
    duals<dim> H(z);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*la*z;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]+F.dx[i]*F.dx[j]*la)*la*z;
    }
    return H;
}

template<ui dim> duals<dim> pow(duals<dim> F,duals<dim> G)
{
    const ldf z=pow(F.x,G.x),lf=log(F.x);
    duals<dim> H(z);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=(G.dx[i]*lf+F.dx[i]*G.x/F.x)*z;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=((G.dx[i]*lf+F.dx[i]*G.x/F.x)*(G.dx[j]*lf+F.dx[j]*G.x/F.x)+G.dxdy[i][j]*lf+(F.dx[i]*G.dx[j]+F.dx[j]*G.dx[i]+F.dxdy[i][j]*G.x-F.dx[i]*F.dx[j]*G.x/F.x)/F.x)*z;
    }
    return H;
}

template<ui dim> duals<dim> log(duals<dim> F)
{
    duals<dim> H(log(F.x));
    const ldf f2=pow(F.x,2);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]/F.x;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]*F.x-F.dx[i]*F.dx[j])/f2;
    }
    return H;
}

template<ui dim> duals<dim> sin(duals<dim> F)
{
    const ldf sf=sin(F.x),cf=cos(F.x);
    duals<dim> H(sf);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*cf;
        for (ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]*cf-F.dx[i]*F.dx[j]*sf;
    }
    return H;
}

template<ui dim> duals<dim> cos(duals<dim> F)
{
    const ldf cf=cos(F.x),sf=sin(F.x);
    duals<dim> H(cf);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=-F.dx[i]*sf;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=-F.dxdy[i][j]*sf-F.dx[i]*F.dx[j]*cf;
    }
    return H;
}

template<ui dim> duals<dim> tan(duals<dim> F)
{
    const ldf tf=tan(F.x);
    duals<dim> H(tf);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*(1.0+tf*tf);
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]+F.dx[i]*F.dx[j]*2.0*tf)*(1.0+tf*tf);
    }
    return H;
}

template<ui dim> duals<dim> sinh(duals<dim> F)
{
    const ldf sf=sinh(F.x),cf=cosh(F.x);
    duals<dim> H(sf);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*cf;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]*cf+F.dx[i]*F.dx[j]*sf;
    }
    return H;
}

template<ui dim> duals<dim> cosh (duals<dim> F)
{
    const ldf cf=cosh(F.x),sf=sinh(F.x);
    duals<dim> H(cf);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*sf;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]*sf+F.dx[i]*F.dx[j]*cf;
    }
    return H;
}

template<ui dim> duals<dim> tanh (duals<dim> F)
{
    const ldf tf=tanh(F.x);
    duals<dim> H(tf);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*(1.0-tf*tf);
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]-F.dx[i]*F.dx[j]*2.0*tf)*(1.0-tf*tf);
    }
    return H;
}

template<ui dim> duals<dim> asin(duals<dim> F)
{   duals<dim> H(asin(F.x));
    const ldf z1=pow(1.0-pow(F.x,2),-0.5),z2=pow(1.0-pow(F.x,2),-1.5);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*z1;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]*z1+F.dx[i]*F.dx[j]*F.x*z2;
    }
    return H;
}

template<ui dim> duals<dim> acos(duals<dim> F)
{
    duals<dim> H(acos(F.x));
    const ldf z1=pow(1.0-pow(F.x,2),-0.5),z2=pow(1.0-pow(F.x,2),-1.5);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=-F.dx[i]*z1;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=-F.dxdy[i][j]*z1-F.dx[i]*F.dx[j]*F.x*z2;
    }
    return H;
}

template<ui dim> duals<dim> atan (duals<dim> F)
{
    duals<dim> H(atan(F.x));
    const ldf z=1.0/(1.0+pow(F.x,2));
    for (ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*z;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]-F.dx[i]*F.dx[j]*2.0*F.x*z)*z;
    }
    return H;
}

template<ui dim> duals<dim> asinh(duals<dim> F)
{
    duals<dim> H(asinh(F.x));
    const ldf z1=pow(1.0+pow(F.x,2),-0.5),z2=-pow(1.0+pow(F.x,2),-1.5);
    for (ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*z1;
        for (ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]*z1+F.dx[i]*F.dx[j]*F.x*z2;
    }
    return H;
}

template<ui dim> duals<dim> acosh(duals<dim> F)
{
    duals<dim> H(acosh(F.x));
    const ldf z1=pow(pow(F.x,2)-1.0,-0.5),z2=-pow(pow(F.x,2)-1,-1.5);
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*z1;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=F.dxdy[i][j]*z1+F.dx[i]*F.dx[j]*F.x*z2;
    }
    return H;
}

template<ui dim> duals<dim> atanh(duals<dim> F)
{
    duals<dim> H(atanh(F.x));
    const ldf z = 1.0/(1.0-pow(F.x,2));
    for(ui i=0;i<dim;i++)
    {
        H.dx[i]=F.dx[i]*z;
        for(ui j=i;j<dim;j++) H.dxdy[i][j]=H.dxdy[j][i]=(F.dxdy[i][j]+F.dx[i]*F.dx[j]*2.0*F.x*z)*z;
    }
    return H;
}

template<ui dim> duals<dim> fabs(duals<dim> F)
{
    return F.x<0.0?-F:F;
}

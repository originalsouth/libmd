#ifndef libmd_h
#include "../libmd.h"
#endif

ldf kdelta(ui i,ui j)
{
    return (i==j)?1.0:0.0;
}

//Flat space
template<class X,ui dim> X FLATSPACE(X x[dim],vector<ldf> *param)
{
    (void) x;
    (void) param;
    return 0.0;
}

//Gaussian bump
template<class X,ui dim> X GAUSSIANBUMP(X x[dim],vector<ldf> *param)
{
    const ldf A=param->at(0);
    const ldf K=param->at(1);
    X retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return A*exp(-K*retval);
}

//Egg carton
template<class X,ui dim> X EGGCARTON(X *x,vector<ldf> *param)
{
    X retval=param->at(0);
    for(ui d=0;d<dim;d++) retval*=cos(x[d]*param->at(d+1));
    return retval;
}

//Mollifier
template<class X,ui dim> X MOLLIFIER(X *x,vector<ldf> *param)
{
    const ldf A=param->at(0);
    const ldf K=param->at(1);
    X retval=0.0,zero=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return (retval<K)?A*exp(-retval/(K-retval))/pow(retval-K,2):zero;
}

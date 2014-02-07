#ifndef libmd_h
#include "../libmd.h"
#endif

template<class X> X COULOMB(X r,vector<ldf> *parameters)
{
    const ldf q=parameters->at(0);
    return q/r;
}

template<class X> X YUKAWA(X r,vector<ldf> *parameters)
{
    const ldf b=parameters->at(0);
    const ldf k=parameters->at(1);
    return b/(r*exp(k*r));
}

template<class X> X HOOKIAN(X r,vector<ldf> *parameters)
{
    const ldf k=parameters->at(0);
    const ldf r0=parameters->at(1);
    return k/2.0*pow(r-r0,2);
}

template<class X> X LJ(X r,vector<ldf> *parameters)
{
    const ldf e=parameters->at(0);
    const ldf s=parameters->at(1);
    return 4.0*e*(pow(s/r,12)-pow(s/r,6));
}

template<class X> X MORSE(X r,vector<ldf> *parameters)
{
    const ldf d=parameters->at(0);
    const ldf a=parameters->at(1);
    const ldf re=parameters->at(3);
    return d*pow(1.0-exp(a*(re-r)),2);
}

template<class X> X FORCEDIPOLE(X r,vector<ldf> *parameters)
{
    const ldf f = parameters->at(0);
    return f*r;
}

template<class X> X HOOKEANFORCEDIPOLE(X r,vector<ldf> *parameters)
{
    const ldf k=parameters->at(0);
    const ldf r0=parameters->at(1);
    const ldf f = parameters->at(2);
    return  k/2.0*pow(r-r0,2)+f*r;
}

template<class X> X ANHARMONICSPRING(X r,vector<ldf> *parameters)
{
    const ldf k=parameters->at(0);
    const ldf r0=parameters->at(1);
    const ldf alpha=parameters->at(2);
    return (k/alpha)*pow(fabs(r-r0),alpha);
}

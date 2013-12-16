#ifndef libmd_h
#include "../libmd.h"
#endif

dual COULOMB(dual r,vector<ldf> *parameters)
{
    const ldf q=parameters->at(0);
    return q/r;
}

dual YUKAWA(dual r,vector<ldf> *parameters)
{
    const ldf b=parameters->at(0);
    const ldf k=parameters->at(1);
    return b/(r*exp(k*r));
}

dual HOOKIAN(dual r,vector<ldf> *parameters)
{
    const ldf k=parameters->at(0);
    const ldf r0=parameters->at(1);
    return k/2.0*pow(r-r0,2);
}

dual LJ(dual r,vector<ldf> *parameters)
{
    const ldf e=parameters->at(0);
    const ldf s=parameters->at(1);
    return 4.0*e*(pow(s/r,12)-pow(s/r,6));
}

dual MORSE(dual r,vector<ldf> *parameters)
{
    const ldf d=parameters->at(0);
    const ldf a=parameters->at(1);
    const ldf re=parameters->at(3);
    return d*pow(1.0-exp(a*(re-r)),2);
}

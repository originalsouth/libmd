template<class X> ldf COULOMB(X r,ldf rsq,vector<ldf> *parameters)
{
    (void) rsq;
    const ldf q=parameters->at(0);
    return q/r;
}

template<class X> ldf YUKAWA(X r,ldf rsq,vector<ldf> *parameters)
{
    (void) rsq;
    const ldf b=parameters->at(0);
    const ldf k=parameters->at(1);
    return b/(r*exp(k*r));
}

template<class X> ldf HOOKIAN(X r,ldf rsq,vector<ldf> *parameters)
{
    (void) rsq;
    const ldf k=parameters->at(0);
    const ldf r0=parameters->at(1);
    return k/2.0*pow(r-r0,2);
}

template<class X> ldf LJ(X r,ldf rsq,vector<ldf> *parameters)
{
    (void) rsq;
    const ldf e=parameters->at(0);
    const ldf s=parameters->at(1);
    return 4.0*e*(pow(s/r,12)-pow(s/r,6));
}

template<class X> ldf MORSE(X r,ldf rsq,vector<ldf> *parameters)
{
    (void) rsq;
    const ldf d=parameters->at(0);
    const ldf a=parameters->at(1);
    const ldf re=parameters->at(3);
    return d*pow(1.0-exp(a*(re-r)),2);
}

ldf dCOULOMBdr(ldf r,ldf rsq,vector<ldf> *parameters)
{
    (void) r;
    const ldf q=parameters->at(0);
    return -q/rsq;
}

ldf dYUKAWAdr(ldf r,ldf rsq,vector<ldf> *parameters)
{
    const ldf b=parameters->at(0);
    const ldf k=parameters->at(1);
    return b*(k*r+1.0)/(rsq*exp(k*r));
}
ldf dHOOKIANdr(ldf r,ldf rsq,vector<ldf> *parameters)
{
    (void) rsq;
    const ldf k=parameters->at(0);
    const ldf r0=parameters->at(1);
    return k*(r-r0);
}

ldf dLJdr(ldf r,ldf rsq,vector<ldf> *parameters)
{
    (void) r;
    const ldf e=parameters->at(0);
    const ldf s=parameters->at(1);
    return (24.0*e/r)*(pow(s,12)/pow(rsq,6)-pow(s,6)/pow(rsq,3));
}

ldf dMORSEdr(ldf r,ldf rsq,vector<ldf> *parameters)
{
    (void) rsq;
    const ldf d=parameters->at(0);
    const ldf a=parameters->at(1);
    const ldf re=parameters->at(3);
    return -2.0*a*d*exp(a*(re-r))*(exp(a*(re-r))-1.0);
}

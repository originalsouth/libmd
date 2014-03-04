#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> mp<dim>::mp()
{
    setmp();
}

template<ui dim> void mp<dim>::setmp(ui i)
{
    switch(i)
    {
        case MP::MP_GAUSSIANBUMP:
            parameters.assign(2,1);
            fmp=&GAUSSIANBUMP<dim>;
            dfmp=&dGAUSSIANBUMP<dim>;
            ddfmp=&ddGAUSSIANBUMP<dim>;
        break;
        default:
            parameters.assign(1,1);
            fmp=&FLATSPACE<dim>;
            dfmp=&dFLATSPACE<dim>;
            ddfmp=&ddFLATSPACE<dim>;
        break;
    }
}

template<ui dim> void mp<dim>::setmp(fmpptr f,dfmpptr df,ddfmpptr ddf)
{
    fmp=f;
    dfmp=df;
    ddfmp=ddf;
}

template<ui dim> ldf mp<dim>::f(ldf x[dim])
{
    return fmp(x,&parameters);
}

template<ui dim> ldf mp<dim>::df(ui i,ldf x[dim])
{
    return dfmp(i,x,&parameters);
}

template<ui dim> ldf mp<dim>::ddf(ui i,ui j,ldf x[dim])
{
    return ddfmp(i,j,x,&parameters);
}

template<ui dim> ldf mp<dim>::g(ui i,ui j,ldf x[dim])
{
    return kdelta(i,j)+df(i,x)*df(j,x);
}

template<ui dim> ldf mp<dim>::ginv(ui i,ui j,ldf x[dim])
{
    ldf det=1.0;
    for(ui d=0;d<dim;d++) det+=pow(df(d,x),2);
    return kdelta(i,j)-(df(i,x)*df(j,x))/det;
}

template<ui dim> ldf mp<dim>::dg(ui s,ui i,ui j,ldf x[dim])
{
    return ddf(s,i,x)*df(j,x)+df(i,x)*ddf(s,j,x);
}

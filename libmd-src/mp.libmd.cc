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
        case MP::MP_MOLLIFIER:
            parameters.assign(2,1.0);
            fmp=&MOLLIFIER<ldf,dim>;
            dfmp=&MOLLIFIER<duals<dim>,dim>;
        break;
        case MP::MP_EGGCARTON:
            parameters.assign(dim+1,1.0);
            fmp=&EGGCARTON<ldf,dim>;
            dfmp=&EGGCARTON<duals<dim>,dim>;
        break;
        case MP::MP_GAUSSIANBUMP:
            parameters.assign(2,1.0);
            fmp=&GAUSSIANBUMP<ldf,dim>;
            dfmp=&GAUSSIANBUMP<duals<dim>,dim>;
        break;
        default:
            parameters.assign(1,1.0);
            fmp=&FLATSPACE<ldf,dim>;
            dfmp=&FLATSPACE<duals<dim>,dim>;
        break;
    }
}

template<ui dim> void mp<dim>::setmp(fmpptr<ldf,dim> f,fmpptr<duals<dim>,dim> df)
{
    fmp=f;
    dfmp=df;
}

template<ui dim> void mp<dim>::calc(ui i,ldf x[dim])
{
    duals<dim> y[dim];
    for(ui d=0;d<dim;d++) y[d]=duals<dim>(x[d],d);
    geometryx[i]=dfmp(y,&parameters);
}

template<ui dim> ldf mp<dim>::f(ldf x[dim])
{
    return fmp(x,&parameters);
}

template<ui dim> ldf mp<dim>::df(ui mu,ldf x[dim])
{
    duals<dim> y[dim];
    for(ui d=0;d<dim;d++) y[d]=duals<dim>(x[d],d);
    return fmp(y,&parameters).dx[mu];
}

template<ui dim> ldf mp<dim>::df(ui mu,ui nu,ldf x[dim])
{
    duals<dim> y[dim];
    for(ui d=0;d<dim;d++) y[d]=duals<dim>(x[d],d);
    return fmp(y,&parameters).dxdy[mu][nu];
}

template<ui dim> ldf mp<dim>::g(ui i,ui mu,ui nu)
{
    return kdelta(mu,nu)+geometryx[i].dx[mu]*geometryx[i].dx[nu];
}

template<ui dim> ldf mp<dim>::gp(ui i,ui mu,ui nu)
{
    return kdelta(mu,nu)+geometryxp[i].dx[mu]*geometryxp[i].dx[nu];
}

template<ui dim> ldf mp<dim>::ginv(ui i,ui mu,ui nu)
{
    ldf det=1.0;
    for(ui d=0;d<dim;d++) det+=pow(geometryx[i].dx[d],2);
    return kdelta(mu,nu)-(geometryx[i].dx[mu]*geometryx[i].dx[nu])/det;
}

template<ui dim> ldf mp<dim>::G(ui i,ui sigma,ui mu,ui nu)
{
    return geometryx[i].dx[nu]*geometryx[i].dxdy[sigma][mu];
}

#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> mp<dim>::mp()
{
    //!
    //! Constructor for the mp<dim> class
    //! By default loads a FLATSPACE Monge patch by calling setmp
    //!
    setmp();
}

template<ui dim> void mp<dim>::setmp(ui i)
{
    //!
    //! Sets a builtin Monge patch
    //! <tt>i</tt> can be a <tt>ui</tt> or <tt>MP enum</tt> type e.g. <tt>MP::FLATSPACE</tt> (default e.g. <tt>MP::FLATSPACE</tt> (default).
    //!
    patch=i;
    switch(i)
    {
        case MP::MOLLIFIER:
            parameters.assign(2,1.0);
            fmp=&MOLLIFIER<ldf,dim>;
            dfmp=&MOLLIFIER<duals<dim>,dim>;
        break;
        case MP::EGGCARTON:
            parameters.assign(dim+1,1.0);
            fmp=&EGGCARTON<ldf,dim>;
            dfmp=&EGGCARTON<duals<dim>,dim>;
        break;
        case MP::GAUSSIANBUMP:
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
    //!
    //! Sets an externaly defined Monge patch
    //!
    fmp=f;
    dfmp=df;
    patch=UI_MAX;
}

template<ui dim> void mp<dim>::calc(ui i,ldf x[dim])
{
    //!
    //! Calculates geometric information for particle <tt>i</tt> at position <tt>x</tt>
    //! The results are stored in the geometryx array
    //!
    duals<dim> y[dim];
    for(ui d=0;d<dim;d++) y[d]=duals<dim>(x[d],d);
    geometryx[i]=dfmp(y,parameters);
}

template<ui dim> void mp<dim>::calc(duals<dim> &z,ldf x[dim])
{
    //!
    //! Calculates geometric information for particle <tt>i</tt> at position <tt>x</tt>
    //! The results are stored in the first argument <tt>z</tt>
    //!
    duals<dim> y[dim];
    for(ui d=0;d<dim;d++) y[d]=duals<dim>(x[d],d);
    z=dfmp(y,parameters);
}

template<ui dim> ldf mp<dim>::f(ldf x[dim])
{
    //!
    //! Calculates value of the (set) Monge function at position <tt>x</tt>
    //!
    return fmp(x,parameters);
}

template<ui dim> ldf mp<dim>::df(ui mu,ldf x[dim])
{
    //!
    //! Calculates the derivative in direction <tt>mu</tt> of the(set) Monge function at position <tt>x</tt>
    //!
    duals<dim> y[dim];
    for(ui d=0;d<dim;d++) y[d]=duals<dim>(x[d],d);
    return fmp(y,parameters).dx[mu];
}

template<ui dim> ldf mp<dim>::ddf(ui mu,ui nu,ldf x[dim])
{
    //!
    //! Calculates the second derivative in the directions <tt>mu nu</tt> of the(set) Monge function at position <tt>x</tt>
    //!
    duals<dim> y[dim];
    for(ui d=0;d<dim;d++) y[d]=duals<dim>(x[d],d);
    return fmp(y,parameters).dxdy[mu][nu];
}

template<ui dim> ldf mp<dim>::g(ui i,ui mu,ui nu)
{
    //!
    //! Calculates Monge metric tensor element <tt>mu</tt><tt>nu</tt> for particle <tt>i</tt>
    //!
    return kdelta(mu,nu)+geometryx[i].dx[mu]*geometryx[i].dx[nu];
}

template<ui dim> ldf mp<dim>::gp(ui i,ui mu,ui nu)
{
    //!
    //! Calculates Monge metric tensor element <tt>mu</tt><tt>nu</tt> for particle <tt>i</tt> at it's previous position <tt>xp</tt>
    //!
    return kdelta(mu,nu)+geometryxp[i].dx[mu]*geometryxp[i].dx[nu];
}

template<ui dim> ldf mp<dim>::ginv(ui i,ui mu,ui nu)
{
    //!
    //! Calculates Monge metric tensor inverse element <tt>mu</tt><tt>nu</tt> for particle <tt>i</tt>
    //!
    ldf det=1.0;
    for(ui d=0;d<dim;d++) det+=std::pow(geometryx[i].dx[d],2);
    return kdelta(mu,nu)-(geometryx[i].dx[mu]*geometryx[i].dx[nu])/det;
}

template<ui dim> ldf mp<dim>::sqrt_ginv(ui i,ui mu,ui nu)
{
    //!
    //! Calculates the square root of Monge metric tensor inverse element <tt>mu</tt><tt>nu</tt> for particle <tt>i</tt>
    //!
    ldf det=1.0;
    for(ui d=0;d<dim;d++) det+=std::pow(geometryx[i].dx[d],2);
    return kdelta(mu,nu)-(geometryx[i].dx[mu]*geometryx[i].dx[nu])/(det+std::sqrt(det));
}

template<ui dim> ldf mp<dim>::A(ui i,ui sigma,ui mu,ui nu)
{
    //!
    //! Calculates A symbol of the van Zuiden integrator element <tt>sigma mu nu</tt> for particle <tt>i</tt>
    //!
    return geometryx[i].dx[nu]*geometryx[i].dxdy[sigma][mu];
}

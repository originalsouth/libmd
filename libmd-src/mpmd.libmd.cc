#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> mpmd<dim>::mpmd():md<dim>()
{

}

template<ui dim> mpmd<dim>::mpmd(ui particlenr):md<dim>(particlenr)
{

}

template<ui dim> ldf mpmd<dim>::embedded_distsq(ui p1,ui p2)
{
    return distsq(p1,p2)+pow(patch.f(particles[p2].x)-patch.f(particles[p1].x),2);
}

template<ui dim> ldf mpmd<dim>::embedded_dd_p1(ui i,ui p1,ui p2)
{
    return dd(i,p1,p2)+((patch.f(particles[p2].x)-patch.f(particles[p1].x))*patch.df(i,particles[p1].x));
}

template<ui dim> ldf mpmd<dim>::embedded_dd_p2(ui i,ui p1,ui p2)
{
    return dd(i,p2,p1)+((patch.f(particles[p1].x)-patch.f(particles[p2].x))*patch.df(i,particles[p2].x));
}

template<ui dim> void mpmd<dim>::zuiden_C(ui i,ldf C[dim])
{
    ldf temp[dim];
    memset(temp,0,dim*sizeof(ldf));
    memset(C,0,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) for(ui mu=0;mu<dim;mu++) temp[d]+=patch.g(d,mu,particles[i].xp)*(particles[i].x[mu]-particles[i].xp[mu]);
    for(ui d=0;d<dim;d++) temp[d]+=pow(integrator.h,2)/particles[i].m*particles[i].F[d];
    for(ui d=0;d<dim;d++) for(ui sigma=0;sigma<dim;sigma++) C[d]+=patch.ginv(d,sigma,particles[i].x)*temp[sigma];
}

template<ui dim> void mpmd<dim>::zuiden_A(ui i,ldf eps[dim])
{
    ldf temp[dim];
    memcpy(temp,eps,dim*sizeof(ldf));
    memset(eps,0,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) for(ui mu=0;mu<dim;mu++) for(ui nu=0;nu<dim;nu++) for(ui sigma=0;sigma<dim;sigma++) eps[d]+=patch.ginv(d,sigma,particles[i].x)*patch.dg(sigma,mu,nu,particles[i].x)*temp[mu]*temp[nu];
}

template<ui dim> void mpmd<dim>::thread_zuiden_wfi(ui i)
{
    ldf eps[dim];
    zuiden_C(i,eps);
    for(ui d=0;d<dim;d++) particles[i].xp[d]=particles[i].x[d];
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden_protect(ui i)
{
    ui count=0;
    ldf val,eps[dim],epsp[dim],C[dim];
    zuiden_C(i,C);
    memcpy(eps,C,dim*sizeof(ldf));
    memcpy(epsp,C,dim*sizeof(ldf));
    do
    {
        val=0.0;
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=C[d];
        for(ui d=0;d<dim;d++) val+=fabs(eps[d]-epsp[d]);
        memcpy(epsp,eps,dim*sizeof(ldf));
        count++;
    }
    while(count<integrator.generations or val>numeric_limits<ldf>::epsilon());
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden(ui i)
{
    ldf val,eps[dim],epsp[dim],C[dim];
    zuiden_C(i,C);
    memcpy(eps,C,dim*sizeof(ldf));
    memcpy(epsp,C,dim*sizeof(ldf));
    do
    {
        val=0.0;
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=C[d];
        for(ui d=0;d<dim;d++) val+=fabs(eps[d]-epsp[d]);
        memcpy(epsp,eps,dim*sizeof(ldf));
    }
    while(val>numeric_limits<ldf>::epsilon());
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_history(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].xp[d]=particles[i].x[d]-particles[i].dx[d]*integrator.h;
}

template<ui dim> void mpmd<dim>::history()
{
    for(ui i=0;i<N;i++) thread_history(i);
}

template<ui dim> void mpmd<dim>::thread_calc_forces(ui i)
{
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>network.skins[i][j].neighbor)
    {
        const ldf rsq=embedded_distsq(i,network.skins[i][j].neighbor);
        if(!network.update or (network.update and rsq<network.rcosq))
        {
            const ldf r=sqrt(rsq);
            const ldf dVdr=v.dr(network.library[network.skins[i][j].interaction].potential,r,&network.library[network.skins[i][j].interaction].parameters);
            for(ui d=0;d<dim;d++)
            {
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]+=embedded_dd_p1(d,i,network.skins[i][j].neighbor)*dVdr/r;
                particles[network.skins[i][j].neighbor].F[d]+=embedded_dd_p2(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]+=embedded_dd_p1(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #pragma omp atomic
                particles[network.skins[i][j].neighbor].F[d]+=embedded_dd_p2(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #else
                particles[i].F[d]+=embedded_dd_p1(d,i,network.skins[i][j].neighbor)*dVdr/r;
                particles[network.skins[i][j].neighbor].F[d]+=embedded_dd_p2(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #endif
            }
        }
    }
}

template<ui dim> void mpmd<dim>::integrate()
{
    switch(integrator.method)
    {
        case MP_INTEGRATOR::MP_VVERLET:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_vverlet_x(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_calc_forces(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_vverlet_dx(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_x(i);
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_calc_forces(i);
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_x(i);
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_calc_forces(i);
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
            #endif
        break;
        case MP_INTEGRATOR::MP_SEULER:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_seuler(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_seuler(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_seuler(i);
            #endif
        break;
        case MP_INTEGRATOR::MP_VZ_WFI:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_zuiden_wfi(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_wfi(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_wfi(i);
            #endif
        break;
        case MP_INTEGRATOR::MP_VZ_P:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_zuiden_protect(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_protect(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_protect(i);
            #endif
        break;
        default:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_zuiden(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden(i);
            #endif
        break;
    }
    periodicity();
}

template<ui dim> ldf mpmd<dim>::thread_T(ui i)
{
    ldf retval=0.0;
    for(ui mu=0;mu<dim;mu++) for(ui nu=0;nu<dim;nu++) retval+=patch.g(mu,nu,particles[i].x)*particles[i].dx[mu]*particles[i].dx[nu];
    return 0.5*particles[i].m*retval;
}

template<ui dim> ldf mpmd<dim>::thread_V(ui i)
{
    ldf retval=0.0;
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i<network.skins[i][j].neighbor)
    {
        const ldf rsq=embedded_distsq(i,network.skins[i][j].neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            retval+=(v(network.library[network.skins[i][j].interaction].potential,r,&network.library[network.skins[i][j].interaction].parameters)-network.library[network.skins[i][j].interaction].vco);
        }
    }
    return retval;
}

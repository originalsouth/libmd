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
    return dd(i,p1,p2)+((patch.geometryx[p2].x-patch.geometryx[p1].x)*patch.geometryx[p1].dx[i]);
}

template<ui dim> ldf mpmd<dim>::embedded_dd_p2(ui i,ui p1,ui p2)
{
    return dd(i,p2,p1)+((patch.geometryx[p1].x-patch.geometryx[p2].x)*patch.geometryx[p2].dx[i]);
}

template<ui dim> void mpmd<dim>::zuiden_C(ui i,ldf ZC[dim])
{
    ldf C[dim]={};
    memset(ZC,0,dim*sizeof(ldf));
    for(ui sigma=0;sigma<dim;sigma++) for(ui mu=0;mu<dim;mu++) C[sigma]+=patch.gp(i,sigma,mu)*(particles[i].x[mu]-particles[i].xp[mu]);
    for(ui sigma=0;sigma<dim;sigma++) C[sigma]+=pow(integrator.h,2)*particles[i].F[sigma]/particles[i].m;
    for(ui rho=0;rho<dim;rho++) for(ui sigma=0;sigma<dim;sigma++) ZC[rho]+=patch.ginv(i,rho,sigma)*C[sigma];
}

template<ui dim> void mpmd<dim>::zuiden_A(ui i,ldf eps[dim])
{
    ldf ZA[dim]={};
    for(ui rho=0;rho<dim;rho++) for(ui sigma=0;sigma<dim;sigma++) for(ui mu=0;mu<dim;mu++) for(ui nu=0;nu<dim;nu++) ZA[rho]+=patch.ginv(i,rho,sigma)*patch.G(i,sigma,mu,nu)*eps[mu]*eps[nu];
    memcpy(eps,ZA,dim*sizeof(ldf));
}

template<ui dim> void mpmd<dim>::thread_zuiden_wfi(ui i)
{
    ldf eps[dim]={};
    zuiden_C(i,eps);
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden_protect(ui i)
{
    ui count=0;
    ldf ZC[dim],eps[dim],epsp[dim],val;
    zuiden_C(i,ZC);
    memcpy(eps,ZC,dim*sizeof(ldf));
    do
    {
        val=0.0;
        count++;
        memcpy(epsp,eps,dim*sizeof(ldf));
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=ZC[d];
        for(ui d=0;d<dim;d++) val+=fabs(epsp[d]-eps[d]);
        DEBUG_3("fixed point iterators cycle: %u",count);
    }
    while(count<integrator.generations and val/dim>numeric_limits<ldf>::epsilon());
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden(ui i)
{
    ldf ZC[dim],eps[dim],epsp[dim],val;
    zuiden_C(i,ZC);
    memcpy(eps,ZC,dim*sizeof(ldf));
    do
    {
        val=0.0;
        memcpy(epsp,eps,dim*sizeof(ldf));
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=ZC[d];
        for(ui d=0;d<dim;d++) val+=fabs(epsp[d]-eps[d]);
    }
    while(val/dim>numeric_limits<ldf>::epsilon());
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
    patch.geometryx.resize(N);
    patch.geometryxp.resize(N);
    for(ui i=0;i<N;i++) patch.calc(i,particles[i].xp);
    for(ui i=0;i<N;i++) patch.geometryxp[i]=patch.geometryx[i];
    for(ui i=0;i<N;i++) patch.calc(i,particles[i].x);
}

template<ui dim> void mpmd<dim>::thread_calc_geometry(ui i)
{
    patch.calc(i,particles[i].x);
}

template<ui dim> void mpmd<dim>::calc_geometry()
{
    if(N!=patch.geometryx.size())
    {
        patch.geometryx.resize(N);
        patch.geometryxp.resize(N);
    }
    for(ui i=0;i<N;i++) thread_calc_geometry(i);
}

template<ui dim> void mpmd<dim>::mp_thread_calc_forces(ui i)
{
    for(auto sij: network.skins[i]) if(i>sij.neighbor)
    {
        const ldf rsq=embedded_distsq(i,sij.neighbor);
        if(!network.update or rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            DEBUG_3("r = %Lf",r);
            const ldf dVdr=v.dr(network.library[sij.interaction].potential,r,&network.library[sij.interaction].parameters);
            DEBUG_3("dV/dr = %Lf",dVdr);
            for(ui d=0;d<dim;d++)
            {
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]+=embedded_dd_p1(d,i,sij.neighbor)*dVdr/r;
                particles[sij.neighbor].F[d]+=embedded_dd_p2(d,i,sij.neighbor)*dVdr/r;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]+=embedded_dd_p1(d,i,sij.neighbor)*dVdr/r;
                #pragma omp atomic
                particles[sij.neighbor].F[d]+=embedded_dd_p2(d,i,sij.neighbor)*dVdr/r;
                #else
                particles[i].F[d]+=embedded_dd_p1(d,i,sij.neighbor)*dVdr/r;
                DEBUG_3("particles[%u].F[d] = %Lf",i,embedded_dd_p1(d,i,sij.neighbor)*dVdr/r);
                particles[sij.neighbor].F[d]+=embedded_dd_p2(d,i,sij.neighbor)*dVdr/r;
                DEBUG_3("particles[%u].F[d] = %Lf",sij.neighbor,embedded_dd_p2(d,i,sij.neighbor)*dVdr/r);
                #endif
            }
        }
    }
    if(network.forcelibrary.size() and network.forces[i].size()) for(ui j=network.forces[i].size()-1;j<numeric_limits<ui>::max();j--)
    {
        ui ftype=network.forces[i][j];
        if(network.forcelibrary[ftype].particles.size() and network.forcelibrary[ftype].particles[i].size()) f(network.forcelibrary[ftype].externalforce,i,&network.forcelibrary[ftype].particles[i],&network.forcelibrary[ftype].parameters,(md<dim>*)this);
        else f(network.forcelibrary[ftype].externalforce,i,nullptr,&network.forcelibrary[ftype].parameters,(md<dim>*)this);
    }
}

template<ui dim> void mpmd<dim>::calc_forces()
{
    DEBUG_2("exec is here");
    avars.export_force_calc=false;
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_clear_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) mp_thread_calc_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    #pragma omp parallel for
    for(ui i=0;i<N;i++) mp_thread_calc_forces(i);
    #else
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    for(ui i=0;i<N;i++) mp_thread_calc_forces(i);
    #endif
}

template<ui dim> void mpmd<dim>::recalc_forces()
{
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) mp_thread_calc_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) mp_thread_calc_forces(i);
    #else
    for(ui i=0;i<N;i++) mp_thread_calc_forces(i);
    #endif

}
template<ui dim> void mpmd<dim>::integrate()
{
    switch(integrator.method)
    {
        case MP_INTEGRATOR::VVERLET:
            WARNING("flatspace integrator");
            DEBUG_2("integrating using flatspace velocity Verlet");
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
            for(ui i=0;i<N;i++) if(!particles[i].fix) mp_thread_calc_forces(i);
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
            #endif
        break;
        case MP_INTEGRATOR::SEULER:
            WARNING("flatspace integrator");
            DEBUG_2("integrating using flatspace symplectic Euler");
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
        case MP_INTEGRATOR::VZ_WFI:
            WARNING("incomplete integrator");
            DEBUG_2("integrating using van Zuiden without fixed point iterations");
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
        case MP_INTEGRATOR::VZ_P:
            DEBUG_2("integrating using van Zuiden with protected (with %u generations) fixed point iterations",integrator.generations);
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
            DEBUG_2("integrating using van Zuiden");
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
    for(ui i=0;i<N;i++) patch.geometryxp[i]=patch.geometryx[i];
    calc_geometry();
}

template<ui dim> ldf mpmd<dim>::thread_T(ui i)
{
    ldf retval=0.0;
    for(ui mu=0;mu<dim;mu++) for(ui nu=0;nu<dim;nu++) retval+=patch.g(i,mu,nu)*particles[i].dx[mu]*particles[i].dx[nu];
    return 0.5*particles[i].m*retval;
}

template<ui dim> ldf mpmd<dim>::thread_V(ui i)
{
    ldf retval=0.0;
    for(auto sij: network.skins[i]) if(i<sij.neighbor)
    {
        const ldf rsq=embedded_distsq(i,sij.neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            retval+=(v(network.library[sij.interaction].potential,r,&network.library[sij.interaction].parameters)-network.library[sij.interaction].vco);
        }
    }
    return retval;
}

#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_clear_forces(ui i)
{
    memset(particles[i].F,0,dim*sizeof(ldf));
}

template<ui dim> void md<dim>::thread_calc_forces(ui i)
{
    ldf rcosq=pow(network.rco,2);
    for(auto sij: network.skins[i]) if(i>sij.neighbor)
    {
        ldf rsq=distsq(i,sij.neighbor);
        if(!network.update or rsq<rcosq)
        {
            ldf r=sqrt(rsq);
            DEBUG_3("r = %Lf",r);
            ldf dVdr=v.dr(network.library[sij.interaction].potential,r,&network.library[sij.interaction].parameters);
            DEBUG_3("dV/dr = %Lf",dVdr);
            for(ui d=0;d<dim;d++)
            {
                ldf F_i=dd(d,i,sij.neighbor)*dVdr/r;
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]+=F_i;
                particles[sij.neighbor].F[d]-=F_i;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]+=F_i;
                #pragma omp atomic
                particles[sij.neighbor].F[d]-=F_i;
                #else
                particles[i].F[d]+=F_i;
                DEBUG_3("particles[%u].F[d] = %Lf",i,F_i);
                particles[sij.neighbor].F[d]-=F_i;
                DEBUG_3("particles[%u].F[d] = %Lf",sij.neighbor,-F_i);
                #endif
            }
        }
    }
    for(auto ftype: network.forces[i])
        f(network.forcelibrary[ftype].externalforce,i,&network.forcelibrary[ftype].particles[i],(!network.forcelibrary[ftype].particles.empty() and !network.forcelibrary[ftype].particles[i].empty())?&network.forcelibrary[ftype].parameters:nullptr,(md<dim>*)this);
}

template<ui dim> void md<dim>::calc_forces()
{
    if(network.update and (avars.reindex or test_index()))
    {
        DEBUG_2("regenerating skinlist");
        index();
    }
    DEBUG_2("exec is here");
    avars.export_force_calc=false;
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_clear_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_calc_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #else
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #endif
}

template<ui dim> void md<dim>::recalc_forces()
{
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_calc_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #else
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #endif
}

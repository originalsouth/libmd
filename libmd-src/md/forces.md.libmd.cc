#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_clear_forces(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].F[d]=0.0;
}

template<ui dim> void md<dim>::thread_calc_forces(ui i)
{
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>network.skins[i][j].neighbor)
    {
        const ldf rsq=distsq(i,network.skins[i][j].neighbor);
        if(!network.update or (network.update and rsq<network.rcosq))
        {
            const ldf r=sqrt(rsq);
            DEBUG_3("r = %Lf",r);
            const ldf dVdr=v.dr(network.library[network.skins[i][j].interaction].potential,r,&network.library[network.skins[i][j].interaction].parameters);
            DEBUG_3("dV/dr = %Lf",dVdr);
            for(ui d=0;d<dim;d++)
            {
                const ldf F_i=dd(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]+=F_i;
                particles[network.skins[i][j].neighbor].F[d]-=F_i;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]+=F_i;
                #pragma omp atomic
                particles[network.skins[i][j].neighbor].F[d]-=F_i;
                #else
                particles[i].F[d]+=F_i;
                DEBUG_3("particles[%u].F[d] = %Lf",i,F_i);
                particles[network.skins[i][j].neighbor].F[d]-=F_i;
                DEBUG_3("particles[%u].F[d] = %Lf",network.skins[i][j].neighbor,-F_i);
                #endif
            }
        }
    }
    if(network.forcelibrary.size() and network.forces[i].size()) for(ui j=network.forces[i].size()-1;j<numeric_limits<ui>::max();j--)
    {
        ui ftype=network.forces[i][j];
        if(network.forcelibrary[ftype].particles.size() and network.forcelibrary[ftype].particles[i].size())
        {
            vector<particle<dim>*> plist;
            uitopptr(&plist,network.forcelibrary[ftype].particles[i]);
            f(network.forcelibrary[ftype].externalforce,&particles[i],&plist,&network.forcelibrary[ftype].parameters,this);
        }
        else f(network.forcelibrary[ftype].externalforce,&particles[i],nullptr,&network.forcelibrary[ftype].parameters,this);
    }
}

template<ui dim> void md<dim>::calc_forces()
{
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

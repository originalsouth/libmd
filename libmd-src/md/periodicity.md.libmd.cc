#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_periodicity_boxshear(ui d,ui i) //FIXME: Jayson
{
    ldf boundaryCrossing=round(particles[i].x[d]/simbox.L[d]);
    if(fabs(boundaryCrossing)>0.1) for(ui k=0;k<dim;k++)
    {
        particles[i].x[k]-=simbox.Lshear[k][d]*boundaryCrossing;
        particles[i].xp[k]-=simbox.Lshear[k][d]*boundaryCrossing;
        particles[i].dx[k]-=simbox.vshear[k][d]*boundaryCrossing;
    }
}

template<ui dim> void md<dim>::thread_periodicity_hard(ui d,ui i)
{
    ldf xnew=simbox.L[d]*(fabs(particles[i].x[d]/simbox.L[d]+0.5-2.0*floor(particles[i].x[d]/(2.0*simbox.L[d])+0.75))-0.5);
    ldf sign=(((int)round(particles[i].x[d]/simbox.L[d]))&1?-1.0:1.0);
    particles[i].xp[d]+=sign*(xnew-particles[i].x[d]);
    particles[i].x[d]=xnew;
    particles[i].dx[d]*=sign;
}

template<ui dim> void md<dim>::periodicity()
{
    if(simbox.bcond) for(ui d=0;d<dim;d++) switch(simbox.bcond[d])
    {
        case BCOND::PERIODIC:
        {
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_periodicity_periodic(d,i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_periodic(d,i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_periodic(d,i);
            #endif
        }
        break;
        case BCOND::BOXSHEAR:
        {
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_periodicity_boxshear(d,i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_boxshear(d,i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_boxshear(d,i);
            #endif
        }
        break;
        case BCOND::HARD:
        {
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_periodicity_hard(d,i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_hard(d,i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_hard(d,i);
            #endif
        }
        break;
    }
}

#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_seuler(ui i)
{
    const ldf o=integrator.h/particles[i].m;
    for(ui d=0;d<dim;d++)
    {
        particles[i].dx[d]+=o*particles[i].F[d];
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=integrator.h*particles[i].dx[d];
    }
}

template<ui dim> void md<dim>::thread_vverlet_x(ui i)
{
    for(ui d=0;d<dim;d++)
    {
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=integrator.h*particles[i].dx[d]+0.5*pow(integrator.h,2)*particles[i].F[d]/particles[i].m;
    }
}

template<ui dim> void md<dim>::thread_vverlet_dx(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].dx[d]+=0.5*integrator.h*particles[i].F[d]/particles[i].m;
}

template<ui dim> void md<dim>::integrate()
{
    avars.export_force_calc=true;
    switch(integrator.method)
    {
        case INTEGRATOR::VVERLET:
            DEBUG_2("integrating using symplectic Velocity Verlet");
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
            recalc_forces();
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
            #endif
        break;
        default:
            DEBUG_2("integrating using symplectic Euler method");
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
    }
    periodicity();
}

template<ui dim> void md<dim>::update_boundaries()
{
    // update box matrix for shear
    for(ui j=0;j<dim;j++) for (ui k=0; k<dim; k++)
    {
        simbox.Lshear[j][k] += simbox.vshear[j][k]*integrator.h;
        // shift by appropriate box lengths so that the off-diagonal entries are in the range -L[j][j]/2 to L[j][j]/2 consistent with the positions
        if (j != k)
        {
            while(simbox.Lshear[j][k]>simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]-=simbox.Lshear[j][j];
            while(simbox.Lshear[j][k]<-simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]+=simbox.Lshear[j][j];
        }
    }
    simbox.invert_box();
}

template<ui dim> void md<dim>::timestep()
{
    if(network.update and test_index())
    {
        DEBUG_2("regenerating skinlist");
        index();
    }
    if(simbox.boxShear)
    {
        DEBUG_2("updating boxshear boundaries");
        update_boundaries();
    }
    calc_forces();
    integrate();
}

template<ui dim> void md<dim>::timesteps(ui k)
{
    for(ui i=0;i<k;i++) timestep();
}

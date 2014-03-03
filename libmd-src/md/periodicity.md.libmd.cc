#ifndef libmd_h
#include "../../libmd.h"
#endif


template<ui dim> void md<dim>::thread_periodicity_periodic(ui d,ui i)
{
    ldf dx=simbox.L[d]*round(particles[i].x[d]/simbox.L[d]);
    particles[i].xp[d]-=dx;
    particles[i].x[d]-=dx;
}


template<ui dim> void md<dim>::thread_periodicity_boxshear(ui d,ui i) 
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
    if (simbox.boxShear) {
        ldf s=0,sdot=0;
        for (ui k=0;k<dim;k++) { 
            s+=simbox.LshearInv[d][k]*particles[i].x[k]; 
            sdot+=simbox.LshearInv[d][k]*particles[i].dx[k]; // note: assumes that vshear does not directly affect dynamics (only indirectly via PBC in shear direction)
        }
        if (fabs(s) > 0.5) { // particle has hit the hard boundary as distorted by the shear
            ldf snew = fabs(s+0.5-2.0*floor(s/2.0+0.75))-0.5;
            ldf sign=(((int)round(s))&1?-1.0:1.0);
            // but now, need to calculate the corrections in all directions due to reflection in s-space
            for (ui j=0;j<dim;j++) {
                particles[i].xp[j]+=sign*simbox.Lshear[j][d]*(snew-s);
                particles[i].x[j]+=simbox.Lshear[j][d]*(snew-s);
                particles[i].dx[j]+=simbox.Lshear[j][d]*sdot*(sign-1);
            }
        }
    }
    else {
        ldf xnew=simbox.L[d]*(fabs(particles[i].x[d]/simbox.L[d]+0.5-2.0*floor(particles[i].x[d]/(2.0*simbox.L[d])+0.75))-0.5);
        ldf sign=(((int)round(particles[i].x[d]/simbox.L[d]))&1?-1.0:1.0);
        particles[i].xp[d]+=sign*(xnew-particles[i].x[d]);
        particles[i].x[d]=xnew;
        particles[i].dx[d]*=sign;
    }
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


template<ui dim> void md<dim>::update_boundaries()
{
    // update box matrix for shear
    for(ui j=0;j<dim;j++) for (ui k=0; k<dim; k++)
    {
        simbox.Lshear[j][k] += simbox.vshear[j][k]*integrator.h;
        // shift by appropriate box lengths so that the off-diagonal entries are in the range -L[j][j]/2 to L[j][j]/2 consistent with the positions
        if (j != k && (simbox.bcond[j]==BCOND::PERIODIC || simbox.bcond[j]==BCOND::BOXSHEAR))
        {
            while(simbox.Lshear[j][k]>simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]-=simbox.Lshear[j][j];
            while(simbox.Lshear[j][k]<-simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]+=simbox.Lshear[j][j];
        }
    }
    simbox.invert_box();
}


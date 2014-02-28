#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::import_pos(ldf *x)
{
    ui d=vvars[0];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].x[d]=x[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
    #else
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
    #endif
}

template<ui dim> template<typename...arg> void md<dim>::import_pos(ldf *x,arg...argv)
{
    ui d=vvars[0];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].x[d]=x[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
    #else
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
    #endif
    import_pos(argv...);
}

template<ui dim> void md<dim>::import_pos(ui i,ldf x)
{
    ui d=vvars[0];
    particles[i].x[d]=x;
}

template<ui dim> template<typename...arg> void md<dim>::import_pos(ui i,ldf x,arg...argv)
{
    ui d=vvars[0];
    particles[i].x[d]=x;
    import_pos(i,argv...);
}

template<ui dim> void md<dim>::import_vel(ldf *dx)
{
    ui d=vvars[1];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].dx[d]=dx[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
    #else
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
    #endif
}

template<ui dim> template<typename...arg> void md<dim>::import_vel(ldf *dx,arg...argv)
{
    ui d=vvars[1];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].dx[d]=dx[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
    #else
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
    #endif
    import_vel(argv...);
}

template<ui dim> void md<dim>::import_vel(ui i,ldf dx)
{
    ui d=vvars[1];
    particles[i].dx[d]=dx;
}

template<ui dim> template<typename...arg> void md<dim>::import_vel(ui i,ldf dx,arg...argv)
{
    ui d=vvars[1];
    particles[i].dx[d]=dx;
    import_vel(i,argv...);
}

template<ui dim> void md<dim>::import_force(ldf *F)
{
    ui d=vvars[2];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].F[d]=F[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
    #else
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
    #endif
}

template<ui dim> template<typename...arg> void md<dim>::import_force(ldf *F,arg...argv)
{
    ui d=vvars[2];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].F[d]=F[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
    #else
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
    #endif
    import_force(argv...);
}

template<ui dim> void md<dim>::import_force(ui i,ldf F)
{
    ui d=vvars[2];
    particles[i].F[d]=F;
}

template<ui dim> template<typename...arg> void md<dim>::import_force(ui i,ldf F,arg...argv)
{
    ui d=vvars[2];
    particles[i].F[d]=F;
    import_force(i,argv...);
}

template<ui dim> void md<dim>::export_pos(ldf *x)
{
    ui d=vvars[3];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) x[i]=particles[i].x[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
    #else
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
    #endif
}

template<ui dim> template<typename...arg> void md<dim>::export_pos(ldf *x,arg...argv)
{
    ui d=vvars[3];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) x[i]=particles[i].x[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
    #else
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
    #endif
    export_pos(argv...);
}

template<ui dim> void md<dim>::export_pos(ui i,ldf &x)
{
    ui d=vvars[3];
    x=particles[i].x[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_pos(ui i,ldf &x,arg...argv)
{
    ui d=vvars[3];
    x=particles[i].x[d];
    export_pos(i,argv...);
}

template<ui dim> void md<dim>::export_vel(ldf *dx)
{
    ui d=vvars[4];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) dx[i]=particles[i].dx[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
    #else
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
    #endif
}

template<ui dim> template<typename...arg> void md<dim>::export_vel(ldf *dx,arg...argv)
{
    ui d=vvars[4];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) dx[i]=particles[i].dx[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
    #else
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
    #endif
    export_vel(argv...);
}

template<ui dim> void md<dim>::export_vel(ui i,ldf &dx)
{
    ui d=vvars[4];
    dx=particles[i].dx[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_vel(ui i,ldf &dx,arg...argv)
{
    ui d=vvars[4];
    dx=particles[i].dx[d];
    export_vel(i,argv...);
}

template<ui dim> void md<dim>::export_force(ldf *F)
{
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) F[i]=particles[i].F[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #else
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #endif
}

template<ui dim> template<typename...arg> void md<dim>::export_force(ldf *F,arg...argv)
{
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) F[i]=particles[i].F[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #else
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #endif
    export_force(argv...);
}

template<ui dim> void md<dim>::export_force(ui i,ldf &F)
{
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    F=particles[i].F[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_force(ui i,ldf &F,arg...argv)
{
    ui d=vvars[5];
    if(!d) calc_forces();
    F=particles[i].F[d];
    export_force(i,argv...);
}

template<ui dim> ldf md<dim>::direct_readout(ui i,uc type)
{
    ui d=vvars[6];
    switch(type)
    {
        case 'v': return particles[i].dx[d]; break;
        case 'F':
        {
            if(avars.export_force_calc) calc_forces();
            return particles[i].F[d];
        }
        break;
        default: return particles[i].x[d]; break;
    }
}

template<ui dim> ldf md<dim>::direct_readout(ui d,ui i,uc type)
{
    switch(type)
    {
        case 'v': return particles[i].dx[d]; break;
        case 'F':
        {
            if(avars.export_force_calc) calc_forces();
            return particles[i].F[d];
        }
        break;
        default: return particles[i].x[d]; break;
    }
}

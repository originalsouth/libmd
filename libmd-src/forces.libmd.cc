#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> void DAMPING(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters,void *sys)
{
    (void) sys;
    (void) particles;
    for(ui d=0;d<dim;d++) p->F[d]+=-parameters->at(0)*p->dx[d];
}

template<ui dim> void DISSIPATION(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters,void *sys)
{
    ldf b = parameters->at(0);
    ui p1=((md<dim>*)sys)->pptrtoui(p);
    for (auto it = particles->begin(); it != particles->end(); it++)
    {
        ui p2=((md<dim>*)sys)->pptrtoui(*it);
        for(ui d=0;d<dim;d++) p->F[d]+=-b*(((md<dim>*)sys)->particles[p1].dx[d]-((md<dim>*)sys)->particles[p2].dx[d]);
    }
}

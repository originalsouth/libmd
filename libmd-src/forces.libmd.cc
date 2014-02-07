#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> void DAMPING(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters)
{
    (void) particles;
    for(ui d=0;d<dim;d++) p->F[d]+=-parameters->at(0)*p->dx[d];
}

template<ui dim> void DISSIPATION(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters)
{
    ldf b = parameters->at(0);
    for (typename vector<particle<dim>*>::iterator it = particles->begin(); it != particles->end(); it++) {
        for(ui d=0;d<dim;d++) p->F[d]+=-b*(p->dx[d]-(*it)->dx[d]);
    }
}

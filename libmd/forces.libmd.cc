#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> void DAMPING(particle<dim> *p,vector<particle<dim>*> particles,vector<ldf> *parameters)
{
    for(ui d=0;d<dim;d++) p->F[d]=-parameters->at(0)*particles->at(0)->dx[d];
}

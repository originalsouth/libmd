#ifndef libmd_h
#include "../libmd.h"
#endif

#define SYS ((md<dim>*)sys)

template<ui dim> void DAMPING(ui i,vector<ui> *particles,vector<ldf> *parameters,void *sys)
{
    (void) particles;
    ldf gamma=parameters->at(0);
    for(ui d=0;d<dim;d++) SYS->particles[i].F[d]+=-gamma*SYS->particles[i].dx[d];
}

template<ui dim> void DISSIPATION(ui i,vector<ui> *particles,vector<ldf> *parameters,void *sys)
{
    ldf b=parameters->at(0);
    for(auto it=particles->begin();it!=particles->end();it++) for(ui d=0;d<dim;d++) {
        SYS->particles[i].F[d]+=b*(SYS->dv(d,i,*it));
    }
}

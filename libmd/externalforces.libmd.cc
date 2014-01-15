#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> externalforces<dim>::externalforces()
{
    add(DAMPING);
}

template<ui dim> ui externalforces<dim>::add(extforceptr<dim> p)
{
    extforces.push_back(p);
    return extforces.size()-1;
}

template<ui dim> void externalforces<dim>::operator()(ui type,ldf force[dim],vector<particle<dim>*> *particles,vector<ldf> *parameters)
{
    extforces[type](force,particles,parameters);
}

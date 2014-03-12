#ifndef libmd_h
#include "../libmd.h"
#endif

forcetype::forcetype(ui noexternalforce,vector<vector<ui>> *plist,vector<ldf> *param)
{
    externalforce=noexternalforce;
    if(plist) particles=*plist;
    else particles.resize(0);
    parameters=*param;
}

template<ui dim> externalforces<dim>::externalforces()
{
    add(DAMPING<dim>);
    add(DISSIPATION<dim>);
}

template<ui dim> ui externalforces<dim>::add(extforceptr<dim> p)
{
    extforces.push_back(p);
    return extforces.size()-1;
}

template<ui dim> void externalforces<dim>::operator()(ui type,particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters)
{
    extforces[type](p,particles,parameters);
}

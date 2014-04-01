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
    extforces.reserve(8);
    add(DAMPING<dim>);
    add(DISSIPATION<dim>);
}

template<ui dim> ui externalforces<dim>::add(extforceptr<dim> p)
{
    extforces.push_back(p);
    return extforces.size()-1;
}

template<ui dim> void externalforces<dim>::operator()(ui type,ui i,vector<ui> *particles,vector<ldf> *parameters,void *sys)
{
    extforces[type](i,particles,parameters,sys);
}

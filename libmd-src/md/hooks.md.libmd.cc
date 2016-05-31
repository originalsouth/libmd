#define __libmd_src_file__
#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_hook(ui hook,std::vector<ldf> &parameters)
{
    hooktype temp(hook,parameters);
    hooks.hookers.push_back(temp);
    return hooks.hookers.size()-1;
}

template<ui dim> bool md<dim>::mod_hook(ui htype,ui hook,std::vector<ldf> &parameters)
{
    if(htype<hooks.hookers.size())
    {
        hooktype temp(hook,parameters);
        hooks.hookers[htype]=temp;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rm_hook(ui htype)
{
    ui pos=hooks.hookers.size();
    if(htype<pos)
    {
        if(htype<pos-1) hooks.hookers[htype]=hooks.hookers.back();
        hooks.hookers.push_back();
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::run_hook(ui htype)
{
    if(htype<hooks.hookers.size())
    {
        hooks.hook(hooks.hookers[htype].hook,hooks.hookers[htype].parameters,this);
        return true;
    }
    else return false;
}

template<ui dim> void md<dim>::run_hooks()
{
    for(hooktype htype:hooks.hookers) hooks.hook(htype.hook,htype.parameters,this);
}

#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> ui t_hook<dim>::add(hookptr<dim> p)
{
    hooks.push_back(p);
    return hooks.size()-1;
}

template<ui dim> void t_hook<dim>::operator()(ui idx,std::vector<ldf> &parameters,void *sys)
{
    hooks[idx](parameters,sys);
}

hooktype::hooktype(ui nohook,std::vector<ldf> &param)
{
    hook=nohook;
    parameters=param;
}

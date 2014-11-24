#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> bcond<dim>::bcond()
{
    bcond_p.reserve(8);
    bcond_x.reserve(8);
    add(BCOND_NONE<dim>,BCOND_NONE<dim>);
    add(BCOND_PERIODIC<dim>,BCOND_PERIODIC<dim>);
    add(BCOND_HARD<dim>,BCOND_HARD<dim>);
    add(BCOND_BOXSHEAR<dim>,BCOND_BOXSHEAR<dim>);
}

template<ui dim> ui bcond<dim>::add(bcondpptr<dim> p,bcondxptr<dim> x)
{
    bcond_p.push_back(p);
    bcond_x.push_back(x);
    return bcond_p.size();
}

template<ui dim> void bcond<dim>::operator()(ui d,ui i,void *sys)
{
    bcond_p[SYS->simbox.bcond[d]](d,i,sys);
}

template<ui dim> void bcond<dim>::operator()(ui d,ldf x[dim],void *sys)
{
    bcond_x[SYS->simbox.bcond[d]](d,x,sys);
}

template<ui dim> void bcond<dim>::operator()(ui k,ui d,ui i,void *sys)
{
    bcond_p[k](d,i,sys);
}

template<ui dim> void bcond<dim>::operator()(ui k,ui d,ldf x[dim],void *sys)
{
    bcond_x[k](d,x,sys);
}

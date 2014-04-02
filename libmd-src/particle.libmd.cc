#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> particle<dim>::particle(ldf mass,ui ptype,bool fixed)
{
    m=mass;
    type=ptype;
    fix=fixed;
    for(ui d=0;d<dim;++d) xsk[d]=numeric_limits<ldf>::infinity();
}

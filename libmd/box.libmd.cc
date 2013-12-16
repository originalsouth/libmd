#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> box<dim>::box()
{
    for(ui d=0;d<dim;d++) bcond[d]=0;
    for(ui d=0;d<dim;d++) L[d]=10.0;
}

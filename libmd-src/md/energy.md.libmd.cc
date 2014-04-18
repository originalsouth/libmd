#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ldf md<dim>::thread_H(ui i)
{
    return thread_T(i)+thread_V(i);
}

template<ui dim> ldf md<dim>::thread_T(ui i)
{
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(particles[i].dx[d],2);
    return 0.5*particles[i].m*retval;
}

template<ui dim> ldf md<dim>::thread_V(ui i)
{
    ldf retval=0.0;
    ldf rcosq=pow(network.rco,2);
    for(auto sij: network.skins[i]) if(i<sij.neighbor)
    {
        const ldf rsq=distsq(i,sij.neighbor);
        if(rsq<rcosq)
        {
            const ldf r=sqrt(rsq);
            retval+=v(network.library[sij.interaction].potential,r,&network.library[sij.interaction].parameters);
            if(network.update) retval-=network.library[sij.interaction].vco;
        }
    }
    return retval;
}

template<ui dim> ldf md<dim>::H()
{
    return T()+V();
}

template<ui dim> ldf md<dim>::T()
{
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=thread_T(i);
    return retval/N;
}

template<ui dim> ldf md<dim>::V()
{
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=thread_V(i);
    return retval/N;
}

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
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i<network.skins[i][j].neighbor)
    {
        const ldf rsq=distsq(i,network.skins[i][j].neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            retval+=(v(network.library[network.skins[i][j].interaction].potential,r,&network.library[network.skins[i][j].interaction].parameters)-network.library[network.skins[i][j].interaction].vco);
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

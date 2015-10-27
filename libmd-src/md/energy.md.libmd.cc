#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ldf md<dim>::thread_H(ui i)
{
    //!
    //! Returns the total energy (kinetic and potential) of particle <tt>i</tt>.
    //!
    return T(i)+V(i);
}

template<ui dim> ldf md<dim>::thread_T(ui i)
{
    //!
    //! Returns the kinetic energy of particle <tt>i</tt>.
    //!
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(particles[i].dx[d],2);
    return 0.5*particles[i].m*retval;
}

template<ui dim> ldf md<dim>::thread_V(ui i,bool higher_index_only)
{
    //!
    //! Returns the potential energy of particle <tt>i</tt>.
    //! If <tt>higher_index_only</tt> is <tt>true</tt>, it only incorporates interactions with particles with a higher index.
    //!
    ldf retval=0.0;
    ldf rcosq=pow(network.rco,2);
    for(auto sij: network.skins[i]) if(!higher_index_only or i<sij.neighbor)
    {
        const ldf rsq=distsq(i,sij.neighbor);
        if(rsq<rcosq)
        {
            const ldf r=sqrt(rsq);
            retval+=v(network.library[sij.interaction].potential,r,network.library[sij.interaction].parameters);
            if(network.update) retval-=network.library[sij.interaction].vco;
        }
    }
    return retval;
}

template<ui dim> ldf md<dim>::H()
{
    //!
    //! Returns the total energy (kinetic and potential) of the system
    //!
    return T()+V();
}

template<ui dim> ldf md<dim>::T()
{
    //!
    //! Returns the total kinetic energy of the system
    //!
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=thread_T(i);
    return retval;
}

template<ui dim> ldf md<dim>::V()
{
    //!
    //! Returns the total potential energy of the system
    //!
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=thread_V(i,true);
    return retval;
}

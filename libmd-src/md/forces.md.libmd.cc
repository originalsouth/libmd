#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_clear_forces(ui i)
{
    //!
    //! Clear the forces variable of particle <tt>i</tt>
    //!
    memset(particles[i].F,0,dim*sizeof(ldf));
}

template<ui dim> void md<dim>::thread_calc_pot_forces(ui i)
{
    //!
    //! This function calculates the forces induced by the potentials acting on particle <tt>i</tt>
    //!
    for(auto sij: network.skins[i]) if(i>sij.neighbor)
    {
        const ldf rsq=distsq(i,sij.neighbor);
        if(!network.update or rsq<pow(get_rco(sij.interaction),2))
        {
            const ldf r=sqrt(rsq);
            DEBUG_3("r = " F_LDF,r);
            const ldf dVdr=v.dr(network.library[sij.interaction].potential,r,&network.library[sij.interaction].parameters);
            DEBUG_3("dV/dr = " F_LDF,dVdr);
            for(ui d=0;d<dim;d++)
            {
                ldf F_i=dd(d,i,sij.neighbor)*dVdr/r;
                particles[i].F[d]+=F_i;
                DEBUG_3("particles[" F_UI "].F[" F_UI "] = " F_LDF,i,d,F_i);
                particles[sij.neighbor].F[d]-=F_i;
                DEBUG_3("particles[" F_UI "].F[" F_UI "] = " F_LDF,sij.neighbor,d,-F_i);
            }
        }
    }
}

template<ui dim> void md<dim>::thread_calc_ext_forces(ui i)
{
    //!
    //! This function calculates the forces induced by external forces acting on particle <tt>i</tt>
    //!
    for(auto ftype: network.forces[i]) f(network.forcelibrary[ftype].externalforce,i,!network.forcelibrary[ftype].particles.empty()?&network.forcelibrary[ftype].particles[i]:nullptr,&network.forcelibrary[ftype].parameters,(md<dim>*)this);
}

template<ui dim> void md<dim>::calc_forces()
{
    //!
    //! This function clears and then calculates all the forces in the system by indexing md<dim>::index if nessecary and then looping over md<dim>::thread_calc_forces.
    //!
    if(network.update and (avars.reindex or test_index()))
    {
        DEBUG_2("regenerating skinlist");
        index();
    }
    DEBUG_2("exec is here");
    avars.export_force_calc=false;
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    recalc_forces();
}

template<ui dim> void md<dim>::recalc_forces()
{
    //!
    //! This function recalculates all the forces in the nessecary for the Velocity Verlet integrator.
    //! Unlike md<dim>::calc_forces this function does not clear nor index.
    //!
    DEBUG_3("exec is here");
    if(!network.library.empty()) for(ui i=0;i<N;i++) thread_calc_pot_forces(i);
    if(!network.forcelibrary.empty()) for(ui i=0;i<N;i++) thread_calc_ext_forces(i);
}

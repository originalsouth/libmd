#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> md<dim>::md(ui particlenr)
{
    //!
    //! Constructor for the md structure.
    //! The default number of particles is zero.
    //! Calls init().
    //!
    init(particlenr);
}

template<ui dim> void md<dim>::init(ui particlenr)
{
    //!
    //! Initialize ::md structure for a given number of particles
    //! specified by \c particlenr. Resizes all lists of structures that
    //! require one element per particle.
    //!
    N=particlenr;
    DEBUG_1("creating md<" F_UI "> with " F_UI " particles",dim,N);
    if(N)
    {
        particles.resize(N);
        network.skins.resize(N);
        network.forces.resize(N);
        network.spid.resize(N);
        for(ui i=0;i<N;i++) network.spid[i]=UI_MAX;
    }
    avars.export_force_calc=true;
}

#include "md/interaction.md.libmd.cc"
#include "md/distances.md.libmd.cc"
#include "md/forces.md.libmd.cc"
#include "md/index.md.libmd.cc"
#include "md/periodicity.md.libmd.cc"
#include "md/integrator.md.libmd.cc"
#include "md/setget.md.libmd.cc"
#include "md/energy.md.libmd.cc"
#include "md/particles.md.libmd.cc"
#include "md/sp.md.libmd.cc"
#include "md/importexport.md.libmd.cc"
#include "md/bonds.md.libmd.cc"

template<ui dim> void md<dim>::clear()
{
    //!
    //! Remove all particles from the ::md structure, and clear all
    //! data types storing interactions and superparticle data. Leaves
    //! the system box and the boundary conditions unchanged.
    //!
    N=0;
    particles.clear();
    network.skins.clear();
    network.library.clear();
    network.lookup.clear();
    network.spid.clear();
    network.superparticles.clear();
    network.sptypes.clear();
    network.forcelibrary.clear();
    network.forces.clear();
    network.free_library_slots.clear();
}

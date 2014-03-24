#ifndef libmd_h
#include "../libmd.h"
#endif


template<ui dim> md<dim>::md(ui particlenr)
{
    //!
    //! Constructor for the md structure.
    //! The default number of particles is zero.
    //! Calls init to actually do the procedure.
    //!
    init(particlenr);
}

template<ui dim> void md<dim>::init(ui particlenr)
{
    N=particlenr;
    DEBUG_1("creating md<%u> with %u particles",dim,N);
    if(N)
    {
        particles.resize(N);
        network.skins.resize(N);
        network.forces.resize(N);
        network.spid.resize(N);
        for(ui i=0;i<N;i++) network.spid[i]=std::numeric_limits<ui>::max();
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





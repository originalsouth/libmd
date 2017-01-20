#define __libmd_src_file__
#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::set_damping(ldf coefficient)
{
    //!
    //! This function adds damping with the given coefficient, or changes the coefficient if damping was already set.
    //!
    if(avars.noftypedamping==UI_MAX)
    {
        DEBUG_2("activating damping with damping coefficient: " F_LDF,coefficient);
        std::vector<ldf> parameters(1,coefficient);
        avars.noftypedamping=add_forcetype(EXTFORCE::DAMPING,parameters);
        assign_all_forcetype(avars.noftypedamping);
    }
    else
    {
        DEBUG_2("modifying damping coefficient to: " F_LDF,coefficient);
        network.forcelibrary[avars.noftypedamping].parameters[0]=coefficient;
    }
}

template<ui dim> bool md<dim>::unset_damping()
{
    //!
    //! This function disables damping and returns whether it was set.
    //!
    if(avars.noftypedamping<UI_MAX)
    {
        DEBUG_2("disabling damping");
        rem_forcetype(avars.noftypedamping);
        avars.noftypedamping=UI_MAX;
        return true;
    }
    else
    {
        WARNING("no damping set");
        return false;
    }
}

template<ui dim> void md<dim>::set_langevin(ldf T,ldf gamma)
{
    //!
    //! This function adds a Langevin thermostat with the given coefficient T and gamma, or changes the coefficient if it was already set.
    //!
    if(avars.noftypelangevin==UI_MAX)
    {
        DEBUG_2("activating langevin thermostat with temperature: " F_LDF " and damping coefficient: " F_LDF,T,gamma);
        std::vector<ldf> parameters={T,gamma};
        avars.noftypelangevin=add_forcetype(EXTFORCE::LANGEVIN,parameters);
        assign_all_forcetype(avars.notypelangevin);
        set_damping(gamma);
    }
    else
    {
        DEBUG_2("modifying langevin coefficient T to: " F_LDF " and gamma to: " F_LDF,T,gamma);
        network.forcelibrary[avars.noftypedamping].parameters[0]=T;
        network.forcelibrary[avars.noftypedamping].parameters[0]=gamma;
        set_damping(gamma);
    }
}

template<ui dim> bool md<dim>::unset_langevin()
{
    //!
    //! This function disables Langevin thermostat and returns whether it was set.
    //!
    if(avars.noftypelangevin<UI_MAX)
    {
        DEBUG_2("disabling Langevin thermostat");
        rem_forcetype(avars.noftypelangevin);
        rem_forcetype(avars.noftypedamping);
        avars.noftypelangevin=UI_MAX;
        avars.noftypedamping=UI_MAX;
        return true;
    }
    else
    {
        WARNING("no Langevin thermostat set");
        return false;
    }
}

template<ui dim> ldf md<dim>::get_rco(ui i,ui j)
{
    auto it=network.lookup.find(network.hash(particles[i].type,particles[j].type));
    if(it!=network.lookup.end()) return get_rco(it->second);
    else return std::numeric_limits<ldf>::quiet_NaN();
}

template<ui dim> ldf md<dim>::get_rco(ui interaction)
{
    return network.library[interaction].rco;
}

template<ui dim> void md<dim>::set_rco(ldf rco)
{
    //!
    //! Sets <tt>network.rco</tt>, the interaction cut-off distance, to <tt>rco</tt>.
    //!
    network.rco=rco;
    for (auto itype: network.library) itype.vco = v(itype.potential,network.rco,itype.parameters);
    if (network.rco > network.ssz)
    {
        WARNING("network.rco is now larger than network.ssz (" F_LDF " > " F_LDF ")",network.rco,network.ssz);
    }
}

template<ui dim> void md<dim>::set_rco(ui interaction,ldf rco)
{
    //!
    //! Sets <tt>network.rco</tt>, the interaction cut-off distance, to <tt>rco</tt>.
    //!
    auto itype=network.library[interaction];
    itype.rco = rco;
    itype.vco = v(itype.potential,rco,itype.parameters);
    if (rco > network.ssz)
    {
        WARNING("this rco is now larger than network.ssz (" F_LDF " > " F_LDF ")",rco,network.ssz);
    }
}

template<ui dim> void md<dim>::set_ssz(ldf ssz)
{
    //!
    //! Sets <tt>network.ssz</tt>, the skinsize, to <tt>ssz</tt> and reserves space for <tt>network.skins[]</tt>.
    //!
    network.ssz=ssz;
    set_reserve(ssz);
    for(auto itype:network.library) if(itype.rco>network.ssz)
    {
        WARNING("rco of interaction %zu is now larger than network.ssz (" F_LDF " > " F_LDF ")",&itype-&network.library[0],itype.rco,network.ssz);
    }
    if(network.rco>network.ssz)
    {
        WARNING("network.rco is now larger than network.ssz (" F_LDF " > " F_LDF ")",network.rco,network.ssz);
    }
}

template<ui dim> void md<dim>::set_reserve(ldf ssz)
{
    //!
    //! This function reserves space for the vectors in <tt>network.skins[]</tt>,
    //! assuming that (approximately) all particles interact with each other.
    //!
    set_reserve(ssz,N);
}

template<ui dim> void md<dim>::set_reserve(ldf ssz,ui M)
{
    //!
    //! This function reserves space for the vectors in <tt>network.skins[]</tt>,
    //! the amount of space being based on the assumption that each particle has at most <tt>M</tt> particles it interacts with
    //! and they are uniformly distributed over the system.
    //!
    ldf area=1.0;
    for(ui d=0;d<dim;d++) area*=simbox.L[d];
    const ldf vol=(std::pow(std::acos(-1.0),((ldf)dim)/2.0)/std::tgamma(1.0+((ldf)dim)/2.0))*std::pow(ssz,dim);
    const ui fed=std::min(N-1,(ui)(((ldf)M)*(vol)/(area))*2+4);
    DEBUG_1("reserved skin size set to " F_UI " skins",fed);
    for(ui i=0;i<N;i++) network.skins[i].reserve(fed);
}

template<ui dim> void md<dim>::set_type(ui p, ui newtype)
{
    //!
    //! Sets the type of particle <tt>p</tt> to <tt>newtype</tt>.
    //!
    particles[p].type=newtype;
}

template<ui dim> void md<dim>::set_index_method(ui method)
{
    //!
    //! Sets the indexing method to <tt>method</tt>.
    //!
    indexdata.method=method;
}

template<ui dim> void md<dim>::set_bcond(uc bcond[dim])
{
    //!
    //! Sets the global boundary conditions
    //!
    memcpy(simbox.bcond,bcond,dim*sizeof(uc));
}

template<ui dim> void md<dim>::set_pbcond(ui i,uc bcond[dim],bool toggle)
{
    //!
    //! Sets particle's boundary conditions
    //!
    particles[i].usepbcond=toggle;
    memcpy(particles[i].pbc,bcond,dim*sizeof(uc));
}

template<ui dim> void md<dim>::set_spbcond(ui spi,uc bcond[dim],bool toggle)
{
    //!
    //! Sets particle's boundary conditions for the superparticles
    //!
    ui K=network.superparticles[spi].backdoor.size(),i;
    for(ui k=0;k<K;k++)
    {
        if((i=network.superparticles[spi].backdoor[k])<UI_MAX)
        {
            particles[i].usepbcond=toggle;
            memcpy(particles[i].pbc,bcond,dim*sizeof(uc));
        }
    }
}

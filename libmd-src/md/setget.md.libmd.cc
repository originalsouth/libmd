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
        vector<ldf> parameters(1,coefficient);
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
        return true;
    }
    else
    {
        WARNING("no damping set");
        return false;
    }
}

template<ui dim> ldf md<dim>::get_rco(ui i,ui j)
{
    auto it=network.lookup.find(network.hash(particles[i].type,particles[j].type));
    if(it!=network.lookup.end()) return get_rco(it->second);
    else return numeric_limits<ldf>::quiet_NaN();
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
        WARNING("rco of interaction %zu is now larger than network.ssz (" F_LDF " > " F_LDF ")",itype-&network.library[0],itype.rco,network.ssz);
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
    const ldf vol=(pow(acos(-1.0),((ldf)dim)/2.0)/tgamma(1.0+((ldf)dim)/2.0))*pow(ssz,dim);
    const ui fed=min(N-1,(ui)(((ldf)M)*(vol)/(area))*2+4);
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

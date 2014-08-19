#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::interactions(ui i,vector<pair<ui,ui>> &table)
{
    //!
    //! This function puts a list of all interaction neighbors of particle <tt>i</tt> in <tt>table</tt>,
    //! as pairs of particles with <tt>i</tt> being the first.
    //!
    table.clear();
    ldf rcosq=pow(network.rco,2);
    for(auto sij: network.skins[i]) if(distsq(i,sij.neighbor)<rcosq) table.push_back(pair<ui,ui>(i,sij.neighbor));
}

template<ui dim> void md<dim>::all_interactions(vector<pair<ui,ui>> &table)
{
    //!
    //! This function puts a list of all interaction neighbors in <tt>table</tt>, as pairs of particles.
    //!
    table.clear();
    ldf rcosq=pow(network.rco,2);
    for(ui i=0;i<N;i++) for(auto sij: network.skins[i]) if(i>sij.neighbor and distsq(i,sij.neighbor)<rcosq) table.push_back(pair<ui,ui>(i,sij.neighbor));
}

template<ui dim> ui md<dim>::add_interaction(ui potential,vector<ldf> *parameters)
{
    //!
    //! This function adds a new interaction, of the given type and with the given parameters, to <tt>network.library[]</tt> and returns its index.
    //!
    interactiontype itype(potential,parameters,v(potential,network.rco,parameters));
    if(network.free_library_slots.empty())
    {
        network.library.push_back(itype);
        return network.library.size()-1;
    }
    else
    {
        ui retval=*network.free_library_slots.begin();
        network.free_library_slots.erase(network.free_library_slots.begin());
        network.library[retval]=itype;
        return retval;
    }
}

template<ui dim> bool md<dim>::mod_interaction(ui interaction,ui potential,vector<ldf> *parameters)
{
    //!
    //! This function replaces the interaction in <tt>network.library[]</tt> with index <tt>interaction</tt>
    //! with a potential of the given type and with the given parameters.
    //! It returns whether the given library element exists.
    //!
    if(interaction>=network.library.size())
    {
        WARNING("Interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("Interaction %d was previously removed", interaction);
        return false;
    }
    else
    {
        interactiontype itype(potential,parameters,v(potential,network.rco,parameters));
        network.library[interaction]=itype;
        avars.reindex=true;
        return true;
    }
}

template<ui dim> bool md<dim>::rem_interaction(ui interaction)
{
    //!
    //! This function removes the interaction with index <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! It returns whether the given library element existed.
    //!
    if(interaction>=network.library.size())
    {
        WARNING("Interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("Interaction %d was previously removed", interaction);
        return false;
    }
    else
    {
        network.free_library_slots.insert(interaction);
        for(auto it=network.lookup.begin();it!=network.lookup.end();) (it->second==interaction)?network.lookup.erase(it++):it++;
        for(ui i=network.sptypes.size()-1;i<UI_MAX;i--) for(auto it=network.sptypes[i].splookup.begin();it!=network.sptypes[i].splookup.end();) (it->second==interaction)?network.sptypes[i].splookup.erase(it++):it++;
        avars.reindex=true;
        return true;
    }
}

template<ui dim> bool md<dim>::add_typeinteraction(ui type1,ui type2,ui interaction)
{
    //!
    //! This function adds a type interaction between the given types, using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! It returns whether the given interaction exists and the two types do not already have an interaction.
    //! If the two types already have an interaction, it is not modified.
    //!
    pair<ui,ui> id=network.hash(type1,type2);
    if(network.lookup.count(id))
        return false;
    else if(interaction>=network.library.size())
    {
        WARNING("Interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("Interaction %d was previously removed", interaction);
        return false;
    }
    else
    {
        network.lookup[id]=interaction;
        avars.reindex=true;
        return true;
    }
}

template<ui dim> bool md<dim>::mod_typeinteraction(ui type1,ui type2,ui interaction)
{
    //!
    //! This function modifies the type interaction between the given types, using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! It returns whether the given interaction exists and the two types already had an interaction.
    //! If the two types do not already have an interaction, it is not added.
    //!
    pair<ui,ui> id=network.hash(type1,type2);
    if(!network.lookup.count(id))
        return false;
    else if(interaction>=network.library.size())
    {
        WARNING("Interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("Interaction %d was previously removed", interaction);
        return false;
    }
    else
    {
        network.lookup[id]=interaction;
        avars.reindex=true;
        return true;
    }
}

template<ui dim> void md<dim>::mad_typeinteraction(ui type1,ui type2,ui interaction)
{
    //!
    //! This function assigns a type interaction to the given pair of types, using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! It does not perform any checks.
    //!
    network.lookup[network.hash(type1,type2)]=interaction;
    avars.reindex=true;
}

template<ui dim> bool md<dim>::add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    //!
    //! This function adds a type interaction between the given types, using a new interaction of the given type and with the given parameters.
    //! It returns whether the two types do not already have an interaction.
    //! If the two types already have an interaction, it is not modified and the interaction is not added to <tt>network.library[]</tt>.
    //!
    pair<ui,ui> id=network.hash(type1,type2);
    if(!network.lookup.count(id))
    {
        network.lookup[id]=add_interaction(potential,parameters);
        avars.reindex=true;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    //!
    //! This function modifies the type interaction between the given types, using a new interaction of the given type and with the given parameters.
    //! It returns whether the two types already had an interaction.
    //! If the two types do not already have an interaction, it is not added and the interaction is not added to <tt>network.library[]</tt>.
    //!
    pair<ui,ui> id=network.hash(type1,type2);
    if(!network.lookup.count(id)) return false;
    else
    {
        network.lookup[id]=add_interaction(potential,parameters);
        avars.reindex=true;
        return true;
    }
}

template<ui dim> void md<dim>::mad_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    //!
    //! This function assigns a type interaction to the given pair of types, using a new interaction of the given type and with the given parameters.
    //! It does not perform any checks.
    //!
    network.lookup[network.hash(type1,type2)]=add_interaction(potential,parameters);
    avars.reindex=true;
}

template<ui dim> bool md<dim>::rem_typeinteraction(ui type1,ui type2)
{
    //!
    //! This function removes the type interaction between the given types. It does not remove the interaction from<tt>network.library[]</tt>.
    //! It returns whether the two types had an interaction.
    //!
    avars.reindex=true;
    return network.lookup.erase(network.hash(type1,type2));
}


/*** Superparticle interactions ***/

template<ui dim> ui md<dim>::add_sptype()
{
    //!
    //! This function adds an empty superparticletype to <tt>network.sptypes[]</tt> and returns its index.
    //!
    network.sptypes.push_back(superparticletype());
    return network.sptypes.size()-1;
}

template<ui dim> bool md<dim>::rem_sptype(ui spt)
{
    //!
    //! This function removes element <tt>spt</tt> from <tt>network.sptypes[]</tt>.
    //! It returns whether this element existed.
    //!
    if(spt<network.sptypes.size())
    {
        ui spn=network.sptypes.size()-1;
        for(ui i=network.superparticles.size()-1;i<UI_MAX;i--)
        {
            if(network.superparticles[i].sptype==spt) network.superparticles[i].sptype=UI_MAX;
            if(network.superparticles[i].sptype==spn) network.superparticles[i].sptype=spt;
        }
        iter_swap(network.superparticles.begin()+spt,network.superparticles.rbegin());
        network.sptypes.pop_back();
        avars.reindex=true;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::add_sp_interaction(ui spt,ui p1,ui p2,ui interaction)
{
    //!
    //! This function adds an interaction between particle numbers <tt>p1</tt> and <tt>p2</tt> of superparticles of type <tt>spt</tt>,
    //! using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! It returns whether the given superparticletype and interaction exist and the particle numbers do not already have an interaction.
    //! If the two particle numbers already have an interaction, it is not modified.
    //!
    pair<ui,ui> id=network.hash(p1,p2);
    if(spt>=network.sptypes.size())
    {
        WARNING("Superparticletype %d does not exist", spt);
        return false;
    }
    else if(network.sptypes[spt].splookup.count(id))
        return false;
    else if(interaction>=network.library.size())
    {
        WARNING("Interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("Interaction %d was previously removed", interaction);
        return false;
    }
    else
    {
        network.sptypes[spt].splookup[id]=interaction;
        avars.reindex=true;
        return true;
    }
}

template<ui dim> bool md<dim>::mod_sp_interaction(ui spt,ui p1,ui p2,ui interaction)
{
    //!
    //! This function modifies the interaction between particle numbers <tt>p1</tt> and <tt>p2</tt> of superparticles of type <tt>spt</tt>,
    //! using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! It returns whether the given superparticletype and interaction exist and the particle numbers already had an interaction.
    //! If the two particle numbers do not already have an interaction, it is not added.
    //!
    pair<ui,ui> id=network.hash(p1,p2);
    if(spt>=network.sptypes.size())
    {
        WARNING("Superparticletype %d does not exist", spt);
        return false;
    }
    else if(!network.sptypes[spt].splookup.count(id))
        return false;
    else if(interaction>=network.library.size())
    {
        WARNING("Interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("Interaction %d was previously removed", interaction);
        return false;
    }
    else
    {
        network.sptypes[spt].splookup[id]=interaction;
        avars.reindex=true;
        return true;
    }
}

template<ui dim> ui md<dim>::mad_sp_interaction(ui spt,ui p1,ui p2,ui interaction)
{
    //!
    //! This function assigns an interaction to particle numbers <tt>p1</tt> and <tt>p2</tt> of superparticles of type <tt>spt</tt>,
    //! using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! If the given superparticletype does not exist, a new one is added to <tt>network.sptypes[]</tt>.
    //! It returns the index of the superparticletype (which is <tt>spt</tt> or the index of the newly created one).
    //! It does not perform any checks.
    //!
    if (spt >= network.sptypes.size())
    {   spt = network.sptypes.size();
        network.sptypes.push_back(superparticletype());
    }
    network.sptypes[spt].splookup[network.hash(p1,p2)] = interaction;
    avars.reindex=true;
    return spt;
}

template<ui dim> bool md<dim>::add_sp_interaction(ui spt,ui p1,ui p2,ui potential,vector<ldf> *parameters)
{
    //!
    //! This function adds an interaction between particle numbers <tt>p1</tt> and <tt>p2</tt> of superparticles of type <tt>spt</tt>,
    //! using a new interaction of the given type and with the given parameters.
    //! It returns whether the given superparticletype exists and the particle numbers do not already have an interaction.
    //! If the two particle numbers already have an interaction,
    //! it is not modified and the interaction is not added to <tt>network.library[]</tt>.
    //!
    pair<ui,ui> id=network.hash(p1,p2);
    if(spt>=network.sptypes.size())
    {
        WARNING("Superparticletype %d does not exist", spt);
        return false;
    }
    else if(network.sptypes[spt].splookup.count(id))
        return false;
    else
    {
        network.sptypes[spt].splookup[id] = add_interaction(potential,parameters);
        avars.reindex=true;
        return true;
    }
}

template<ui dim> bool md<dim>::mod_sp_interaction(ui spt,ui p1,ui p2,ui potential,vector<ldf> *parameters)
{
    //!
    //! This function modifies the interaction between particle numbers <tt>p1</tt> and <tt>p2</tt> of superparticles of type <tt>spt</tt>,
    //! using a new interaction of the given type and with the given parameters.
    //! It returns whether the given superparticletype and interaction exist and the particle numbers already had an interaction.
    //! If the two particle numbers do not already have an interaction,
    //! it is not added and the interaction is not added to <tt>network.library[]</tt>.
    //!
    pair<ui,ui> id=network.hash(p1,p2);
    if(spt>=network.sptypes.size())
    {
        WARNING("Superparticletype %d does not exist", spt);
        return false;
    }
    else if(!network.sptypes[spt].splookup.count(id))
        return false;
    else
    {
        network.sptypes[spt].splookup[id] = add_interaction(potential,parameters);
        avars.reindex=true;
        return true;
    }
}

template<ui dim> ui md<dim>::mad_sp_interaction(ui spt,ui p1,ui p2,ui potential,vector<ldf> *parameters)
{
    //!
    //! This function assigns an interaction to particle numbers <tt>p1</tt> and <tt>p2</tt> of superparticles of type <tt>spt</tt>,
    //! using a new interaction of the given type and with the given parameters.
    //! If the given superparticletype does not exist, a new one is added to <tt>network.sptypes[]</tt>.
    //! It returns the index of the superparticletype (which is <tt>spt</tt> or the index of the newly created one).
    //! It does not perform any checks.
    //!
    if (spt >= network.sptypes.size())
    {
        spt = network.sptypes.size();
        network.sptypes.push_back(superparticletype());
    }
    network.sptypes[spt].splookup[network.hash(p1,p2)] = add_interaction(potential,parameters);
    avars.reindex=true;
    return spt;
}

template<ui dim> bool md<dim>::rem_sp_interaction(ui spt,ui p1,ui p2)
{
    //!
    //! This function removes the interaction between particle numbers <tt>p1</tt> and <tt>p2</tt> of superparticles of type <tt>spt</tt>.
    //! It does not remove the interaction from<tt>network.library[]</tt>.
    //! It returns whether the given superparticletype exists and the two particle numbers had an interaction.
    //!
    if(spt>=network.sptypes.size())
    {
        WARNING("Superparticletype %d does not exist", spt);
        return false;
    }
    else
    {
        avars.reindex=true;
        return network.sptypes[spt].splookup.erase(network.hash(p1,p2));
    }
}


/*** Forcetypes ***/

template<ui dim> ui md<dim>::add_forcetype(ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters)
{
    //!
    //! This function adds a new forcetype, of the given type and with the given parameters, to <tt>network.forcelibrary[]</tt> and returns its index.
    //!
    forcetype temp(force,noparticles,parameters);
    network.forcelibrary.push_back(temp);
    return network.forcelibrary.size()-1;
}

template<ui dim> bool md<dim>::mod_forcetype(ui ftype,ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters)
{
    //!
    //! This function replaces the forcetype in <tt>network.forcelibrary[]</tt> with index <tt>ftype</tt>
    //! with a forcetype of the given type and with the given parameters.
    //! It returns whether the given forcelibrary element exists.
    //!
    if(ftype<network.forcelibrary.size())
    {
        forcetype temp(force,noparticles,parameters);
        network.forcelibrary[ftype]=temp;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_forcetype(ui ftype)
{
    //!
    //! This function removes the forcetype with index <tt>ftype</tt> from <tt>network.forcelibrary[]</tt>.
    //! It returns whether the given forcelibrary element existed.
    //!
    ui pos=network.forcelibrary.size();
    if(ftype<pos)
    {
        if(ftype<pos-1)
        {
            network.forcelibrary[ftype]=network.forcelibrary[pos-1];
            for(ui i=0;i<N;i++) for(ui j=network.forces[i].size()-1;j<UI_MAX;j--) if(network.forces[i][j]==pos-1) network.forces[i][j]=ftype; //FIXME
        }
        network.forcelibrary.pop_back();
        for(ui i=0;i<N;i++) for(ui j=network.forces[i].size()-1;j<UI_MAX;j--) if(network.forces[i][j]==ftype)
        {
            network.forces[i][j]=network.forces[i].back();
            network.forces[i].pop_back();
            break;
        }
    }
    else return false;
}

template<ui dim> bool md<dim>::assign_forcetype(ui i,ui ftype)
{
    //!
    //! This function assigns a forcetype to particle <tt>i</tt>, using element <tt>ftype</tt> from <tt>network.forcelibrary[]</tt>.
    //! It returns whether the given forcetype was not already assigned to the given particle.
    //!
    for(ui j=network.forces[i].size()-1;j<UI_MAX;j--) if(network.forces[i][j]==ftype) return false;
    network.forces[i].push_back(ftype);
    return true;
}

template<ui dim> void md<dim>::assign_all_forcetype(ui ftype)
{
    //!
    //! This function assigns the forcetype <tt>network.forcelibrary[ftype]</tt> to all particles.
    //!
    for(ui i=0;i<N;i++) assign_forcetype(i,ftype);
}

template<ui dim> bool md<dim>::unassign_forcetype(ui i,ui ftype)
{
    //!
    //! This function removes the forcetype <tt>network.forcelibrary[ftype]</tt> from particle <tt>i</tt>.
    //! It returns whether the given forcetype was assigned to the given particle.
    //!
    for(ui j=network.forces[i].size()-1;j<UI_MAX;j--) if(network.forces[i][j]==ftype)
    {
        network.forces[i][j]=network.forces[i].back();
        network.forces[i].pop_back();
        return true;
    }
    return false;
}

template<ui dim> void md<dim>::unassign_all_forcetype(ui ftype)
{
    //!
    //! This function removes the forcetype <tt>network.forcelibrary[ftype]</tt> from all particles.
    //!
    for(ui i=0;i<N;i++) unassign_forcetype(i,ftype);
}

template<ui dim> void md<dim>::clear_all_assigned_forcetype()
{
    //!
    //! This function removes all forcetypes from all particles
    //!
    for(ui i=0;i<N;i++) network.forces[i].clear();
}

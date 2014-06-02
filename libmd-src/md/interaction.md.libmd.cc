#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::all_interactions(vector<pair<ui,ui>> &table)
{
    table.clear();
    ldf rcosq=pow(network.rco,2);
    for(ui i=0;i<N;i++) for(auto sij: network.skins[i]) if(i>sij.neighbor and distsq(i,sij.neighbor)<rcosq) table.push_back(pair<ui,ui>(i,sij.neighbor));
}

template<ui dim> ui md<dim>::add_interaction(ui potential,vector<ldf> *parameters)
{
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
    network.lookup[network.hash(type1,type2)]=interaction;
    avars.reindex=true;
}

template<ui dim> bool md<dim>::add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
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
    network.lookup[network.hash(type1,type2)]=add_interaction(potential,parameters);
    avars.reindex=true;
}

template<ui dim> bool md<dim>::rem_typeinteraction(ui type1,ui type2)
{
    avars.reindex=true;
    return network.lookup.erase(network.hash(type1,type2));
}


/*** Superparticle interactions ***/

template<ui dim> ui md<dim>::add_sptype()
{
    network.sptypes.push_back(superparticletype());
    return network.sptypes.size()-1;
}

template<ui dim> bool md<dim>::rem_sptype(ui spt)
{
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
    forcetype temp(force,noparticles,parameters);
    network.forcelibrary.push_back(temp);
    return network.forcelibrary.size()-1;
}

template<ui dim> bool md<dim>::mod_forcetype(ui ftype,ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters)
{
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
    ui pos=network.forcelibrary.size();
    if(ftype<pos)
    {
        if(ftype<pos-1)
        {
            network.forcelibrary[ftype]=network.forcelibrary[pos-1];
            for(ui i=0;i<N;i++) for(ui j=network.forces[i].size()-1;j<UI_MAX;j--) if(network.forces[i][j]==pos-1) network.forces[i][j]=ftype;
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
    for(ui j=network.forces[i].size()-1;j<UI_MAX;j--) if(network.forces[i][j]==ftype) return false;
    network.forces[i].push_back(ftype);
    return true;
}

template<ui dim> void md<dim>::assign_all_forcetype(ui ftype)
{
    for(ui i=0;i<N;i++) assign_forcetype(i,ftype);
}

template<ui dim> void md<dim>::unassign_forcetype(ui i,ui ftype)
{
    for(ui j=network.forces[i].size()-1;j<UI_MAX;j--) if(network.forces[i][j]==ftype)
    {
        network.forces[i][j]=network.forces[i].back();
        network.forces[i].pop_back();
        break;
    }
}

template<ui dim> void md<dim>::unassign_all_forcetype(ui ftype)
{
    for(ui i=0;i<N;i++) unassign_forcetype(i,ftype);
}

template<ui dim> void md<dim>::clear_all_assigned_forcetype()
{
    for(ui i=0;i<N;i++) network.forces[i].clear();
}

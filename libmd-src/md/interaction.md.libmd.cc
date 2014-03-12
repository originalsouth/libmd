#ifndef libmd_h
#include "../../libmd.h"
#endif

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
        ui retval=network.free_library_slots.top();
        network.free_library_slots.pop();
        network.library[retval]=itype;
        return retval;
    }
}

template<ui dim> bool md<dim>::mod_interaction(ui interaction,ui potential,vector<ldf> *parameters)
{
    if(interaction<network.library.size())
    {
        interactiontype itype(potential,parameters,v(potential,network.rco,parameters));
        network.library[interaction]=itype;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_interaction(ui interaction)
{
    if(interaction<network.library.size())
    {
        network.free_library_slots.push(interaction);
        for(auto it=network.lookup.begin();it!=network.lookup.end();) (it->second==interaction)?network.lookup.erase(it++):it++;
        ui M=network.sptypes.size();
        for(ui i=0;i<M;i++) for(auto it=network.sptypes[i].splookup.begin();it!=network.sptypes[i].splookup.end();) (it->second==interaction)?network.sptypes[i].splookup.erase(it++):it++;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::add_typeinteraction(ui type1,ui type2,ui interaction)
{
    pair<ui,ui> id=network.hash(type1,type2);
    if(!network.lookup.count(id))
    {
        network.lookup[id]=interaction;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::mod_typeinteraction(ui type1,ui type2,ui interaction)
{
    pair<ui,ui> id=network.hash(type1,type2);
    if(!network.lookup.count(id)) return false;
    else
    {
        network.lookup[id]=interaction;
        return true;
    }
}

template<ui dim> void md<dim>::mad_typeinteraction(ui type1,ui type2,ui interaction)
{
    network.lookup[network.hash(type1,type2)]=interaction;
}

template<ui dim> ui md<dim>::add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    if(!network.lookup.count(id))
    {
        network.lookup[id]=add_interaction(potential,parameters);
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
        return true;
    }
}

template<ui dim> bool md<dim>::rem_typeinteraction(ui type1,ui type2)
{
    return network.lookup.erase(network.hash(type1,type2));
}

template<ui dim> ui md<dim>::add_sp_interaction(ui spt,ui p1,ui p2,ui interaction)
{
    if(spt<network.sptypes.size())
    {
        pair<ui,ui> id=network.hash(p1,p2);
        network.sptypes[spt].splookup[id]=interaction;
        return spt;
    }
    else
    {
        spt=network.sptypes.size();
        superparticletype sptype;
        pair<ui,ui> id=network.hash(p1,p2);
        sptype.splookup[id]=interaction;
        network.sptypes.push_back(sptype);
        return spt;
    }
}

template<ui dim> bool md<dim>::mod_sp_interaction(ui spt,ui p1,ui p2,ui interaction)
{
    if(spt<network.sptypes.size())
    {
        add_sp_interaction(spt,p1,p2,interaction);
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_sp_interaction(ui spt,ui p1,ui p2)
{
    if(spt<network.sptypes.size())
    {
        pair<ui,ui> id=network.hash(p1,p2);
        if(network.sptypes[spt].splookup.find(id)!=network.sptypes[spt].splookup.end())
        {
            network.sptypes[spt].splookup.erase(id);
            return true;
        }
        else return false;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_sp_interaction(ui spt)
{
    if(spt<network.sptypes.size())
    {
        ui spn=network.sptypes.size()-1;
        for(ui i=network.superparticles.size();i<numeric_limits<ui>::max();i--)
        {
            if(network.superparticles[i].sptype==spt) network.superparticles[i].sptype=numeric_limits<ui>::max();
            if(network.superparticles[i].sptype==spn) network.superparticles[i].sptype=spt;
        }
        iter_swap(network.superparticles.begin()+spt,network.superparticles.end());
        network.sptypes.pop_back();
        return true;
    }
    else return false;
}

template<ui dim> ui md<dim>::add_forcetype(ui force,vector<ui> *noparticles,vector<ldf> *parameters)
{
    forcetype temp(force,noparticles,parameters);
    network.forcelibrary.push_back(temp);
    return network.forcelibrary.size()-1;
}

template<ui dim> bool md<dim>::mod_forcetype(ui notype,ui force,vector<ui> *noparticles,vector<ldf> *parameters)
{
    if(notype<network.forcelibrary.size())
    {
        forcetype temp(force,noparticles,parameters);
        network.forcelibrary[notype]=temp;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_forcetype(ui notype)
{
    ui pos=network.forcelibrary.size();
    if(notype<pos)
    {
        if(notype==pos-1)
        {
             network.forcelibrary.erase(network.forcelibrary.begin()+notype);
             for(ui i=0;i<N;i++) for(ui j=network.forces[i].size()-1;j<numeric_limits<ui>::max();j--) if(network.forces[i][j]==notype) network.forces[i].erase(network.forces[i].begin()+j);
        }
        else
        {
            network.forcelibrary[notype]=network.forcelibrary[pos-1];
            network.forcelibrary.erase(network.forcelibrary.begin()+pos-1);
            for(ui i=0;i<N;i++) for(ui j=network.forces[i].size()-1;j<numeric_limits<ui>::max();j--)
            {
                if(network.forces[i][j]==notype) network.forces[i].erase(network.forces[i].begin()+j);
                if(network.forces[i][j]==pos-1) network.forces[i][j]=notype;
            }
        }
    }
    else return false;
}

template<ui dim> void md<dim>::assign_forcetype(ui particlenr,ui ftype)
{
    network.forces[particlenr].push_back(ftype);
}

template<ui dim> void md<dim>::assign_all_forcetype(ui ftype)
{
    for(ui i=0;i<N;i++) network.forces[i].push_back(ftype);
}

template<ui dim> void md<dim>::unassign_forcetype(ui particlenr,ui ftype)
{
    for(ui i=network.forces[particlenr].size()-1;i<numeric_limits<ui>::max();i--) if(network.forces[particlenr][i]==ftype)
    {
        network.forces[particlenr].erase(network.forces[particlenr].begin()+i);
        break;
    }
}

template<ui dim> void md<dim>::unassign_all_forcetype(ui ftype)
{
    for(ui j=0;j<N;j++) for(ui i=network.forces[j].size()-1;i<numeric_limits<ui>::max();i++) if(network.forces[j][i]==ftype)
    {
        network.forces[j].erase(network.forces[j].begin()+i);
        break;
    }
}

template<ui dim> void md<dim>::clear_all_assigned_forcetype()
{
    for(ui i=0;i<N;i++) network.forces[i].clear();
}
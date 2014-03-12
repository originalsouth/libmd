#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_interaction(ui potential,vector<ldf> *parameters)
{
    return numeric_limits<ui>::max();
}

template<ui dim> bool md<dim>::mod_interaction(ui interaction,ui potential,vector<ldf> *parameters)
{
    return true;
}

template<ui dim> bool md<dim>::rem_interaction(ui interaction)
{
    return true;
}

template<ui dim> ui md<dim>::add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,parameters));
    if(network.lookup.find(id)==network.lookup.end())
    {
        network.library.push_back(itype);
        network.lookup[id]=network.library.size()-1;
        network.backdoor.push_back(id);
        return network.library.size()-1;
    }
    else return numeric_limits<ui>::max();
}

template<ui dim> bool md<dim>::mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,network.rcosq,parameters));
    if(network.lookup.find(id)==network.lookup.end()) return false;
    else
    {
        network.library[network.lookup[id]]=itype;
        return true;
    }
}

template<ui dim> bool md<dim>::rem_typeinteraction(ui type1,ui type2)
{
    pair<ui,ui> id=network.hash(type1,type2);
    if(network.lookup.find(id)!=network.lookup.end())
    {
        ui pos=network.lookup[id];
        network.library[pos]=network.library.back();
        network.backdoor[pos]=network.backdoor.back();
        network.lookup[network.backdoor[pos]]=pos;
        network.library.pop_back();
        network.backdoor.pop_back();
        network.lookup.erase(id);
        return true;
    }
    else return false;
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

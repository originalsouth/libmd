#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_sp(ui sptype)
{
    ui spi = network.superparticles.size();
    network.superparticles.push_back(superparticle());
    network.superparticles[spi].sptype = sptype;
    DEBUG_2("added superparticle #%u with type %u",spi,sptype);
    return spi;
}

template<ui dim> bool md<dim>::rem_sp(ui spi)
{
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle %d does not exist", spi);
        return false;
    }
    else
    {
        ui spn=network.superparticles.size()-1;
        if(spi<spn)
        {
            for(auto m: network.superparticles[spn].particles) network.spid[m.first]=spi;
            iter_swap(network.superparticles.begin()+spi,network.superparticles.rbegin());
        }
        for(auto m: network.superparticles[spn].particles) network.spid[m.first]=UI_MAX;
        network.superparticles.pop_back();
        DEBUG_2("removed superparticle #%u",spi);
        return true;
    }
}

template<ui dim> bool md<dim>::rem_sp_particles(ui spi)
{
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle %d does not exist", spi);
        return false;
    }
    else
    {
        for(auto m: network.superparticles[spi].particles) rem_particle(m.first);
        network.superparticles[spi].particles.clear();
        network.superparticles[spi].backdoor.clear();
        DEBUG_2("removed particles of superparticle #%u",spi);
        return true;
    }
}

template<ui dim> ui md<dim>::sp_ingest(ui spi,ui i,ui idx)
{
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle %d does not exist", spi);
        return UI_MAX;
    }
    else if(network.spid[i] < UI_MAX)
    {
        WARNING("particle %u is already in superparticle %u",i,network.spid[i]);
        return UI_MAX;
    }
    else
    {   
        if(idx>=network.superparticles[spi].backdoor.size())
        {
            if(idx==UI_MAX) idx=network.superparticles[spi].backdoor.size();
            network.superparticles[spi].backdoor.resize(idx+1);
        }
        else if(network.superparticles[spi].backdoor[idx]<UI_MAX)
        {
            WARNING("superparticle #%u already contains index #%u",spi,idx);
            return UI_MAX;
        }
        network.spid[i]=spi;
        network.superparticles[spi].backdoor[idx]=i;
        DEBUG_2("particle #%u is ingested by superparticle #%u",i,spi);
        return network.superparticles[spi].particles[i]=idx;
    }
}

template<ui dim> bool md<dim>::sp_dispose(ui i)
{
    ui spi=network.spid[i];
    if(spi==UI_MAX) return false;
    else
    {
        network.spid[i]=UI_MAX;
        DEBUG_2("particle #%u is removed from superparticle #%u",i,spi);
        ui j=network.superparticles[spi].particles[i];
        if(j==network.superparticles[spi].backdoor.size()-1)
        {
            network.superparticles[spi].backdoor.pop_back();
            while(!network.superparticles[spi].backdoor.empty() and network.superparticles[spi].backdoor.back()==UI_MAX) network.superparticles[spi].backdoor.pop_back();
        }
        else network.superparticles[spi].backdoor[j]=UI_MAX;
        return network.superparticles[spi].particles.erase(i);
    }
}

template<ui dim> bool md<dim>::sp_dispose_idx(ui spi,ui idx)
{
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle #%u does not exist",spi);
        return false;
    }
    else if(idx>=network.superparticles[spi].backdoor.size())
    {
        WARNING("superparticle #%u does not contain index #%u",spi,idx);
        return false;
    }
    else
    {
        ui i=network.superparticles[spi].backdoor[idx];
        if(i<UI_MAX)
        {
            sp_dispose(i);
            DEBUG_2("particle #%u is removed from superparticle #%u",i,spi);
            return true;
        }
        else
        {
            WARNING("superparticle #%u does not contain index #%u",spi,idx);
            return false;
        }
    }
}

template<ui dim> ui md<dim>::sp_pid(ui spi,ui idx)
{
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle #%u does not exist",spi);
        return UI_MAX;
    }
    else if(idx>=network.superparticles[spi].backdoor.size())
    {
        WARNING("superparticle #%u does not contain index #%u",spi,idx);
        return UI_MAX;
    }
    else return network.superparticles[spi].backdoor[idx];
}

template<ui dim> void md<dim>::fix_sp(ui spi,bool fix)
{
    DEBUG_2("fixing(%d) super particle particle %u.",fix,spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) particles[it->first].fix=fix;
}

template<ui dim> ui md<dim>::clone_sp(ui spi,ldf x[dim])
{
    DEBUG_2("cloning super particle #%u",spi);
    ui retval=add_sp(network.superparticles[spi].sptype);
    network.superparticles[retval].backdoor=network.superparticles[spi].backdoor;
    for(auto m: network.superparticles[spi].particles)
    {
        ui p=clone_particle(m.first,x);
        network.spid[p]=retval;
        network.superparticles[retval].particles[p]=m.second;
        network.superparticles[retval].backdoor[m.second]=p;
    }
    for(auto f:network.forcelibrary) for(auto u:f.particles) for(auto &v:u) if(network.spid[v]==spi) v=network.superparticles[retval].backdoor[network.superparticles[spi].particles[v]];
    return retval;
}

template<ui dim> void md<dim>::translate_sp(ui spi,ldf x[dim])
{
    DEBUG_2("translating super particle #%u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) translate_particle(it->first,x);
}

template<ui dim> void md<dim>::drift_sp(ui spi,ldf dx[dim])
{
    DEBUG_2("drifting Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) drift_particle(it->first,dx);
}

template<ui dim> void md<dim>::heat_sp(ui spi,ldf lambda)
{
    DEBUG_2("drifting Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) drift_particle(it->first,lambda);
}

template<ui dim> void md<dim>::set_position_sp(ui spi,ldf x[dim])
{
    DEBUG_2("drifting Translating super particle particle %u.",spi);
    avars.export_force_calc=true;
    ldf delx[dim];
    get_position_sp(spi,delx);
    for(ui d=0;d<dim;d++) delx[d]=x[d]-delx[d];
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) translate_particle(it->first,delx);
}

template<ui dim> void md<dim>::set_velocity_sp(ui spi,ldf dx[dim])
{
    DEBUG_2("drifting Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) memcpy(particles[it->first].dx,dx,dim*sizeof(ldf));
}

template<ui dim> void md<dim>::get_position_sp(ui spi,ldf x[dim])
{
    DEBUG_2("calculating center of mass super particle particle %u.",spi);
    ldf m=0.0;
    for(ui d=0;d<dim;d++) x[d]=0.0;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++)
    {
        for(ui d=0;d<dim;d++) x[d]+=particles[it->first].m*particles[it->first].x[d];
        m+=particles[it->first].m;
    }
    for(ui d=0;d<dim;d++) x[d]/=m;
}

template<ui dim> void md<dim>::get_velocity_sp(ui spi,ldf dx[dim])
{
    DEBUG_2("calculating velocity of center of mass super particle particle %u.",spi);
    ldf m=0.0;
    for(ui d=0;d<dim;d++) dx[d]=0.0;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++)
    {
        for(ui d=0;d<dim;d++) dx[d]+=particles[it->first].m*particles[it->first].dx[d];
        m+=particles[it->first].m;
    }
    for(ui d=0;d<dim;d++) dx[d]/=m;
}

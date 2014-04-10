#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_superparticle(ui sptype)
{
    ui spi = network.superparticles.size();
    network.superparticles.push_back(superparticle());
    network.superparticles[spi].sptype = sptype;
    DEBUG_2("added superparticle #%u with type %u",spi,sptype);
    return spi;
}

template<ui dim> bool md<dim>::rem_superparticle(ui spi)
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
        for(auto m: network.superparticles[spn].particles) network.spid[m.first]=numeric_limits<ui>::max();
        network.superparticles.pop_back();
        DEBUG_2("removed superparticle #%u",spi);
        return true;
    }
}

template<ui dim> ui md<dim>::sp_ingest(ui spi,ui i)
{
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle %d does not exist", spi);
        return numeric_limits<ui>::max();
    }
    else if(network.spid[i] < numeric_limits<ui>::max())
    {
        WARNING("particle %u is already in superparticle %u",i,network.spid[i]);
        return numeric_limits<ui>::max();
    }
    else
    {   
        network.spid[i]=spi;
        ui n=network.superparticles[spi].particles.size();
        DEBUG_2("particle #%u is ingested by superparticle #%u",i,spi);
        return network.superparticles[spi].particles[i]=n;
    }
}

template<ui dim> bool md<dim>::sp_dispose(ui i)
{
    ui spi=network.spid[i];
    if(spi==numeric_limits<ui>::max())
        return false;
    else if(network.superparticles[spi].particles.size()<2)
    {
        DEBUG_2("particle #%u is removed from superparticle #%u; removed superparticle",i,spi);
        return rem_superparticle(spi);
    }
    else
    {
        network.spid[i]=numeric_limits<ui>::max();
        DEBUG_2("particle #%u is removed from superparticle #%u",i,spi);
        return network.superparticles[spi].particles.erase(i);
    }
}

template<ui dim> ui md<dim>::sp_pid(ui spi,ui idx)
{
    for(auto m: network.superparticles[spi].particles) if(m.second==idx) return m.first;
    return numeric_limits<ui>::max();
}

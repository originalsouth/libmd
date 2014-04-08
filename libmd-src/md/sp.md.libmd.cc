#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::sp_ingest(ui spi,ui i)
{
    if(spi<network.superparticles.size())
    {
        if(network.spid[i]==spi)
        {
            WARNING("Praticle %u is already in super_particle %u.",i,spi);
            return spi;
        }
        network.spid[i]=spi;
        network.superparticles[spi].particles[i]=network.superparticles[spi].particles.size();
    }
    else
    {
        spi=network.superparticles.size();
        network.spid[i]=spi;
        superparticle sp;
        sp.particles[i]=sp.particles.size();
        sp.sptype=numeric_limits<ui>::max();
        network.superparticles.push_back(sp);
    }
    DEBUG_2("particle #%u is ingested by super partice #%u",i,spi);
    return spi;
}

template<ui dim> ui md<dim>::sp_ingest(ui spi,ui sptype,ui i)
{
    if(spi<network.superparticles.size())
    {
        network.spid[i]=spi;
        ui n=network.superparticles[spi].particles.size();
        network.superparticles[spi].particles[i]=n;
        network.superparticles[spi].sptype=sptype;
    }
    else
    {
        spi=network.superparticles.size();
        network.spid[i]=spi;
        superparticle sp;
        sp.particles[i]=0;
        sp.sptype=sptype;
        network.superparticles.push_back(sp);
    }
    DEBUG_2("particle #%u is ingested by super partice #%u with super particle type %u",i,spi,sptype);
    return spi;
}

template<ui dim> void md<dim>::sp_dispose(ui spi)
{
    if(spi<network.superparticles.size())
    {
        ui spn=network.superparticles.size()-1;
        if(spi<spn)
        {
            for(auto it=network.superparticles[spn].particles.begin();it!=network.superparticles[spn].particles.end();it++) network.spid[it->first]=spi;
            iter_swap(network.superparticles.begin()+spi,network.superparticles.end());
        }
        for(auto it=network.superparticles[spn].particles.begin();it!=network.superparticles[spn].particles.end();it++) network.spid[it->first]=numeric_limits<ui>::max();
        network.superparticles.pop_back();
    }
}

template<ui dim> void md<dim>::sp_p_dispose(ui i)
{
    if(i<N)
    {
        ui spi=network.spid[i];
        if(spi<network.superparticles.size())
        {
            if(network.superparticles[spi].particles.size()<2) sp_dispose(spi);
            else
            {
                network.spid[i]=numeric_limits<ui>::max();
                network.superparticles[spi].particles.erase(i);
            }
        }
    }
}

#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_particle(ldf mass,ui ptype,bool fixed)
{
    N++;
    particles.push_back(particle<dim>(mass,ptype,fixed));
    network.spid.push_back(numeric_limits<ui>::max());
    network.skins.resize(N);
    network.forces.resize(N);
    network.usedtypes[ptype].insert(N-1);
    index();
    return N-1;
}

template<ui dim> ui md<dim>::add_particle(ldf x[dim],ldf mass,ui ptype,bool fixed)
{
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("Created particle #%u",i);
    for(ui d=0;d<dim;d++)
    {
        particles[i].x[d]=x[d];
        particles[i].dx[d]=0.0;
    }
    return i;
}

template<ui dim> ui md<dim>::add_particle(ldf x[dim],ldf dx[dim],ldf mass,ui ptype,bool fixed)
{
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("Created particle #%u",i);
    for(ui d=0;d<dim;d++)
    {
        particles[i].x[d]=x[d];
        particles[i].dx[d]=dx[d];
    }
    return i;
}

template<ui dim> void md<dim>::rem_particle(ui i)
{
    if(network.spid[i]<N) sp_dispose(network.spid[i],i);
    N--;
    // store particle types of deleted particle and particle that will replace it
    ui deleted_ptype = particles[i].type;
    ui last_ptype = particles.rbegin()->type;
    // swap particle to delete with last particle, to prevent changing index of all particles after particlenr, and then delete it
    std::iter_swap(particles.begin()+i, particles.rbegin());
    particles.pop_back();
    // update usedtypes dictionary
    network.usedtypes[deleted_ptype].erase(i);
    network.usedtypes[last_ptype].erase(N);
    network.usedtypes[last_ptype].insert(i);
    // update the network
    std::iter_swap(network.skins.begin()+i, network.skins.rbegin());
    network.skins.pop_back();
    std::iter_swap(network.forces.begin()+i, network.forces.rbegin());
    network.forces.pop_back();
    index();
}

template<ui dim> void md<dim>::fix_particle(ui i,bool fix)
{
    DEBUG_2("Fixing(%d) particle %u.",fix,i);
    particles[i].fix=fix;
}

template<ui dim> void md<dim>::fix_particles(ui spi,bool fix)
{
    DEBUG_2("Fixing(%d) super particle particle %u.",fix,spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) particles[it->first].fix=fix;
}

template<ui dim> ui md<dim>::clone_particle(ui i,ldf x[dim])
{
    ui retval=add_particle(particles[i].x,particles[i].dx,particles[i].mass,particles[i].type,particles[i].fix);
    memcpy(particles[retval].xp,particles[i].xp,dim*sizeof(ldf));
    translate_particle(retval,x);
    return retval;
}

template<ui dim> ui md<dim>::clone_particles(ui spi,ldf x[dim])
{
    ui retval=network.superparticles.size(),dummy;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++)
    {
        dummy=clone_particle(it->first,x);
        sp_ingest(retval,network.superparticles[spi].sptype,dummy);
    }
    return retval;
}

template<ui dim> void md<dim>::translate_particle(ui i,ldf x[dim])
{
    DEBUG_2("Translating particle particle %u.",i);
    for(ui d=0;d<dim;d++) particles[i].x[d]+=x[d];
}

template<ui dim> void md<dim>::translate_particles(ui spi,ldf x[dim])
{
    DEBUG_2("Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].x[d]+=x[d];
}

template<ui dim> void md<dim>::drift_particle(ui i,ldf dx[dim])
{
    DEBUG_2("Drifting particle particle %u.",i);
    for(ui d=0;d<dim;d++) particles[i].dx[d]+=dx[d];
}

template<ui dim> void md<dim>::drift_particles(ui spi,ldf dx[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].dx[d]+=dx[d];
}

template<ui dim> void md<dim>::set_position_particles(ui spi,ldf x[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    ldf dx[dim];
    get_position_particles(spi,dx);
    for(ui d=0;d<dim;d++) dx[d]=x[d]-dx[d];
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].x[d]+=dx[d];
}

template<ui dim> void md<dim>::set_velocity_particles(ui spi,ldf dx[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].dx[d]=dx[d];
}

template<ui dim> void md<dim>::get_position_particles(ui spi,ldf x[dim])
{
    DEBUG_2("Calculating center of mass super particle particle %u.",spi);
    ldf m=0.0;
    for(ui d=0;d<dim;d++) x[d]=0.0;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++)
    {
        x[d]+=particles[it->first].m*particles[it->first].x[d];
        m+=particles[it->first].m;
    }
    for(ui d=0;d<dim;d++) x[d]/=m;
}

template<ui dim> void md<dim>::get_velocity_particles(ui spi,ldf dx[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    for(ui d=0;d<dim;d++) dx[d]=0.0;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) dx[d]+=particles[it->first].dx[d];
    for(ui d=0;d<dim;d++) dx[d]/=network.superparticles[spi].particles.size();
}

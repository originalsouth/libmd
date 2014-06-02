#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_particle(ldf mass,ui ptype,bool fixed)
{
    avars.export_force_calc=true;
    N++;
    particles.push_back(particle<dim>(mass,ptype,fixed));
    network.spid.push_back(UI_MAX);
    network.skins.resize(N);
    network.forces.resize(N);
    for(auto lib: network.forcelibrary) lib.particles.resize(N);
    avars.reindex=true;
    return N-1;
}

template<ui dim> ui md<dim>::add_particle(ldf x[dim],ldf mass,ui ptype,bool fixed)
{
    ldf tempx[dim];
    memcpy(tempx,x,dim*sizeof(ldf));
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("created particle #%u with given position",i);
    memcpy(particles[i].x,tempx,dim*sizeof(ldf));
    memset(particles[i].dx,0,dim*sizeof(ldf));
    return i;
}

template<ui dim> ui md<dim>::add_particle(ldf x[dim],ldf dx[dim],ldf mass,ui ptype,bool fixed)
{
    ldf tempx[dim],tempdx[dim];
    memcpy(tempx,x,dim*sizeof(ldf));
    memcpy(tempdx,dx,dim*sizeof(ldf));
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("created particle #%u with given position and velocity",i);
    memcpy(particles[i].x,tempx,dim*sizeof(ldf));
    memcpy(particles[i].dx,tempdx,dim*sizeof(ldf));
    return i;
}

template<ui dim> void md<dim>::rem_particle(ui i)
{
    DEBUG_2("removing particle %u.",i);
    if(network.spid[i]<UI_MAX)
        network.superparticles[network.spid[i]].particles.erase(i);
    ui j, k, p;
    if (i < N-1)
    {   // swap particle to delete with last particle, to prevent changing index of all particles after particlenr, and then delete it
        std::iter_swap(particles.begin()+i, particles.rbegin());
        // update the network
        std::iter_swap(network.skins.begin()+i, network.skins.rbegin());
        std::iter_swap(network.forces.begin()+i, network.forces.rbegin());
        // Modify skins
        for (j = network.skins[i].size()-1; j < UI_MAX; j--)
        {   p = network.skins[i][j].neighbor;
            if (p == i)
                p = N-1;
            for (k = network.skins[p].size()-1; k < UI_MAX && network.skins[p][k].neighbor != N-1; k--);
            if (k > N)
            {   ERROR("(formerly) last particle not found in skinlist of particle %d", p);
                return;
            }
            network.skins[p][k].neighbor = i;
        }
        // Modify superparticle
        if ((network.spid[i] = network.spid[N-1]) < UI_MAX)
        {   auto it = network.superparticles[network.spid[i]].particles.find(N-1);
            if (it == network.superparticles[network.spid[i]].particles.end())
            {   ERROR("particle #%u not found in superparticle #%u", N-1, network.spid[i]);
                return;
            }
            network.superparticles[network.spid[i]].particles[i] = it->second;
            network.superparticles[network.spid[i]].particles.erase(it);
        }
    }
    // Modify skins
    for (j = network.skins[N-1].size()-1; j < UI_MAX; j--)
    {   p = network.skins[N-1][j].neighbor;
        for (k = network.skins[p].size()-1; k < UI_MAX && network.skins[p][k].neighbor != i; k--);
        if (k > N)
        {   ERROR("Particle to be deleted not found in skinlist of particle %d", p);
            return;
        }
        std::iter_swap(network.skins[p].begin()+k, network.skins[p].rbegin());
        network.skins[p].pop_back();
    }
    // Remove last particle
    N--;
    particles.pop_back();
    network.skins.pop_back();
    network.forces.pop_back();
    network.spid.pop_back();
    avars.reindex=true;
}

template<ui dim> void md<dim>::fix_particle(ui i,bool fix)
{
    DEBUG_2("fixing(%d) particle %u.",fix,i);
    particles[i].fix=fix;
}

template<ui dim> ui md<dim>::clone_particle(ui i,ldf x[dim])
{
    DEBUG_2("cloning particle #%u",i);
    ui retval=add_particle();
    particles[retval]=particles[i];
    translate_particle(retval,x);
    network.forces[retval]=network.forces[i];
    for(auto j:network.forces[i]) network.forcelibrary[j].particles[retval]=network.forcelibrary[j].particles[i];
    for(auto f:network.forcelibrary) for(auto u:f.particles) for(auto v:u) if(v==i) u.push_back(retval);
    return retval;
}

template<ui dim> void md<dim>::translate_particle(ui i,ldf x[dim])
{
    DEBUG_2("translating particle #%u.",i);
    avars.export_force_calc=true;
    for(ui d=0;d<dim;d++)
    {
        particles[i].x[d]+=x[d];
        particles[i].xp[d]+=x[d];
    }
    thread_periodicity(i);
    avars.reindex=true;
}

template<ui dim> void md<dim>::drift_particle(ui i,ldf dx[dim])
{
    DEBUG_2("drifting particle particle %u.",i);
    for(ui d=0;d<dim;d++) particles[i].dx[d]+=dx[d];
}

template<ui dim> void md<dim>::heat_particle(ui i,ldf lambda)
{
    DEBUG_2("drifting particle particle %u.",i);
    for(ui d=0;d<dim;d++) particles[i].dx[d]*=lambda;
}

template<ui dim> void md<dim>::uitopptr(vector<particle<dim>*> *x,vector<ui> i)
{
    ui Ni=i.size();
    for(ui j=0;j<Ni and j<N;j++) x->push_back(&particles[i[j]]);
};

template<ui dim> ui md<dim>::pptrtoui(particle<dim> *x)
{
    return x-&particles[0];
};

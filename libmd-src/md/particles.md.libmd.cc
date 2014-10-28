#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_particle(ldf mass,ui ptype,bool fixed)
{
    //!
    //! This function adds a particle to the system and returns its index.
    //! Optionally, provide its mass (default: 1.0), type (default: 0) and/or whether it is fixed (default: false).
    //! The position and velocity of the new particle are not set.
    //!
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
    //!
    //! This function adds a particle to the system at position <tt>x[]</tt> and returns its index.
    //! Optionally, provide its mass (default: 1.0), type (default: 0) and/or whether it is fixed (default: false).
    //! The velocity of the new particle is not set.
    //!
    ldf tempx[dim];
    memcpy(tempx,x,dim*sizeof(ldf));
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("created particle #" F_UI " with given position",i);
    memcpy(particles[i].x,tempx,dim*sizeof(ldf));
    memset(particles[i].dx,0,dim*sizeof(ldf));
    return i;
}

template<ui dim> ui md<dim>::add_particle(ldf x[dim],ldf dx[dim],ldf mass,ui ptype,bool fixed)
{
    //!
    //! This function adds a particle to the system at position <tt>x[]</tt> with velocity <tt>dx[]</tt> and returns its index.
    //! Optionally, provide its mass (default: 1.0), type (default: 0) and/or whether it is fixed (default: false).
    //!
    ldf tempx[dim],tempdx[dim];
    memcpy(tempx,x,dim*sizeof(ldf));
    memcpy(tempdx,dx,dim*sizeof(ldf));
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("created particle #" F_UI " with given position and velocity",i);
    memcpy(particles[i].x,tempx,dim*sizeof(ldf));
    memcpy(particles[i].dx,tempdx,dim*sizeof(ldf));
    return i;
}

template<ui dim> void md<dim>::rem_particle(ui i)
{
    //!
    //! This function removes particle <tt>i</tt> and all references to it from the system.<br>
    //! Note: particle <tt>N-1</tt> becomes the new particle <tt>i</tt> (if \f$i<N-1\f$).
    //!
    DEBUG_2("removing particle #" F_UI "",i);
    if(network.spid[i]<UI_MAX)
        network.superparticles[network.spid[i]].particles.erase(i);
    ui j, k, p;
    if (i < N-1)
    {   // swap particle to delete with last particle, to prevent changing index of all particles after particlenr, and then delete it
        iter_swap(particles.begin()+i, particles.rbegin());
        // update the network
        iter_swap(network.skins.begin()+i, network.skins.rbegin());
        iter_swap(network.forces.begin()+i, network.forces.rbegin());
        // Modify skins
        for (j = network.skins[i].size()-1; j < UI_MAX; j--)
        {   p = network.skins[i][j].neighbor;
            if (p == i)
                p = N-1;
            for (k = network.skins[p].size()-1; k < UI_MAX && network.skins[p][k].neighbor != N-1; k--);
            if (k > N)
            {   ERROR("(formerly) last particle not found in skinlist of particle #%d", p);
                return;
            }
            network.skins[p][k].neighbor = i;
        }
        // Modify superparticle
        if ((network.spid[i] = network.spid[N-1]) < UI_MAX)
        {   auto it = network.superparticles[network.spid[i]].particles.find(N-1);
            if (it == network.superparticles[network.spid[i]].particles.end())
            {   ERROR("particle #" F_UI " not found in superparticle #" F_UI "", N-1, network.spid[i]);
                return;
            }
            network.superparticles[network.spid[i]].particles[i] = it->second;
            network.superparticles[network.spid[i]].particles.erase(it);
        }
        for(auto ftype:network.forcelibrary) if(!ftype.particles.empty()) iter_swap(ftype.particles.begin()+i,ftype.particles.rbegin());
    }
    for(auto ftype:network.forcelibrary) if(!ftype.particles.empty())
    {
        ftype.particles.pop_back();
        for(j=0;j<N-1;j++) for(k=ftype.particles[j].size()-1;k<UI_MAX;k--) if(ftype.particles[j][k]==i)
        {
            iter_swap(ftype.particles[j].begin()+k,ftype.particles[j].rbegin());
            ftype.particles[j].pop_back();
        }
    }
    if(i<N-1) for(auto ftype:network.forcelibrary) if(!ftype.particles.empty()) for(j=0;j<N-1;j++) for(k=ftype.particles[j].size()-1;k<UI_MAX;k--) if(ftype.particles[j][k]==N-1) ftype.particles[j][k]=i;
    // Modify skins
    for (j = network.skins[N-1].size()-1; j < UI_MAX; j--)
    {   p = network.skins[N-1][j].neighbor;
        for (k = network.skins[p].size()-1; k < UI_MAX && network.skins[p][k].neighbor != i; k--);
        if (k > N)
        {   ERROR("particle to be deleted not found in skinlist of particle #%d", p);
            return;
        }
        iter_swap(network.skins[p].begin()+k, network.skins[p].rbegin());
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
    //!
    //! This function fixes (<tt>fix=true</tt>) or unfixes (<tt>fix=false</tt>) particle <tt>i</tt>
    //!
    DEBUG_2("%sfixing particle #" F_UI "",fix?"":"un",i);
    particles[i].fix=fix;
}

template<ui dim> ui md<dim>::clone_particle(ui i,ldf x[dim])
{
    //!
    //! This function creates a new particle that is a copy of particle <tt>i</tt>.
    //! The position of the new particle is translated by the vector <tt>x[]</tt> with respect to the original one.
    //! It returns the index of the new particle.
    //!
    DEBUG_2("cloning particle #" F_UI "",i);
    ui retval=add_particle();
    particles[retval]=particles[i];
    translate_particle(retval,x);
    network.forces[retval]=network.forces[i];
    for(auto j:network.forces[i]) if(!network.forcelibrary[j].particles.empty()) network.forcelibrary[j].particles[retval]=network.forcelibrary[j].particles[i];
    for(auto f:network.forcelibrary) for(auto u:f.particles) for(auto v:u) if(v==i) u.push_back(retval);
    return retval;
}

template<ui dim> void md<dim>::translate_particle(ui i,ldf x[dim])
{
    //!
    //! This function translates particle <tt>i</tt> by the vector <tt>x[]</tt>.
    //!
    DEBUG_2("translating particle #" F_UI "",i);
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
    //!
    //! This function adds the vector <tt>dx[]</tt> to the velocity of particle <tt>i</tt>.
    //!
    DEBUG_2("drifting particle #" F_UI "",i);
    for(ui d=0;d<dim;d++) particles[i].dx[d]+=dx[d];
}

template<ui dim> void md<dim>::heat_particle(ui i,ldf lambda)
{
    //!
    //! This function increases the velocity of particle <tt>i</tt> by a factor of <tt>lambda</tt>.
    //!
    DEBUG_2("heating particle #" F_UI "",i);
    for(ui d=0;d<dim;d++) particles[i].dx[d]*=lambda;
}

#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_sp(ui sptype)
{
    //!
    //! This function adds a superparticle of type <tt>sptype</tt> to <tt>network.superparticles[]</tt> and returns its index
    //!
    ui spi = network.superparticles.size();
    network.superparticles.push_back(superparticle());
    network.superparticles[spi].sptype = sptype;
    DEBUG_2("added superparticle #" F_UI " with type " F_UI "",spi,sptype);
    return spi;
}

template<ui dim> bool md<dim>::rem_sp(ui spi)
{
    //!
    //! This function removes the superparticle with index <tt>spi</tt> from <tt>network.superparticles[]</tt>
    //! and resets the superparticle reference (<tt>network.spid[]</tt>) of all its particles.
    //! It returns whether the given superparticle existed.<br>
    //! Note: the last superparticle in <tt>network.superparticles[]</tt> takes the place of the removed one
    //!
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle #%d does not exist", spi);
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
        DEBUG_2("removed superparticle #" F_UI "",spi);
        return true;
    }
}

template<ui dim> bool md<dim>::rem_sp_particles(ui spi)
{
    //!
    //! This function removes all the particles of the superparticle with index <tt>spi</tt> from the system
    //! and clears the superparticle (but does not remove it).
    //! It returns whether the given superparticle exists.
    //!
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle #%d does not exist", spi);
        return false;
    }
    else
    {
        ui k=0;
        vector<ui> target(network.superparticles[spi].particles.size());
        for(auto m: network.superparticles[spi].particles) target[k++]=m.first;
        for(auto i: target) rem_particle(i);
        network.superparticles[spi].particles.clear();
        network.superparticles[spi].backdoor.clear();
        DEBUG_2("removed particles of superparticle #" F_UI "",spi);
        return true;
    }
}

template<ui dim> ui md<dim>::sp_ingest(ui spi,ui i,ui idx)
{
    //!
    //! This function puts particle <tt>i</tt> in slot <tt>idx</tt> of the superparticle with index <tt>spi</tt>.
    //! If <tt>idx</tt> = <tt>UI_MAX</tt>, the index is set to highest index plus one.
    //! It returns the index.
    //! If the given superparticle does not exist, or the given particle is already in a superparticle, or the given slot is already taken,
    //! nothing is done and UI_MAX is returned.
    //!
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle #%d does not exist", spi);
        return UI_MAX;
    }
    else if(network.spid[i] < UI_MAX)
    {
        WARNING("particle #" F_UI " is already in superparticle #" F_UI "",i,network.spid[i]);
        return UI_MAX;
    }
    else
    {
        if(idx>=network.superparticles[spi].backdoor.size())
        {
            if(idx==UI_MAX) idx=network.superparticles[spi].backdoor.size();
            network.superparticles[spi].backdoor.resize(idx+1,UI_MAX);
        }
        else if(network.superparticles[spi].backdoor[idx]<UI_MAX)
        {
            WARNING("superparticle #" F_UI " already contains index #" F_UI "",spi,idx);
            return UI_MAX;
        }
        network.spid[i]=spi;
        network.superparticles[spi].backdoor[idx]=i;
        DEBUG_2("particle #" F_UI " is ingested by superparticle #" F_UI "",i,spi);
        avars.reindex=true;
        return network.superparticles[spi].particles[i]=idx;
    }
}

template<ui dim> bool md<dim>::sp_dispose(ui i)
{
    //!
    //! This function removes particle <tt>i</tt> from its superparticle (but not from the system)
    //! It returns whether the given particle was in a superparticle.
    //!
    ui spi=network.spid[i];
    if(spi==UI_MAX) return false;
    else
    {
        network.spid[i]=UI_MAX;
        DEBUG_2("particle #" F_UI " is removed from superparticle #" F_UI "",i,spi);
        ui j=network.superparticles[spi].particles[i];
        if(j==network.superparticles[spi].backdoor.size()-1)
        {
            network.superparticles[spi].backdoor.pop_back();
            while(!network.superparticles[spi].backdoor.empty() and network.superparticles[spi].backdoor.back()==UI_MAX) network.superparticles[spi].backdoor.pop_back();
        }
        else network.superparticles[spi].backdoor[j]=UI_MAX;
        network.superparticles[spi].particles.erase(i);
        avars.reindex=true;
        return true;
    }
}

template<ui dim> bool md<dim>::sp_dispose_idx(ui spi,ui idx)
{
    //!
    //! This function removes the particle from slot <tt>idx</tt> of the superparticle with index <tt>spi</tt> from the superparticle
    //! (but not from the system).
    //! It returns whether the given superparticle exists and has a particle in the given slot.
    //!
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle #" F_UI " does not exist",spi);
        return false;
    }
    else if(idx>=network.superparticles[spi].backdoor.size())
    {
        WARNING("superparticle #" F_UI " does not contain index #" F_UI "",spi,idx);
        return false;
    }
    else
    {
        ui i=network.superparticles[spi].backdoor[idx];
        if(i<UI_MAX)
        {
            sp_dispose(i);
            DEBUG_2("particle #" F_UI " is removed from superparticle #" F_UI "",i,spi);
            return true;
        }
        else
        {
            WARNING("superparticle #" F_UI " does not contain index #" F_UI "",spi,idx);
            return false;
        }
    }
}

template<ui dim> ui md<dim>::sp_pid(ui spi,ui idx)
{
    //!
    //! This function returns the id of the particle in slot <tt>idx</tt> of the superparticle with index <tt>spi</tt>.
    //! If the superparticle does not exist or it does not have a particle in the given slot, it returns UI_MAX.
    //!
    if(spi>=network.superparticles.size())
    {
        WARNING("superparticle #" F_UI " does not exist",spi);
        return UI_MAX;
    }
    else if(idx>=network.superparticles[spi].backdoor.size())
    {
        WARNING("superparticle #" F_UI " does not contain index #" F_UI "",spi,idx);
        return UI_MAX;
    }
    else return network.superparticles[spi].backdoor[idx];
}

template<ui dim> void md<dim>::fix_sp(ui spi,bool fix)
{
    //!
    //! This function fixes (<tt>fix=true</tt>) or unfixes (<tt>fix=false</tt>) all the particles
    //! belonging to the superparticle with index <tt>spi</tt>.
    //!
    DEBUG_2("%sfixing superparticle #" F_UI "",fix?"":"un",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) particles[it->first].fix=fix;
}

template<ui dim> ui md<dim>::clone_sp(ui spi,ldf x[dim])
{
    //!
    //! This function creates a new superparticle that is a copy of the superparticle with index <tt>spi</tt>, with new particles.
    //! The position of all the new particles are translated by the vector <tt>x[]</tt> with respect to the original ones.
    //! It returns the index of the new superparticle in <tt>network.superparticles[]</tt>.
    //!
    DEBUG_2("cloning superparticle #" F_UI "",spi);
    ui retval=add_sp(network.superparticles[spi].sptype);
    network.superparticles[retval].backdoor=network.superparticles[spi].backdoor;
    ui K=network.superparticles[spi].backdoor.size(),i;
    for(ui k=0;k<K;k++)
    {
        if((i=network.superparticles[spi].backdoor[k])<UI_MAX)
        {
            ui p=clone_particle(i,x);
            network.spid[p]=retval;
            network.superparticles[retval].particles[p]=k;
            network.superparticles[retval].backdoor[k]=p;
        }
    }
    for(auto f:network.forcelibrary) for(auto u:f.particles) for(auto &v:u) if(network.spid[v]==spi) v=network.superparticles[retval].backdoor[network.superparticles[spi].particles[v]];
    return retval;
}

template<ui dim> void md<dim>::translate_sp(ui spi,ldf x[dim])
{
    //!
    //! This function translates all particles belonging to the superparticle with index <tt>spi</tt> by the vector <tt>x[]</tt>.
    //!
    DEBUG_2("translating superparticle #" F_UI "",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) translate_particle(it->first,x);
}

template<ui dim> void md<dim>::drift_sp(ui spi,ldf dx[dim])
{
    //!
    //! This function adds the vector <tt>dx[]</tt> to the velocities all particles belonging to the superparticle with index <tt>spi</tt>.
    //!
    DEBUG_2("drifting superparticle #" F_UI "",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) drift_particle(it->first,dx);
}

template<ui dim> void md<dim>::heat_sp(ui spi,ldf lambda)
{
    //!
    //! This function increases the velocity of all particles belonging to the superparticle with index <tt>spi</tt> by a factor of <tt>lambda</tt>.
    //!
    DEBUG_2("heating superparticle #" F_UI "",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) drift_particle(it->first,lambda);
}

template<ui dim> void md<dim>::set_position_sp(ui spi,ldf x[dim])
{
    //!
    //! This function translates the superparticle with index <tt>spi</tt> such that its center of mass is at position <tt>x[]</tt>.
    //!
    DEBUG_2("repositioning superparticle #" F_UI "",spi);
    avars.export_force_calc=true;
    ldf delx[dim];
    get_position_sp(spi,delx);
    for(ui d=0;d<dim;d++) delx[d]=x[d]-delx[d];
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) translate_particle(it->first,delx);
}

template<ui dim> void md<dim>::set_velocity_sp(ui spi,ldf dx[dim])
{
    //!
    //! This function sets the velocity of all particles belonging to the superparticle with index <tt>spi</tt> equal to <tt>dx[]</tt>.
    //!
    DEBUG_2("setting velocity of superparticle #" F_UI "",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) memcpy(particles[it->first].dx,dx,dim*sizeof(ldf));
}

template<ui dim> void md<dim>::get_position_sp(ui spi,ldf x[dim])
{
    //!
    //! This function puts the position of the center of mass of the superparticle with index <tt>spi</tt> in the vector <tt>x[]</tt>.
    //!
    DEBUG_2("calculating center of mass of superparticle #" F_UI "",spi);
    ldf m=0.0;
    memset(x,0,dim*sizeof(ldf));
    ui i=network.superparticles[spi].particles.begin()->first;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++)
    {
        for(ui d=0;d<dim;d++) x[d]+=particles[it->first].m*dd(d,i,it->first);
        m+=particles[it->first].m;
    }
    for(ui d=0;d<dim;d++) x[d]=particles[i].x[d]+x[d]/m;
    thread_periodicity(x);
}

template<ui dim> void md<dim>::get_velocity_sp(ui spi,ldf dx[dim])
{
    //!
    //! This function puts the velocity of the center of mass of the superparticle with index <tt>spi</tt> in the vector <tt>dx[]</tt>.
    //!
    DEBUG_2("calculating velocity of center of mass of superparticle #" F_UI "",spi);
    ldf m=0.0;
    memset(dx,0,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) dx[d]=0.0;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++)
    {
        for(ui d=0;d<dim;d++) dx[d]+=particles[it->first].m*particles[it->first].dx[d];
        m+=particles[it->first].m;
    }
    for(ui d=0;d<dim;d++) dx[d]/=m;
}

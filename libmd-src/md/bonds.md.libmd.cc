#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::update_skins(ui p1, ui p2)
{
    //!
    //! Updates the skin lists of <tt>p1</tt> and <tt>p2</tt>. This function
    //! is called by functions that modify pairwise interactions between
    //! specific particle pairs, to rebuild the skins taking into account
    //! new particle ID assignments and interaction definitions.
    //!
    ui interaction, s, t, j, k;
    pair<ui,ui> id;
    if ((s = network.spid[p1]) < UI_MAX && s == network.spid[p2]
        && network.sptypes[t = network.superparticles[s].sptype].splookup.count(id = network.hash(network.superparticles[s].particles[p1], network.superparticles[s].particles[p2])))
        interaction = network.sptypes[t].splookup[id];
    else if (network.lookup.count(id = network.hash(particles[p1].type, particles[p2].type)))
        interaction = network.lookup[id];
    else
        interaction = UI_MAX;
    for (j = network.skins[p1].size()-1; j < UI_MAX && network.skins[p1][j].neighbor != p2; j--);
    if (j < UI_MAX)
    {   for (k = network.skins[p2].size()-1; network.skins[p2][k].neighbor != p1; k--);
        if (interaction < UI_MAX)
            network.skins[p1][j].interaction = network.skins[p2][k].interaction = interaction;
        else
        {   std::iter_swap(network.skins[p1].begin()+j, network.skins[p1].rbegin());
            network.skins[p1].pop_back();
            std::iter_swap(network.skins[p2].begin()+k, network.skins[p2].rbegin());
            network.skins[p2].pop_back();
        }
    }
    else if (interaction < UI_MAX)
    {   network.skins[p1].push_back(interactionneighbor(p2,interaction));
        network.skins[p2].push_back(interactionneighbor(p1,interaction));
    }
}

template<ui dim> bool md<dim>::add_bond(ui p1, ui p2, ui interaction)
{   
    //!
    //! This function adds a bond between particles <tt>p1</tt> and <tt>p2</tt> of 
    //! previously-defined interaction type
    //! <tt>interaction</tt>.
    //! (See the [Interactions](@ref md-pairpotentials) documentation for more 
    //! information on interaction types.)
    //! <br><br>
    //! Does nothing and returns false if one of the following holds:
    //! <ul>
    //!    <li>an interaction already exists between <tt>p1</tt> and <tt>p2</tt>,
    //!    <li>the interaction specified by the third argument does not exist or was previously removed.
    //! </ul>
    //! Otherwise, creates a bond of the specified type between <tt>p1</tt> and <tt>p2</tt>, and returns <tt>true</tt>.
    //! 
    //!
    if (network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
        return false;
    else if(interaction>=network.library.size())
    {
        WARNING("interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("interaction %d was previously removed", interaction);
        return false;
    }
    else
    {   assign_unique_types(p1, p2);
        network.lookup[network.hash(particles[p1].type,particles[p2].type)]=interaction;
        update_skins(p1,p2);
        return true;
    }
}

template<ui dim> bool md<dim>::mod_bond(ui p1, ui p2, ui interaction)
{   
    //!
    //! This function modifies any previous interaction between particles <tt>p1</tt> and <tt>p2</tt>,
    //! and turns it into a bond of previously-defined interaction type specified in the third argument.
    //! (See the [Interactions](@ref md-pairpotentials) documentation for more 
    //! information on interaction types.)
    //! <br><br>
    //! Does nothing and returns false if one of the following holds:
    //! <ul>
    //!    <li>no interaction exists between <tt>p1</tt> and <tt>p2</tt>,
    //!    <li>the interaction specified by the third argument does not exist or was previously removed.
    //! </ul>
    //! Otherwise, creates a bond of the specified type between <tt>p1</tt> and <tt>p2</tt>, and returns <tt>true</tt>.
    //! 
    //!
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
        return false;
    else if(interaction>=network.library.size())
    {
        WARNING("interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("interaction %d was previously removed", interaction);
        return false;
    }
    else
    {   assign_unique_types(p1, p2);
        network.lookup[network.hash(particles[p1].type,particles[p2].type)]=interaction;
        update_skins(p1,p2);
        return true;
    }
}

template<ui dim> void md<dim>::mad_bond(ui p1, ui p2, ui interaction)
{   
    //!
    //! Same as md<dim>::add_bond(ui p1, ui p2, ui interaction), but performs no checks.
    //! <br><br>
    //! <b>Warning:</b> Assumes that the specified interaction
    //! type actually exists, and replaces any previously defined interaction
    //! between <tt>p1</tt> and <tt>p2</tt>.
    //! 
    //!
    assign_unique_types(p1, p2);
    network.lookup[network.hash(particles[p1].type,particles[p2].type)]=interaction;
    update_skins(p1,p2);
}

template<ui dim> bool md<dim>::add_bond(ui p1, ui p2, ui potential, vector<ldf> &parameters)
{
    //!
    //! This function adds a bond between particles <tt>p1</tt> and <tt>p2</tt> with
    //! the potential type and parameters specified in the arguments.
    //! (See the [Interactions](@ref md-pairpotentials) documentation for more 
    //! information on potential types and passing parameters.)
    //! <br><br>
    //! Does nothing if an interaction already exists between <tt>p1</tt> and <tt>p2</tt>, and returns <tt>false</tt>.
    //! Otherwise, creates a bond of the specified type between <tt>p1</tt> and <tt>p2</tt>, and returns <tt>true</tt>.
    //! 
    //!
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
    {   mad_bond(p1, p2, add_interaction(potential, parameters));
        return true;
    }
    else
        return false;
}

template<ui dim> bool md<dim>::mod_bond(ui p1, ui p2, ui potential, vector<ldf> &parameters)
{   
    //!
    //! This function modifies any previous interaction between particles <tt>p1</tt> and <tt>p2</tt> 
    //! into a pair-specific bond with the potential type and parameters
    //! specified in the arguments.
    //! (See the [Interactions](@ref md-pairpotentials) documentation for more 
    //! information on potential types and passing parameters.)
    //! <br><br>
    //! Does nothing if an interaction already exists between <tt>p1</tt> and <tt>p2</tt>, and returns <tt>false</tt>.
    //! Otherwise, creates a bond of the specified type between <tt>p1</tt> and <tt>p2</tt>, and returns <tt>true</tt>.
    //! 
    //!
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
        return false;
    else
    {   mad_bond(p1, p2, add_interaction(potential, parameters));
        return true;
    }
}

template<ui dim> void md<dim>::mad_bond(ui p1, ui p2, ui potential, vector<ldf> &parameters)
{   
    //!
    //! Same as md<dim>::add_bond(ui p1, ui p2, ui potential, vector<ldf> &parameters), but performs no checks.
    //! <br><br>
    //! <b>Warning:</b> Replaces any previously defined interaction
    //! between <tt>p1</tt> and <tt>p2</tt>.
    //! 
    //!
    mad_bond(p1, p2, add_interaction(potential, parameters));
}

template<ui dim> bool md<dim>::rem_bond(ui p1, ui p2)
{   
    //!
    //! If a prior interaction between particles <tt>p1</tt> and <tt>p2</tt>
    //! exists, removes it and returns <tt>true</tt>, else returns false.
    //! 
    //!
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
        return false;
    else
    {   assign_unique_types(p1,p2);
        network.lookup.erase(network.hash(particles[p1].type,particles[p2].type));
        update_skins(p1,p2);
        return true;
    }
}

template<ui dim> void md<dim>::assign_unique_types(ui p1, ui p2)
{
    //!
    //! Assign unique particle types to <tt>p1</tt> and <tt>p2</tt>,
    //! that are not shared by any other particles. Preserves all pairwise
    //! interactions that existed between <tt>p1</tt> or <tt>p2</tt> and
    //! any other particles.
    //! 
    //!

    // assign unique types to points
    ui P[2], old_type[2], new_type[2];
    P[0] = p1;
    P[1] = p2;
    ui maxtype = network.lookup.empty() ? 0 : network.lookup.rbegin()->first.first, i, j, q;

    for (j = 0; j < 2; j++)
        new_type[j] = old_type[j] = particles[P[j]].type;
    for (i = 0; i < N; i++)
        if (maxtype < particles[i].type)
            maxtype = particles[i].type;

    for (j = 0; j < 2; j++)
    {   for (i = 0; i < N && (i == P[j] || particles[i].type != old_type[j]); i++); // Check if there is at least one other particle with same type as p1
        if (i < N) // P[j] does not have a unique type; reassign.
            set_type(P[j], new_type[j] = ++maxtype); // use one greater than largest used particle type
        if (old_type[0] == old_type[1])
        {   set_type(P[1], new_type[1] = new_type[0]);
            break;
        }
    }

    for (j = 0; j < 2; j++)
        if (new_type[j] != old_type[j])
        {   // clone previously defined interactions
            map<pair<ui,ui>,ui>::iterator it = (old_type[j]==0 ? network.lookup.begin() : network.lookup.upper_bound(pair<ui,ui>(old_type[j]-1, UI_MAX)));
            for (; it != network.lookup.end() && it->first.first == old_type[j]; it++)
                if ((q = it->first.second) != new_type[1-j])
                    mad_typeinteraction(new_type[j], q, it->second);
            for (; it != network.lookup.end(); it++)
                if (it->first.second == old_type[j] && (q = it->first.first) != new_type[1-j])
                    mad_typeinteraction(new_type[j], q, it->second);
            if (old_type[0] == old_type[1])
                break;
        }
}

template<ui dim> void md<dim>::add_spring(ui p1, ui p2, ldf springconstant, ldf l0)
{   
    //! 
    //! Create a hookean spring [potential type HOOKEAN()] between <tt>p1</tt> and <tt>p2</tt>, with spring constant 
    //! and rest length prescribed by the third and fourth arguments respectively.
    //! 
    /* add a spring between two points with specified springconstant and equilibrium length */
    vector<ldf> params = {springconstant, l0};
    add_bond(p1,p2,POT::HOOKEAN,params);
}

// TODO: document superparticle bond functions

template<ui dim> bool md<dim>::add_sp_bond(ui p1, ui p2, ui interaction)
{   
    //!
    //! This function creates a new superparticle type and assigns it to the superparticle that <tt>p1</tt> and <tt>p2</tt> belong to.
    //! It then adds an interaction between the two particles, using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! If the two particles are not in the same superparticle, or have an interaction already, it does nothing and returns false.
    //!
    ui spi = network.spid[p1];
    if (spi == UI_MAX || spi != network.spid[p2])
        return false;
    ui spt = network.superparticles[spi].sptype;
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    if (network.sptypes[spt].splookup.count(id))
        return false;
    network.sptypes[clone_sptype(spi)].splookup[id] = interaction;
    update_skins(p1,p2);
    return true;
}

template<ui dim> bool md<dim>::mod_sp_bond(ui p1, ui p2, ui interaction)
{   
    //!
    //! This function creates a new superparticle type and assigns it to the superparticle that <tt>p1</tt> and <tt>p2</tt> belong to.
    //! It then modifies the interaction between the two particles, using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! If the two particles are not in the same superparticle, or do not have an interaction already, it does nothing and returns false.
    //!
    ui spi = network.spid[p1];
    if (spi == UI_MAX || spi != network.spid[p2])
        return false;
    ui spt = network.superparticles[spi].sptype;
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    if (!network.sptypes[spt].splookup.count(id))
        return false;
    network.sptypes[clone_sptype(spi)].splookup[id] = interaction;
    update_skins(p1,p2);
    return true;
}

template<ui dim> void md<dim>::mad_sp_bond(ui p1, ui p2, ui interaction)
{   
    //!
    //! This function creates a new superparticle type and assigns it to the superparticle that <tt>p1</tt> and <tt>p2</tt> belong to.
    //! It then assigns an interaction to the two particles, using element <tt>interaction</tt> from <tt>network.library[]</tt>.
    //! It does not perform any checks.
    //!
    ui spi = network.spid[p1];
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    network.sptypes[clone_sptype(spi)].splookup[id] = interaction;
    update_skins(p1,p2);
}

template<ui dim> bool md<dim>::add_sp_bond(ui p1, ui p2, ui potential, vector<ldf> &parameters)
{   
    //!
    //! This function creates a new superparticle type and assigns it to the superparticle that <tt>p1</tt> and <tt>p2</tt> belong to.
    //! It then adds an interaction between the two particles, using a new interaction of the given type and with the given parameters.
    //! If the two particles are not in the same superparticle, or have an interaction already, it does nothing and returns false.
    //!
    ui spi = network.spid[p1];
    if (spi == UI_MAX || spi != network.spid[p2])
        return false;
    ui spt = network.superparticles[spi].sptype;
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    if (network.sptypes[spt].splookup.count(id))
        return false;
    network.sptypes[clone_sptype(spi)].splookup[id] = add_interaction(potential, parameters);
    update_skins(p1,p2);
    return true;
}

template<ui dim> bool md<dim>::mod_sp_bond(ui p1, ui p2, ui potential, vector<ldf> &parameters)
{   
    //!
    //! This function creates a new superparticle type and assigns it to the superparticle that <tt>p1</tt> and <tt>p2</tt> belong to.
    //! It then modifies the interaction between the two particles, using a new interaction of the given type and with the given parameters.
    //! If the two particles are not in the same superparticle, or do not have an interaction already, it does nothing and returns false.
    //!
    ui spi = network.spid[p1];
    if (spi == UI_MAX || spi != network.spid[p2])
        return false;
    ui spt = network.superparticles[spi].sptype;
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    if (!network.sptypes[spt].splookup.count(id))
        return false;
    network.sptypes[clone_sptype(spi)].splookup[id] = add_interaction(potential, parameters);
    update_skins(p1,p2);
    return true;
}

template<ui dim> void md<dim>::mad_sp_bond(ui p1, ui p2, ui potential, vector<ldf> &parameters)
{   
    //!
    //! This function creates a new superparticle type and assigns it to the superparticle that <tt>p1</tt> and <tt>p2</tt> belong to.
    //! It then assigns an interaction to the two particles, using a new interaction of the given type and with the given parameters.
    //! It does not perform any checks.
    //!
    ui spi = network.spid[p1];
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    network.sptypes[clone_sptype(spi)].splookup[id] = add_interaction(potential, parameters);
    update_skins(p1,p2);
}

template<ui dim> bool md<dim>::rem_sp_bond(ui p1, ui p2)
{   
    //!
    //! This function creates a new superparticle type and assigns it to the superparticle that <tt>p1</tt> and <tt>p2</tt> belong to.
    //! It then removes the interaction between the two particles.
    //! If the two particles are not in the same superparticle, or do not have an interaction, it does nothing and returns false.
    //!
    ui spi = network.spid[p1];
    if (spi == UI_MAX || spi != network.spid[p2])
        return false;
    ui spt = network.superparticles[spi].sptype;
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    if (!network.sptypes[spt].splookup.count(id))
        return false;
    network.sptypes[clone_sptype(spi)].splookup.erase(id);
    update_skins(p1,p2);
    return true;
}

template<ui dim> ui md<dim>::clone_sptype(ui spi)
{   
    //!
    //! This function creates a new superparticle type that is exactly the same as that of superparticle <tt>spi</tt>,
    //! assigns it to the superparticle and returns its index in <tt>network.sptypes[]</tt>.
    //! If superparticle <tt>spi</tt> is the only superparticle of its type, this fuction does nothing
    //! and returns the index of its type in <tt>network.sptypes[]</tt>.
    //!
    ui i, spt = network.superparticles[spi].sptype;
    // Check for uniqueness
    for (i = network.superparticles.size()-1; i < UI_MAX && (i == spi || network.superparticles[i].sptype != spt); i--);
    if (i < UI_MAX)
    {   network.sptypes.push_back(network.sptypes[spt]);
        return network.superparticles[spi].sptype = network.sptypes.size()-1;
    }
    else
        return spt;
}

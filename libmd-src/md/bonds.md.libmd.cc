#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::update_skins(ui p1, ui p2)
{
    ui interaction, s, t, j, k;
    pair<ui,ui> id;
    if ((s = network.spid[p1]) < numeric_limits<ui>::max() && s == network.spid[p2]
        && network.sptypes[t = network.superparticles[s].sptype].splookup.count(id = network.hash(network.superparticles[s].particles[p1], network.superparticles[s].particles[p2])))
        interaction = network.sptypes[t].splookup[id];
    else if (network.lookup.count(id = network.hash(particles[p1].type, particles[p2].type)))
        interaction = network.lookup[id];
    else
        interaction = numeric_limits<ui>::max();
    for (j = network.skins[p1].size()-1; j < numeric_limits<ui>::max() && network.skins[p1][j].neighbor != p2; j--);
    if (j < numeric_limits<ui>::max())
    {   for (k = network.skins[p2].size()-1; network.skins[p2][k].neighbor != p1; k--);
        if (interaction < numeric_limits<ui>::max())
            network.skins[p1][j].interaction = network.skins[p2][k].interaction = interaction;
        else
        {   std::iter_swap(network.skins[p1].begin()+j, network.skins[p1].rbegin());
            network.skins[p1].pop_back();
            std::iter_swap(network.skins[p2].begin()+k, network.skins[p2].rbegin());
            network.skins[p2].pop_back();
        }
    }
    else if (interaction < numeric_limits<ui>::max())
    {   network.skins[p1].push_back(interactionneighbor(p2,interaction));
        network.skins[p2].push_back(interactionneighbor(p1,interaction));
    }
}

template<ui dim> bool md<dim>::add_bond(ui p1, ui p2, ui interaction)
{
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
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
        return false;
    else if(interaction>=network.library.size())
    {
        WARNING("Interaction %d does not exist", interaction);
        return false;
    }
    else if(network.free_library_slots.count(interaction))
    {
        WARNING("Interaction %d was previously removed", interaction);
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
    assign_unique_types(p1, p2);
    network.lookup[network.hash(particles[p1].type,particles[p2].type)]=interaction;
    update_skins(p1,p2);
}

template<ui dim> bool md<dim>::add_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
    {   mad_bond(p1, p2, add_interaction(potential, parameters));
        return true;
    }
    else
        return false;
}

template<ui dim> bool md<dim>::mod_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
        return false;
    else
    {   mad_bond(p1, p2, add_interaction(potential, parameters));
        return true;
    }
}

template<ui dim> void md<dim>::mad_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    mad_bond(p1, p2, add_interaction(potential, parameters));
}

template<ui dim> bool md<dim>::rem_bond(ui p1, ui p2)
{
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
    /* Remove bond-style interaction between particles p1 and p2. does not affect other interactions.
     * Parameter force is false by default: it then checks whether there is a bond to begin with.
     * NOTE: forces p1 and p2 to have unique particle types. Replicates former interactions experienced between
     * p1 or p2 and other particle types. */

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
            map<pair<ui,ui>,ui>::iterator it = (old_type[j]==0 ? network.lookup.begin() : network.lookup.upper_bound(pair<ui,ui>(old_type[j]-1, numeric_limits<ui>::max())));
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
    /* add a spring between two points with specified springconstant and equilibrium length */
    vector<ldf> params = {springconstant, l0};
    add_bond(p1,p2,POT::HOOKEAN,&params);
}

template<ui dim> bool md<dim>::add_sp_bond(ui p1, ui p2, ui interaction)
{   ui spi = network.spid[p1];
    if (spi == numeric_limits<ui>::max() || spi != network.spid[p2])
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
{   ui spi = network.spid[p1];
    if (spi == numeric_limits<ui>::max() || spi != network.spid[p2])
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
{   ui spi = network.spid[p1];
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    network.sptypes[clone_sptype(spi)].splookup[id] = interaction;
    update_skins(p1,p2);
}

template<ui dim> bool md<dim>::add_sp_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{   ui spi = network.spid[p1];
    if (spi == numeric_limits<ui>::max() || spi != network.spid[p2])
        return false;
    ui spt = network.superparticles[spi].sptype;
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    if (network.sptypes[spt].splookup.count(id))
        return false;
    network.sptypes[clone_sptype(spi)].splookup[id] = add_interaction(potential, parameters);
    update_skins(p1,p2);
    return true;
}

template<ui dim> bool md<dim>::mod_sp_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{   ui spi = network.spid[p1];
    if (spi == numeric_limits<ui>::max() || spi != network.spid[p2])
        return false;
    ui spt = network.superparticles[spi].sptype;
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    if (!network.sptypes[spt].splookup.count(id))
        return false;
    network.sptypes[clone_sptype(spi)].splookup[id] = add_interaction(potential, parameters);
    update_skins(p1,p2);
    return true;
}

template<ui dim> void md<dim>::mad_sp_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{   ui spi = network.spid[p1];
    pair<ui,ui> id = network.hash(network.superparticles[spi].particles[p1], network.superparticles[spi].particles[p2]);
    network.sptypes[clone_sptype(spi)].splookup[id] = add_interaction(potential, parameters);
    update_skins(p1,p2);
}

template<ui dim> bool md<dim>::rem_sp_bond(ui p1, ui p2)
{
    ui spi = network.spid[p1];
    if (spi == numeric_limits<ui>::max() || spi != network.spid[p2])
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
{   ui i, spt = network.superparticles[spi].sptype;
    // Check for uniqueness
    for (i = network.superparticles.size()-1; i < numeric_limits<ui>::max() && (i == spi || network.superparticles[i].sptype != spt); i--);
    if (i < numeric_limits<ui>::max())
    {   network.sptypes.push_back(network.sptypes[spt]);
        return network.superparticles[spi].sptype = network.sptypes.size()-1;
    }
    else
        return spt;
}

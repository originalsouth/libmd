#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> bool md<dim>::add_bond(ui p1, ui p2, ui interaction)
{
    if (network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
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
    {   rem_bond(p1, p2, true); // "Remove" bond by force (assign unique types)
        network.lookup[network.hash(particles[p1].type,particles[p2].type)]=interaction;
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
    {   rem_bond(p1, p2, true); // "Remove" bond by force (assign unique types)
        network.lookup[network.hash(particles[p1].type,particles[p2].type)]=interaction;
        return true;
    }
}

template<ui dim> void md<dim>::mad_bond(ui p1, ui p2, ui interaction)
{
    rem_bond(p1, p2, true); // "Remove" bond by force (assign unique types)
    network.lookup[network.hash(particles[p1].type,particles[p2].type)]=interaction;
}

template<ui dim> bool md<dim>::add_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    if (!network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
    {   rem_bond(p1, p2, true); // "Remove" bond by force (assign unique types)
        mad_typeinteraction(particles[p1].type, particles[p2].type, potential, parameters); // Add by force
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
    {   rem_bond(p1, p2, true); // Remove bond
        mad_typeinteraction(particles[p1].type, particles[p2].type, potential, parameters); // Add by force
        return true;
    }
}

template<ui dim> void md<dim>::mad_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    rem_bond(p1, p2, true); // "Remove" bond by force (assign unique types)
    mad_typeinteraction(particles[p1].type, particles[p2].type, potential, parameters); // Add by force
}

template<ui dim> bool md<dim>::rem_bond(ui p1, ui p2, bool force)
{
    /* Remove bond-style interaction between particles p1 and p2. does not affect other interactions.
     * Parameter force is false by default: it then checks whether there is a bond to begin with.
     * NOTE: forces p1 and p2 to have unique particle types. Replicates former interactions experienced between
     * p1 or p2 and other particle types. */

    if (!force)
    {   // Check if there is a bond
        if (network.lookup.count(network.hash(particles[p1].type,particles[p2].type)))
        {   // Remove from skins (if present)
            ui j;
            for (j = network.skins[p1].size()-1; j < numeric_limits<ui>::max() && network.skins[p1][j].neighbor != p2; j--);
            if (j < numeric_limits<ui>::max())
            {   std::iter_swap(network.skins[p1].begin()+j, network.skins[p1].rbegin());
                network.skins[p1].pop_back();
            }
            for (j = network.skins[p2].size()-1; j < numeric_limits<ui>::max() && network.skins[p2][j].neighbor != p1; j--);
            if (j < numeric_limits<ui>::max())
            {   std::iter_swap(network.skins[p2].begin()+j, network.skins[p2].rbegin());
                network.skins[p2].pop_back();
            }
        }
        else
            return false;
    }

    // assign unique types to points
    ui P[2], old_type[2], new_type[2];
    P[0] = p1;
    P[1] = p2;
    ui maxtype = 0, i, j, q;

    for (j = 0; j < 2; j++)
        new_type[j] = old_type[j] = particles[P[j]].type;
    for (i = 0; i < N; i++)
        maxtype = max(maxtype, particles[i].type);
    for (auto it = network.lookup.begin(); it != network.lookup.end(); it++)
        maxtype = max(maxtype, it->first.second);

    for (j = 0; j < 2; j++)
    {   for (i = 0; i < N && (i == P[j] || particles[i].type != old_type[j]); i++); // Check if there is at least one other particle with same type as p1
        if (i < N)
        {   // P[j] does not have a unique type; reassign.
            new_type[j] = ++maxtype; // use one greater than largest used particle type
            set_type(P[j], new_type[j]);
        }
    }

    for (j = 0; j < 2; j++)
        if (new_type[j] != old_type[j])
        {   // clone previously defined interactions
            map<pair<ui,ui>,ui>::iterator it, end1 = network.lookup.lower_bound(make_pair(old_type[j],0));
            for (it = network.lookup.begin(); it != end1; it++)
                if (it->first.second == old_type[j] && (q = it->first.first) != new_type[1-j])
                    mad_typeinteraction(new_type[j], q, it->second);
            for (; it != network.lookup.end() && it->first.first == old_type[j]; it++)
                if ((q = it->first.second) != new_type[1-j])
                    mad_typeinteraction(new_type[j], q, it->second);
        }

    rem_typeinteraction(new_type[0], new_type[1]); // removes any old interaction between the unique ids new_p1type and new_p2type
    return true;
}

template<ui dim> void md<dim>::add_spring(ui p1, ui p2, ldf springconstant, ldf l0)
{
    /* add a spring between two points with specified springconstant and equilibrium length */
    vector<ldf> params = {springconstant, l0};
    add_bond(p1,p2,POT::POT_HOOKEAN,&params);
}


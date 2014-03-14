#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::add_bond(ui p1, ui p2, ui itype, vector<ldf> *params)
{
    rem_bond(p1, p2);
    mad_typeinteraction(particles[p1].type, particles[p2].type, itype, params);
}

template<ui dim> void md<dim>::add_spring(ui p1, ui p2, ldf springconstant, ldf l0)
{
    /* add a spring between two points with specified springconstant and equilibrium length */
    vector<ldf> params = {springconstant, l0};
    add_bond(p1,p2,POT::POT_HOOKEAN,&params);
}

template<ui dim> bool md<dim>::share_bond(ui p1, ui p2)
{
    /* Check whether particles p1 and p2 share a bond.
     *  NOTE: Does not take indexing into account, i.e. does
     *  not check if particles are within interaction range.
     *  */


    // 1. Do the particles have unique types?
    if (network.usedtypes[particles[p1].type].size() > 1 || network.usedtypes[particles[p2].type].size() > 1) return false;

    // 2. Do the unique types have an interaction entry?
    pair<ui,ui> id=network.hash(particles[p1].type,particles[p2].type);
    if(network.lookup.find(id)==network.lookup.end()) return false;

    return true;
}

template<ui dim> bool md<dim>::rem_bond(ui p1, ui p2)
{
    /* remove bond-style interaction between particles p1 and p2. does not affect other interactions. */
    //if (!share_bond(p1, p2)) return false;
    //return rem_typeinteraction(particles[p1].type, particles[p2].type);

    /* add a 'bond' i.e. a specific interaction between two particles, of type itype and with parameter params */
    /* NOTE: forces p1 and p2 to have unique particle types. Replicates former interactions experienced between
     * p1 or p2 and other particle types. */

    // assign unique types to points
    ui P[2];
    P[0] = p1;
    P[1] = p2;
    ui old_type[2], new_type[2];
    ui j;
    for (j = 0; j < 2; j++)
        new_type[j] = old_type[j] = particles[P[j]].type;
    
    ui maxtype = 0;
    ui i;
    for (i = 0; i < N; i++)
        maxtype = max(maxtype, particles[i].type);
    for (auto it = network.lookup.begin(); it != network.lookup.end(); it++)
        maxtype = max(maxtype, it->second);
    
    for (j = 0; j < 2; j++)
    {
        for (i = 0; i < N && (i == P[j] || particles[i].type != old_type[j]); i++); // Check if there is at least one other particle with same type as p1
        if (i < N)
        {   // P[j] does not have a unique type; reassign.
            new_type[j] = ++maxtype; // use one greater than largest used particle type
            set_type(P[j], new_type[j]);
        }
    }
    
    for (j = 0; j < 2; j++)
        if (new_type[j] != old_type[j])
        {
            // clone previously defined interactions
            for (map<pair<ui,ui>,ui>::iterator it = network.lookup.begin(); it != network.lookup.end() && it->first.first <= old_type[j]; it++) {
                pair<ui, ui> ipair = it->first;
                if (ipair.first == old_type[j] || ipair.second == old_type[j])
                {   ui q = (ipair.first == old_type[j] ? ipair.second : ipair.first);
                    if (q != new_type[1-j])
                    {   interactiontype old_interaction = network.library[network.lookup[network.hash(old_type[j],q)]];
                        vector<ldf> old_params = old_interaction.parameters;
                        add_typeinteraction(new_type[j], q, old_interaction.potential, &old_params);
                    }
                }
            }
        }


    // now remove the interaction
    rem_typeinteraction(new_type[0], new_type[1]); // removes any old interaction between the unique ids new_p1type and new_p2type

    return true;
}

template<ui dim> bool md<dim>::mod_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    /* modify bond-style interaction between particles p1 and p2. does not affect other interactions. */
    if (!share_bond(p1, p2)) return false;
    return mod_typeinteraction(particles[p1].type, particles[p2].type, potential, parameters);
}

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

    set<ui> partners_of_p1, partners_of_p2;

    // assign unique types to points
    ui old_p1type = particles[p1].type;
    ui new_p1type = old_p1type;
    
    ui maxtype = 0;
    ui i;
    for (i = 0; i < N; i++)
        maxtype = max(maxtype, particles[i].type);
    for (auto it = network.lookup.begin(); it != network.lookup.end(); it++)
        maxtype = max(maxtype, it->second);
    
    for (i = 0; i < N && (i == p1 || particles[i].type != old_p1type); i++); // Check if there is at least one other particle with same type as p1
    if (i < N) {
        // p1 does not have a unique type; reassign.
        new_p1type = ++maxtype; // use one greater than largest used particle type
        set_type(p1, new_p1type);

        // keep track of previously defined interactions
        for (map<pair<ui,ui>,ui>::iterator it = network.lookup.begin(); it != network.lookup.end() && it->first.first <= old_p1type; it++) {
            pair<ui, ui> ipair = it->first;
            if (ipair.first == old_p1type) partners_of_p1.insert(ipair.second);
            else if (ipair.second == old_p1type) partners_of_p1.insert(ipair.first);
        }
    }

    ui old_p2type = particles[p2].type;
    ui new_p2type = old_p2type;

    for (i = 0; i < N && (i == p2 || particles[i].type != old_p2type); i++); // Check if there is at least one other particle with same type as p1
    if (i < N) {
        // p2 does not have a unique type; reassign.
        new_p2type = ++maxtype; // use one greater than largest used particle type
        set_type(p2, new_p2type);

        // keep track of previously defined interactions
        for (map<pair<ui,ui>,ui>::iterator it = network.lookup.begin(); it != network.lookup.end() && it->first.first <= old_p2type; it++) {
            pair<ui, ui> ipair = it->first;
            if (ipair.first == old_p2type) partners_of_p2.insert(ipair.second);
            else if (ipair.second == old_p2type) partners_of_p2.insert(ipair.first);
        }
    }

    // now add the interaction
    rem_typeinteraction(new_p1type, new_p2type); // removes any old interaction between the unique ids new_p1type and new_p2type
    //add_typeinteraction(new_p1type, new_p2type, itype, params);

    // loop through previously defined interactions and clone them so that they are preserved
    if (partners_of_p1.size() > 0) {
        for (set<ui>::iterator it = partners_of_p1.begin(); it != partners_of_p1.end(); it++) {
            if (*it != new_p2type) {
                interactiontype old_interaction = network.library[network.lookup[network.hash(old_p1type,*it)]];
                vector<ldf> old_params = old_interaction.parameters;
                add_typeinteraction(new_p1type, *it, old_interaction.potential, &old_params);
            }
        }
    }
    if (partners_of_p2.size() > 0) {
        for (set<ui>::iterator it = partners_of_p2.begin(); it != partners_of_p2.end(); it++) {
            if (*it != new_p1type) {
                interactiontype old_interaction = network.library[network.lookup[network.hash(old_p2type,*it)]];
                vector<ldf> old_params = network.library[network.lookup[network.hash(old_p2type,*it)]].parameters;
                add_typeinteraction(new_p2type, *it, old_interaction.potential, &old_params);
            }
        }
    }
    return true;
}

template<ui dim> bool md<dim>::mod_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    /* modify bond-style interaction between particles p1 and p2. does not affect other interactions. */
    if (!share_bond(p1, p2)) return false;
    return mod_typeinteraction(particles[p1].type, particles[p2].type, potential, parameters);
}

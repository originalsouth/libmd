Interactions 
============


Interactions in libmd               {#md-interactions}
=====================

Interactions are of two types: pair and forcetypes.

Pair potentials                     {#md-pairpotentials}
---------------

Pairwise interactions are the primary means to define interactions among 
particles in libmd. These are forces due to pair potentials that are a 
function only of the distance between pair of particles. 

TODO: add detail here (Jayson)

Each <tt>md<dim></tt> instance has a library of 
interactions, stored in <tt>network.library[]</tt>.  


### Bonds
A bond is a pair interaction specific to two particles. If two particles share 
a bond, each has a unique particle type that is not shared with any 
other particle in the system. As a result, the particles could have a unique 
interaction that is not shared by any other pair of particles in the system. 
(This is however not a requirement -- a bond could be created between 
particles with an interaction type that is also present between other particle 
types.)

Bonds provide a framework to assign specific interactions to pairs of 
particles based on the particle ID rather than the particle type. For 
instance, a disordered spring network with every spring having a unique bond 
length or spring constant can be implemented by adding bonds among connected 
particles, each of which is a Hookean interaction with the appropriate 
parameters. Many functions exist to automate the task of adding, removing, or 
modifying bonds between pairs of particles without the user having to keep 
track of particle type assignments. These are summarized here. All such 
functions take two particle indices as their first two arguments. The order of 
the particle indices is unimportant.

Functions to create/modify bonds fall into two classes. The first class
assigns a predefined interaction type, referenced by interaction ID, to the 
particle pair:

- md<dim>::add_bond(ui pi, ui p2, ui interaction)
- md<dim>::mod_bond(ui pi, ui p2, ui interaction)
- md<dim>::mad_bond(ui pi, ui p2, ui interaction)

The second class creates a new interaction from a specified potential type 
and a parameter list, and assigns this newly created interaction type to the 
pair of particles:

- md<dim>::add_bond(ui pi, ui p2, ui potential, vector<ldf> *parameters)
- md<dim>::mod_bond(ui pi, ui p2, ui potential, vector<ldf> *parameters)
- md<dim>::mad_bond(ui pi, ui p2, ui potential, vector<ldf> *parameters)

The function md<dim>::rem_bond removes any interaction 
(including a pair interaction of non-bond type) between a specified pair of 
particles. A special convenience function md<dim>::add_spring enables easy 
creation of harmonic springs of specified spring constant and rest length 
between specific particle pairs.

Note that calling functions to add, remove, or modify bonds between particle 
pairs typically changes the particle types associated with each member of the 
pair.

### Bonds within superparticles
TODO: add after Thomas has documented the sp_bond functions in 
bonds.md.libmd.cc


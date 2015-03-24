Interactions 
============


Interactions in libmd               {#md-interactions}
=====================

Interactions are of two types: pair potentials and forcetypes. In molecular 
dynamics, the <em>forces</em> on each particle are computed at each time step 
as a function of neighbouring particles, external constraints, etc. Pair 
porentials are a concise way of defining central forces that can be written as
a gradient of a potential energy, whereas forcetypes encompass all other 
possibilities.

Conservative forces                     {#md-pairpotentials}
-------------------

Molecular dynamics simulations often include conservative forces between 
pairs of particles, which depend only on the distance between the two particles
and act along the vector joining them. They can be written as the gradient of 
a potential energy \f$V(r)\f$ that depends on the separation distance. If \f$\mathbf{r}_i\f$ and 
\f$\mathbf{r}_j\f$ are the position vectors of two particles, then the 
pairwise central forces on each particle are
\f[
    \mathbf{F}_i = -\nabla_{\mathbf{r}_i} V(|\mathbf{r}_i-\mathbf{r}_j|) = 
    -\mathbf{F}_j
\f]

\c libmd uses the framework of <em>particle types</em> to define pairwise
interactions. Every \ref particle instance has an associated \c type, an <tt>unsigned 
int</tt> stored in the \c particle.type variable. Interactions are defined 
between pairs of particle types. For instance, a system of identical particles 
interacting via a Yukawa potential can be set up by assigning the same 
particle type, say \c 0, to all particles in the system, and defining a 
Yukawa-type interaction between the pair of particle types <tt>(0,0)</tt>.

The particle type is set while creating a \ref particle instance via the \c 
ptype argument of the particle<dim>::particle constructor, and can be updated 
by calling the md<dim>::set_type(ui p, ui newtype) function. 

Defining a pairwise interaction between two particle types requires making a 
distinction between <em>pair potentials</em> and <em>interactions</em>:

- A <em>pair potential</em> is the definition of a potential function 
\f$V(r)\f$ that may depend on some parameters in addition to the particle 
separation \f$r\f$. The md<dim>::v structure stores information about pair 
potentials in the simulation. See the subsection on [pair 
potentials](#md-pairpotentialdef) for more information.

- An <em>interaction</em> is a combination of a pair potential and a specific 
set of parameters. The md<dim>::network structure stores information 
about the interactions deriving from pair potentials in the system. See the 
subsection on [interactions](#md-interactiondef) for more information. 

### Pair potentials                             {#md-pairpotentialdef}
Pair potential functions are defined outside the \c md<dim>() structure, and 
added to the simulation using function pointers. A pair potential function takes two 
arguments: the separation distance \c r and a pointer <tt>vector<ldf> 
*params</tt> to a vector of <tt>float</tt>s that contains the parameters 
needed to compute the interaction, and returns the potential energy.
\c libmd uses [automatic differentiation](http://en.wikipedia.org/wiki/Automatic_differentiation)
implemented in autodiff.libmd.cc to calculate the forces from the potential 
definition. This mandates that the potential function be defined with the 
first argument (the distance \c r) and the return value having a templated type:
\code{.cpp}
    template<class X> X my_potential(X r,vector<ldf> *parameters) {...}
\endcode
However, the function itself can be written by treating  \c r and the return 
value as \c ldf variables. See  potentials.libmd.cc for example definitions of potential functions.

Once the potential function has been defined, a function pointer pointing to 
it is added to the list of potential functions 
<tt>md<dim>::v.potentials</tt> through a call to the pairpotentials::add() 
member function of md<dim>::v. Each added potential is given a unique index, 
which is the return value of <tt>md<dim>::v.add()</tt>. 
The 
following command adds the potential function \c my_potential() defined 
above to an <tt>md\<dim\></tt> instance named \c sys, and stores the 
associated index in \c pidx:
\code{.cpp}
    md<dim> sys;
    ...
    ui pidx = sys.v.add(my_potential<dual>);
\endcode
(The <tt>\<dual\></tt> template argument is needed for the automatic 
differentiation system.) This index is 
then used to define interactions between pairs of particle types through 
calls to md<dim>::add_interaction and md<dim>::add_typeinteraction functions 
catalogued in the subsection on [interactions](#md-interactiondef). 



#### Predefined pair potentials
\c libmd comes with some predefined pair potentials that can be combined with
any choices of parameters to easily define interactions between particle types.
A complete listing is available in the documentation for the potentials.libmd.cc
file. The potentials are pre-loaded into the <tt>md<dim>::v.potentials</tt>
list of every \ref md instance, and are to be referenced by their positions in this list. For ease of 
indexing, a global enum structure #POT has been defined, with 
<tt>POT::\<potential name\></tt> providing the appropriate pair potential 
index. For instance, #POT::COULOMB is the index of the predefined COULOMB() 
potential.


### Interactions                                {#md-interactiondef}
Once a potential has been added to the \ref md instance, it can be combined 
with specific values of parameters to create distinct \a interactions, which 
are instances of the \ref interactiontype structure. Each 
\ref md instance has a library of interactions, stored in
<tt>md<dim>::network.library[]</tt> which is a vector of \ref interactiontype 
instances. Different interactions are indexed by their position in this vector.

Entries in the interaction library are added, modified and removed using the 
following set of functions:

- md<dim>::add_interaction(ui pidx, vector<ldf> *parameters) creates an 
interaction from the pair potential indexed by \c pidx with the given 
parameters. The return value is the index of the interaction in 
<tt>md<dim>::network.library[]</tt>, which we will call \c iidx.
- md<dim>::mod_interaction(ui iidx, ui pidx,vector<ldf> *parameters) replaces 
the potential function and parameters of the interaction indexed by \c iidx.
- md<dim>::rem_interaction(ui iidx) removes the interaction indexed by 
\c iidx.

The following set of functions assigns a particular interaction indexed by \c 
iidx to act between particles of types \c type1 and \c type2:

- md<dim>::add_typeinteraction(ui type1, ui type2, ui iidx)
- md<dim>::mod_typeinteraction(ui type1, ui type2, ui iidx)
- md<dim>::mad_typeinteraction(ui type1, ui type2, ui iidx)

The functions differ in their behaviour when an interaction between the pair 
of types has or has not been defined.

There is also a set of functions which allows the creation of an interaction and 
its assignment to act between particle types \c type1 and \c type2 in a single 
command:

- md<dim>::add_typeinteraction(ui type1, ui type2, ui pidx, vector<ldf> *parameters)
- md<dim>::mod_typeinteraction(ui type1, ui type2, ui pidx, vector<ldf> *parameters)
- md<dim>::mad_typeinteraction(ui type1, ui type2, ui pidx, vector<ldf> *parameters)


The following code snippets are therefore equivalent:
\code{.cpp}
    /* Uses predefined Hookean spring potential, indexed by POT::HOOKEAN */
   
    vector<ldf> params = {1.0,2.0};     // Vector of two parameters for a Hookean spring 
                                        // (spring constant and rest length)
                                        
    ui onespring = md<dim>::add_interaction(POT::HOOKEAN, &params); 
                                        // Create a Hookean interaction with the given parameters
                                        // and store its index in onespring
                                        
    md<dim>::add_typeinteraction(2,8,onespring);  
                                        // Add the interaction indexed by onespring between 
                                        // particle types 2 and 8
\endcode

\code{.cpp}
    /* Uses predefined Hookean spring potential, indexed by POT::HOOKEAN  */
   
    vector<ldf> params = {1.0,2.0};     // Vector of two parameters for a Hookean spring
                                        // (spring constant and rest length)
                                        
    md<dim>::add_typeinteraction(2,8,POT::HOOKEAN,&params);  
                                        // Add a Hookean interaction between particles of 
                                        // type 2 and type 8 with defined params
\endcode

#### Cutoff radius
An important parameter in limiting unnecessary computations of  pair potentials 
is the *cutoff radius*, a value of the particle separation beyond which the 
pair potential and resultant force are assumed to be zero. Each interaction 
(i.e. each \ref interactiontype instance in \c md<dim>::network.library[]) has 
a unique cutoff radius stored in the \c rco member variable. By default, this 
is set to be equal to the value of \c md<dim>::network.rco, which is the case 
for the examples listed above. However, every 
function that creates a new interaction also has a version in which the cutoff 
radius can be explicitly specified as an additional parameter, as follows:

- md<dim>::add_interaction(ui pidx, ldf rco, vector<ldf> *parameters)
- md<dim>::mod_interaction(ui iidx, ui pidx, ldf rco, vector<ldf> *parameters)
- md<dim>::add_typeinteraction(ui type1, ui type2, ui pidx, ldf rco, vector<ldf> *parameters)
- md<dim>::mod_typeinteraction(ui type1, ui type2, ui pidx, ldf rco, vector<ldf> *parameters)
- md<dim>::mad_typeinteraction(ui type1, ui type2, ui pidx, ldf rco, vector<ldf> *parameters)


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
assigns a predefined interaction type, referenced by its index \c iidx in 
<tt>md<dim>::network.library[]</tt>, to the particle pair:

- md<dim>::add_bond(ui p1, ui p2, ui iidx)
- md<dim>::mod_bond(ui p1, ui p2, ui iidx)
- md<dim>::mad_bond(ui p1, ui p2, ui iidx)

The second class creates a new interaction from a specified potential type 
(referenced by its index \c pidx in <tt>md<dim>::v.potentials[]</tt>)
and a parameter list, and assigns this newly created interaction type to the 
pair of particles:

- md<dim>::add_bond(ui p1, ui p2, ui pidx, vector<ldf> *parameters)
- md<dim>::mod_bond(ui p1, ui p2, ui pidx, vector<ldf> *parameters)
- md<dim>::mad_bond(ui p1, ui p2, ui pidx, vector<ldf> *parameters)

The function md<dim>::rem_bond removes any interaction 
(including a pair interaction of non-bond type) between a specified pair of 
particles. A special convenience function md<dim>::add_spring enables easy 
creation of harmonic springs of specified spring constant and rest length 
between specific particle pairs.

**Note:** calling functions to add, remove, or modify bonds between particle 
pairs typically **changes the particle types** associated with each member of the 
pair. Therefore, it is best to use bond functions after all particle type and 
potential assignments have been completed (the bond functions will preserve 
type interactions through the particle reassignment). 


### Interactions and bonds within superparticles
Interactions and bonds have a special meaning in the context of 
[superparticles](#md-superparticles). 

A pairwise *interaction* can be defined 
between two particle subtypes \c p1 and \c p2 within a superparticle type 
\c spt; this interaction will then be present in all copies of \c spt in the 
system. The pair potentials and interactions are shared with the 
ordinary particles (i.e. the same libraries and indices are used), but the \c 
add/mod/mad_typeinteraction() functions described in the 
section on [interactions](#md-interactiondef) are replaced by the following 
functions:

- md<dim>::add_sp_interaction(ui spt, ui p1, ui p2, ui iidx)
- md<dim>::add_sp_interaction(ui spt, ui p1, ui p2, ui pidx, vector<ldf> *parameters) 
- md<dim>::add_sp_interaction(ui spt, ui p1, ui p2, ui pidx, ldf rco, vector<ldf> *parameters) 

(and similarly for \c mod_sp_interaction and \c mad_sp_interaction).

A *bond* within a superparticle is created by specifying two *ordinary* 
particle indices (not subtypes) \c p1 and \c p2 which must belong to the 
same superparticle. An interaction is then created between these two 
particles, which only exists in that particular superparticle instance, making 
it distinct from all other superparticles in the system. Successful 
superparticle bond creation always 
gives rise to a new superparticle type specific to the 
superparticle instance that has been targeted. See the following functions for 
more details:

- md<dim>::add_sp_bond(ui p1, ui p2, ui iidx)
- md<dim>::add_sp_bond(ui p1, ui p2, ui pidx, vector<ldf> *parameters) 

(and similarly for \c mod_sp_bond and \c mad_sp_bond).

Finally the functions md<dim>::rem_sp_interaction() and md<dim>::rem_sp_bond() 
remove superparticle interactions and bonds respectively. 





Non-conservative forces
---------------------

Many applications of molecular dynamics involve forces that cannot be defined 
as gradients of [pair potentials](#md-pairpotentialdef); e.g. dissipative drag 
on a particle. We use the \ref forcetype framework to define such forces on 
particles. This framework is completely general, and can be used to define any 
number of operations (not restricted to force calculation) on individual particles that can be calculated as a function 
of particle positions, velocities, or external parameters. 

As with conservative forces, a distinction is to be made between *external 
force functions* and *forcetypes*, which are the equivalent of *pair 
potentials* and *type interactions* respectively. Information about external 
force functions is stored in the md<dim>::f structure, whereas the forcetypes 
deriving from these force functions  are stored in the   md<dim>::network 
structure (specifically \c md<dim>::network.forces and 
\c md<dim>::network.forcelibrary). 

### External force functions

The \ref forcetype framework requires operations on particles to be represented 
as *external force functions*, which live outside the \c md<dim>() structure 
and fits the following prototype:
\code{.cpp}
    template <ui dim> void my_external_force(ui i, vector<ui> *partners, vector<ldf> *parameters, void *sys) {...}
\endcode
where \c i is a particle index, \c partners points to a vector of particle 
indices which might influence particle \c i, and \c parameters can be 
specified. During execution, a pointer to the  
\ref md object itself is passed to the external  
force function as the last argument. Under normal use, a call to \c my_external_force() 
function would compute a force on particle \c i that is a function 
of the positions and velocities of \c i and \c partners, and update the force 
vector \c ((md<dim>*) sys)->particles[i].F with this force. 
However, the external force function is not restricted to calculating and 
updating forces, but could modify the \ref md object in any way. 


The macro \ref SYS eases the referencing of member variables within \c sys 
within the \c my_external_force() function. 
For instance,  \c ((md<dim>*) sys)->particles[i].F can be rewritten as \c 
SYS->particles[i].F. 

Once an external force function has been defined, a function pointer pointing to 
it is added to the list of external force functions 
<tt>md<dim>::f.extforces</tt> through a call to the externalforces::add() 
member function of md<dim>::f. Each added external force function is given a unique index, 
which is the return value of \c md<dim>::f.add(). The 
following command adds the function \c my_external_force() to an <tt>md\<dim\></tt> instance named \c sys, and stores the 
associated index in \c fidx:
\code{.cpp}
    md<dim> sys;
    ...
    ui fidx = sys.f.add(my_external_force<dim>);
\endcode
*Particular* forcetypes are 
obtained by combining an external force function indexed by \c fidx with a specific set of 
parameters, using the md<dim>::add_forcetype function 
catalogued in the [following subsection](#md-forcetypedef). 

#### Predefined external force functions  

\c libmd comes with two predefined external force functions, which illustrate 
the concept. 

- The \ref DAMPING() function updates \c SYS->particles[i].F with a viscous drag force that is a 
function solely of the velocity of particle \c i. The corresponding force 
index is \c EXTFORCE::DAMPING.
- The \ref DISSIPATION() function updates \c SYS->particles[i].F with a force proportional to 
the velocity difference between \c i and each particle index in \c particles. 
The function is indexed by EXTFORCE::DISSIPATION.


### Forcetypes                              {#md-forcetypedef}
An external force function indexed by \c fidx is combined 
with specific values of parameters and a specific list of particles associated 
with each particle to create distinct \a forcetypes, which 
are instances of the \ref forcetype structure. These are stored in
<tt>md<dim>::network.forcelibrary[]</tt>, a vector of \ref forcetype 
instances. Different forcetypes are indexed by their position in this vector.

Entries in the forcetype library are added, modified and removed using the 
following functions, which returns a \c bool indicating the success of the 
operation:

- md<dim>::add_forcetype(ui fidx, vector<vector<ui>> *partnerlist, vector<ldf> 
*parameters) creates a forcetype from the external force function indexed by \c fidx with the given 
parameters. The \c partnerlist is either empty, or points to a list of lists 
of partner indices, so that \c (*partnerlist)[i] contains the partners of 
particle \c i to be passed to the external force function. The return value is 
the index of the forcetype in  
<tt>md<dim>::network.forcelibrary[]</tt>, which we will call \c ftype.
- md<dim>::mod_forcetype(ui ftype, ui fidx, vector<vector<ui>> *partnerlist, 
vector<ldf> *parameters) replaces the external force function and partner list 
of the forcetype indexed by \c ftype.
- md<dim>::rem_forcetype(ui ftype) removes the interaction indexed by 
\c ftype.


Once a forcetype entry has been created, it still needs to be assigned to a 
particle to influence the dynamics of that particle. Forcetype assignment is handled by the 
following set of functions:

- md<dim>::assign_forcetype(ui i, ui ftype) assigns \c ftype to particle \c i.
- md<dim>::assign_all_forcetype(ui ftype) assigns \c ftype to all particles.
- md<dim>::unassign_forcetype(ui ftype) and 
md<dim>::unassign_all_forcetype(ui ftype, ui i) remove the assignment of \c ftype 
from one or all particles.

Finally, the function md<dim>::clear_all_assigned_forcetype() is useful to 
clear all assignments to forcetypes from all particles, although it does not 
remove the forcetypes themselves. 

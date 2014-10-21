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
Pair potential functions are defined outside the md<dim>() structure, and 
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
which is the return value of <tt>md<dim>::v.add()</tt>. This index is 
then used to define interactions between pairs of particle types through 
calls to md<dim>::add_interaction and md<dim>::add_typeinteraction functions 
catalogued in the subsection on [interactions](#md-interactiondef). The 
following command adds the potential function \c my_potential() defined 
above to an <tt>md\<dim\></tt> instance named \c sys, and stores the 
associated index in \c pidx:
\code{.cpp}
    md<dim> sys;
    ...
    ui pidx = sys.v.add(my_potential<dual>);
\endcode
The <tt>\<dual\></tt> template argument is needed for the automatic 
differentiation system. 


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

- md<dim>::add_interaction(ui pidx,vector<ldf> *parameters) creates an 
interaction from the pair potential indexed by \c pidx with the given 
parameters. The return value is the index of the interaction in <tt>md<dim>::network.library[]</tt>.
- md<dim>::mod_interaction(ui iidx, ui pidx,vector<ldf> *parameters) replaces 
the potential function and parameters of the interaction indexed by \c iidx.
- md<dim>::rem_interaction(ui interaction) removes the interaction indexed by 
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
   
    vector<ldf> params = {1.0,2.0};                          // Vector of two parameters for a Hookean spring -- spring constant and rest length
    ui onespring = md<dim>::add_interaction(POT::HOOKEAN, &params); // Create a Hookean interaction with the given parameters and store its index in onespring
    md<dim>::add_typeinteraction(2,8,onespring);  // Add a Hookean interaction between particles of type 2 and type 8 with defined params
\endcode

\code{.cpp}
    /* Uses predefined Hookean spring potential, indexed by POT::HOOKEAN  */
   
    vector<ldf> params = {1.0,2.0};                          // Vector of two parameters for a Hookean spring -- spring constant and rest length
    md<dim>::add_typeinteraction(2,8,POT::HOOKEAN,&params);  // Add a Hookean interaction between particles of type 2 and type 8 with defined params
\endcode



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
<tt>md<dim>::network.library[]</tt>, to the 
particle pair:

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

Note that calling functions to add, remove, or modify bonds between particle 
pairs typically changes the particle types associated with each member of the 
pair.

Non-conservative forces
-----------------------



### Bonds within superparticles
TODO: add after Thomas has documented the sp_bond functions in 
bonds.md.libmd.cc


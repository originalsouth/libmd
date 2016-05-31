///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of libmd HEADER file                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef libmd_h
#define libmd_h

#if __cplusplus < 201103L
#error "C++11 not detetected: libmd requires C++11 to work (update compiler)."
#endif

#include <cstdio>                                                       //< Standard input output (faster than IOstream and also threadsafe) (C)
#include <cstdlib>                                                      //< Standard library (C)
#include <cmath>                                                        //< Standard math library  (C)
#include <cstring>                                                      //< Memcpy and memmove support (C)
#include <vector>                                                       //< Vector support (C++)
#include <map>                                                          //< Map support (C++)
#include <unordered_map>                                                //< Map support (C++)
#include <list>                                                         //< List support (C++)
#include <stack>                                                        //< Stack support (C++)
#include <set>                                                          //< Set support (C++)
#include <unordered_set>                                                //< Unordered set support (C++11)
#include <utility>                                                      //< Pair support (C++)
#include <limits>                                                       //< Limits of types (C++)
#include <algorithm>                                                    //< Algorithm support (C++)
#include <functional>                                                   //< Functional support (C++11)
#include <chrono>                                                       //< Timing support (C++11)

#ifdef FE
#include <fenv.h>                                                       //< Floating point exception handling (C)
#endif

#include "libmd-src/macros.libmd.h"                                     //< Implementation of libmd (preproccessor) macros

#ifdef LIBMD__LONG_DOUBLE__                                             //< user wants to use long double precision
typedef long double ldf;                                                //< long double is now aliased as ldf
#define F_LDF "%.17Lg"                                                  //< defines the printf format for ldf as long double
#define F_LDFs "%Lg"                                                    //< defines the scanf format for ldf as long double
#elif defined LIBMD__FLOAT__                                            //< user wants to use float precision
typedef float ldf;                                                      //< float is now aliased as ldf
#define F_LDF "%.17g"                                                   //< defines the printf format for ldf as float
#define F_LDFs "%g"                                                     //< defines the scanf format for ldf as float
#else                                                                   //< user wants to use double precision (default)
typedef double ldf;                                                     //< double is now aliased as ldf
#define F_LDF "%.17lg"                                                  //< defines the printf format for ldf as double
#define F_LDFs "%lg"                                                    //< defines the scanf format for ldf as double
#endif

typedef unsigned int ui;                                                //< unsigned int is now aliased as ui
#define F_UI "%u"                                                       //< defines the printf format for ui
typedef unsigned char uc;                                               //< unsigned char is now aliased as uc
#define F_UC "%c"                                                       //< defines the printf format for uc

#include "libmd-src/enums.libmd.h"                                     //< Implementation of enums defined in libmd

const ui UI_MAX=std::numeric_limits<ui>::max();                         //< UI_MAX is defined as the largest ui (unsigned integer)

//These functions defined outside of a libmd structure
void __libmd__info();                                                   ///< Basic libmd comilation info
ldf TicToc();                                                           ///< High precision timer
template<ui dim> ldf dotprod (ldf A[], ldf B[]);

template<ui dim> void BCOND_NONE(ui d,ui i,void *sys);
template<ui dim> void BCOND_NONE(ui d,ldf x[dim],void *sys);
template<ui dim> void BCOND_PERIODIC(ui d,ui i,void *sys);
template<ui dim> void BCOND_PERIODIC(ui d,ldf x[dim],void *sys);
template<ui dim> void BCOND_HARD(ui d,ui i,void *sys);
template<ui dim> void BCOND_HARD(ui d,ldf x[dim],void *sys);
template<ui dim> void BCOND_BOXSHEAR(ui d,ui i,void *sys);
template<ui dim> void BCOND_BOXSHEAR(ui d,ldf x[dim],void *sys);

template<class X> X COULOMB(X r,std::vector<ldf> &parameters);          ///< Coulomb potential functions
template<class X> X YUKAWA(X r,std::vector<ldf> &parameters);           ///< Yukawa potential functions
template<class X> X HOOKEAN(X r,std::vector<ldf> &parameters);          ///< Hookean potential functions
template<class X> X LJ(X r,std::vector<ldf> &parameters);               ///< The famous Lennard-Jones potential functions
template<class X> X MORSE(X r,std::vector<ldf> &parameters);            ///< Morse potential functions
template<class X> X FORCEDIPOLE(X r,std::vector<ldf> &parameters);      ///< Force dipole potential functions
template<class X> X HOOKEANFORCEDIPOLE(X r,std::vector<ldf> &parameters); ///< Hookean force dipole potential functions
template<class X> X ANHARMONICSPRING(X r,std::vector<ldf> &parameters); ///< Anharmonic spring potential functions

template<ui dim> void DAMPING(ui i,std::vector<ui> &particles,std::vector<ldf> &parameters,void *sys); ///< Damping external force functions
template<ui dim> void DISSIPATION(ui i,std::vector<ui> &particles,std::vector<ldf> &parameters,void *sys); ///< Dissipation external force functions

ldf kdelta(ui i,ui j);                                                  ///< Kronecker delta function

template<class X,ui dim> X FLATSPACE(X x[dim],std::vector<ldf> &param); ///< Flat space Monge function
template<class X,ui dim> X GAUSSIANBUMP(X x[dim],std::vector<ldf> &param);///< Gaussian bump Monge function
template<class X,ui dim> X EGGCARTON(X x[dim],std::vector<ldf> &param); ///< Egg carton bump Monge function
template<class X,ui dim> X MOLLIFIER(X x[dim],std::vector<ldf> &param); ///< Mollifier bump Monge function

//Function pointers used by libmd
template<ui dim> using bcondpptr=void (*)(ui d,ui i,void *sys);         ///< Function pointer to particle bcond function is now called perodicitypptr
template<ui dim> using bcondxptr=void (*)(ui d,ldf x[dim],void *sys);   ///< Function pointer to position bcond function is now called perodicityxptr
template<class X> using potentialptr=X (*)(X,std::vector<ldf> &);       ///< Function pointer to potential functions is now called potentialptr
template<ui dim> using extforceptr=void (*)(ui,std::vector<ui> &,std::vector<ldf> &,void *); ///< Function pointer to external force functions is now called extforceptr
template<class X,ui dim> using fmpptr=X (*)(X x[dim],std::vector<ldf> &param); ///< Monge patch function pointer
template<ui dim> using hookptr=void (*)(std::vector<ldf> &,void *);     ///< Function pointer to external force functions is now called extforceptr

/// This structure handles errors/warnings/debug levels
struct t_error
{
    ui term_level;                                                      ///< Terminate level for libmd. The default value is 1.
    FILE *error_file;                                                   ///< libmd error output file (default stderr)
    FILE *warning_file;                                                 ///< libmd warning output file (default stderr)
    FILE *debug_1_file;                                                 ///< libmd debug[1] output file (default stdout)
    FILE *debug_2_file;                                                 ///< libmd debug[2] output file (default stdout)
    FILE *debug_3_file;                                                 ///< libmd debug[3] output file (default stdout)
    FILE *debug_timer_file;                                             ///< libmd debug[timer] output file (default stdout)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    t_error();                                                          ///< Constructor
    ~t_error();                                                         ///< Destructor (to close the files)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void set_error_file(const char *fname);                             ///< Sets the error output file
    void set_warning_file(const char *fname);                           ///< Sets the warning output file
    void set_debug_1_file(const char *fname);                           ///< Sets the debug[1] output file
    void set_debug_2_file(const char *fname);                           ///< Sets the debug[2] output file
    void set_debug_3_file(const char *fname);                           ///< Sets the debug[3] output file
    void set_debug_timer_file(const char *fname);                       ///< Sets the debug[timer] output file
    void print_error(char *buffer);                                     ///< Prints a error to the error output file (for internal use)
    void print_warning(char *buffer);                                   ///< Prints a warning to the warning output file (for internal use)
    void print_debug_1(char *buffer);                                   ///< Prints debug[1] message to the debug[1] output file (for internal use)
    void print_debug_2(char *buffer);                                   ///< Prints debug[2] message to the debug[2] output file (for internal use)
    void print_debug_3(char *buffer);                                   ///< Prints debug[3] message to the debug[3] output file (for internal use)
    void print_debug_timer(char *buffer);                               ///< Prints debug[timer] message to the debug[timer] output file (for internal use)
    void terminate(ui term);                                            ///< Terminate if termlevel allows it (for internal use)
} error;

/// This structure contains all the information for a single particle
template<ui dim> struct particle
{
    ldf m;                                                              ///< Mass
    ldf x[dim];                                                         ///< Position
    ldf xp[dim];                                                        ///< Previous particle position
    ldf xsk[dim];                                                       ///< Position when building the skinlist
    ldf dx[dim];                                                        ///< Velocity
    ldf F[dim];                                                         ///< Forces on particle
    ui type;                                                            ///< This particle has type number
    bool fix;                                                           ///< Can this particle move
    bool usepbcond;                                                     ///< Use per particle boundary conditions
    uc pbc[dim];                                                        ///< Particle boundary conditions
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    particle(ldf mass=1.0,ui ptype=0,bool fixed=false,bool pbconded=false); ///< Constructor
};

/// This structure contains information about the simulation box
template<ui dim> struct box
{
    ldf L[dim];                                                         ///< Box size
    bool useLshear;                                                     ///< Use sheared box matrix
    ldf vshear[dim][dim];                                               ///< Shear velocity vshear[i][j] is shear velocity in direction i of boundary with normal in direction j. currently vshear[i][i] != 0 results in undefined behaviour.
    ldf Lshear[dim][dim];                                               ///< Box matrix that is updated at each time step. Used to compute distances for shear, in lieu of simbox.L
    ldf LshearInv[dim][dim];                                            ///< Inverse of Lshear[][]
    uc bcond[dim];                                                      ///< Boundary conditions in different dimensions NONE/PERIODIC/HARD/BOXSHEAR
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    box();                                                              ///< Constructor
    void shear_boundary(ui i, ui j, ldf velocity);                      ///< Set up boundary shear velocity in direction i of boundary with normal direction j
    void skew_boundary(ui i, ui j, ldf displacement);                   ///< Skew the simulation box by moving boundary with normal direction j by amount 'displacement' in direction i
    void invert_box();                                                  ///< Invert the Lshear[][] box matrix
};

template<ui dim> struct bcond
{
    std::vector<bcondpptr<dim>> bcond_p;                                ///< Vector of bcond particle function pointers
    std::vector<bcondxptr<dim>> bcond_x;                                ///< Vector of bcond position function pointers
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bcond();                                                            ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(bcondpptr<dim> p,bcondxptr<dim> x);                          ///< Add bcond functions to their respective vectors
    void operator()(ui d,ui i,void *sys);                               ///< Periodicity operator
    void operator()(ui d,ldf x[dim],void *sys);                         ///< Periodictty overloaded operator
    void operator()(ui k,ui d,ui i,void *sys);                          ///< Periodicity operator
    void operator()(ui k,ui d,ldf x[dim],void *sys);                    ///< Periodictty overloaded operator
};

/// This structure saves the particle type interactions and calculates the the potentials
struct interactiontype
{
    ui potential;                                                       ///< Type of potential
    std::vector<ldf> parameters;                                        ///< Parameters of potential
    ldf rco;                                                            ///< R_cuttoff radius
    ldf vco;                                                            ///< Cuttoff potential energy
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interactiontype(ui ppot,std::vector<ldf> &param,ldf Rco,ldf Vco);   ///< Constructor
};

/// This struct saves the neighboring particle number and the interaction type library number
struct interactionneighbor
{
    ui neighbor;                                                        ///< Neighbor number
    ui interaction;                                                     ///< Interaction number
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interactionneighbor(ui noneighbor,ui nointeraction);                ///< Constructor
};

/// This structure saves the external force functions and calculates them
struct forcetype
{
    ui externalforce;                                                   ///< External force type
    std::vector<std::vector<ui>> particles;                             ///< Interacting particle list
    std::vector<ldf> parameters;                                        ///< Parameters for the external force
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    forcetype(ui noexternalforce,std::vector<ldf> &param);              ///< Constructor
    forcetype(ui noexternalforce,std::vector<std::vector<ui>> &plist,std::vector<ldf> &param); ///< Constructor
};

/// This structure introduces "super_particles" i.e. particles that consist of (sub_)particles
struct superparticle
{
    std::unordered_map<ui,ui> particles;                                ///< Particles in super particles
    std::vector<ui> backdoor;                                           ///< Super particle index to particle id
    ui sptype;                                                          ///< Super particle type
    bool center_bcond;                                                  ///< Use the boundary conditions of its particle calculate about its center of mass
};

/// This structure caries a lookup device for a specific super particle type
struct superparticletype
{
    std::map<std::pair<ui,ui>,ui> splookup;                             ///< This is the interaction lookup device
};

/// This structure stores all interactions and their types
struct interact
{
    bool update;                                                        ///< Should we update the network
    ldf rco;                                                            ///< Default R_cutoff radius
    ldf ssz;                                                            ///< Skin radius
    std::vector<std::vector<ui>> forces;                                ///< List of external forces acting on the particles
    std::vector<forcetype> forcelibrary;                                ///< Library of external forces
    std::vector<std::vector<interactionneighbor>> skins;                ///< Particle skin by index (array of vector)
    std::vector<interactiontype> library;                               ///< This is the interaction library
    std::unordered_set<ui> free_library_slots;                          ///< Stores free library slots
    std::map<std::pair<ui,ui>,ui> lookup;                               ///< This is the interaction lookup device
    std::vector<ui> spid;                                               ///< Super particle identifier array
    std::vector<superparticle> superparticles;                          ///< Actual super particle array
    std::vector<superparticletype> sptypes;                             ///< Super particle type array
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interact();                                                         ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::pair<ui,ui> hash(ui type1,ui type2);                           ///< Hash function
    bool probe(ui type1,ui type2);                                      ///< Check if a typeinteraction exists between two types
};

/// This structure automatically differentiates first order
struct dual
{
    ldf x;                                                              ///< Function value
    ldf dx;                                                             ///< Function derivative value
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dual();                                                             ///< Constructor
    dual(ldf f,ldf fx=0.0);                                             ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dual operator=(dual G);                                             ///< Assign operator
    template<class X> dual operator=(X a);                              ///< Assign foreign type operator
    template<class X> operator X();                                     ///< Cast overload
};

/// This structure takes care of pair potentials (who live outside of the class)
struct pairpotentials
{
    std::vector<potentialptr<dual>> potentials;                         ///< Pair potential vector
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pairpotentials();                                                   ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(potentialptr<dual> p);                                       ///< Add a potentials
    ldf operator()(ui type,ldf r,std::vector<ldf> &parameters);         ///< Pair potential executer
    ldf dr(ui type,ldf r,std::vector<ldf> &parameters);                 ///< Pair potential d/dr executer
};

/// This structure takes care of additional (external) forces acting on particles
template<ui dim> struct externalforces
{
    std::vector<extforceptr<dim>> extforces;                            ///< External forces function container
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    externalforces();                                                   ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(extforceptr<dim> p);                                         ///< Add an external force function
    void operator()(ui type,ui i,std::vector<ui> &particles,std::vector<ldf> &parameters,void *sys); ///< Execute external force function
};

/// This structure defines and saves integration metadata
struct integrators
{
    ldf h;                                                              ///< Timestep size
    uc method;                                                          ///< Type of integration
    ui generations;                                                     ///< Maximum generations of timestep
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    integrators();                                                      ///< Constructor
};

/// This structure is specific for the indexer
template<ui dim> struct indexer
{
    uc method;                                                          ///< Method of indexing
    struct celldatatype
    {
        ui Q[dim];                                                      ///< Number of cells per dimension
        ui nCells;                                                      ///< Total number of cells (= prod(Q))
        ui totNeighbors;                                                ///< Total number of (potential) neighboring cells to check (= (3^d-1)/2)
        ldf CellSize[dim];                                              ///< Length of cell in each dimension
        int (*IndexDelta)[dim];                                         ///< Indices of neighboring cells relative to cell
        std::vector<ui> *Cells;                                         ///< List of particles per cell
        bool *OutsideBox;                                               ///< Indicates for each particle whether it is outside the simbox
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        celldatatype();                                                 ///< Constructor
        ~celldatatype();                                                ///< Destructor
    };
    celldatatype celldata;                                              ///< Cell data object
    struct kdtreedatatype
    {
        ui (*Idx);                                                      ///< Indices of particles, ordered by tree-structure
        ui DivideByDim[30];                                             ///< Dimension to divide by at each recursion level (assuming N <= 2^30)
        ldf (*Pmin)[dim],(*Pmax)[dim];                                  ///< Minimum and maximum value of each coordinate, for every subtree
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        kdtreedatatype();                                               ///< Constructor
        ~kdtreedatatype();                                              ///< Destructor
    };
    kdtreedatatype kdtreedata;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    indexer();                                                          ///< Constructor
};

/// This structure
template<ui dim> struct t_hook
{
    std::vector<hookptr<dim>> hooks;                                    ///<
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(hookptr<dim> p);                                             ///<
    void operator()(ui idx,std::vector<ldf> &parameters,void *sys);     ///<
};

/// This structure
struct hooktype
{
    ui hook;                                                            ///< Hooktype
    std::vector<ldf> parameters;                                        ///< Hook parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    hooktype(ui nohook,std::vector<ldf> &param);                        ///< Constructor
};

/// This structure
template<ui dim> struct hooker
{
    t_hook<dim> hook;                                                   ///<
    std::vector<hooktype> hookers;                                      ///<
};

/// This structure stores some cyclic variables for the variadic functions
template<ui dim> struct variadic_vars
{
    std::vector<ui> vvars;                                              ///< Container of variables
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    variadic_vars();                                                    ///< Initialize variables (set everyting to zero)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void reset();                                                       ///< Reset all variadic_vars to zero
    void reset(ui i);                                                   ///< Reset a certain variadic_vars to zero
    ui operator[](ui i);                                                ///< Rotate and return previous for the ith variable
};

/// This structure stores additional variables
template<ui dim> struct additional_vars
{
    ui noftypedamping=UI_MAX;                                           ///< This variable stores the number of the ftype that damps/drags the system
    bool export_force_calc=false;                                       ///< This variable tells export_force if the forces have been calculated for this output
    bool reindex=true;                                                  ///< This variable tells if the system needs to be reindexed
};

/// This structure defines the molecular dynamics simulation
template<ui dim> struct md
{
    ui N;                                                               ///< Number of particles
    box<dim> simbox;                                                    ///< Simulation box
    bcond<dim> boundary;                                                ///< Boundary conditions functor
    std::vector<particle<dim>> particles;                               ///< Particle array
    interact network;                                                   ///< Interaction network
    indexer<dim> indexdata;                                             ///< Data structure for indexing
    pairpotentials v;                                                   ///< Pair potential functor
    externalforces<dim> f;                                              ///< External forces functor
    hooker<dim> hooks;                                                  ///< Hook functor
    integrators integrator;                                             ///< Integration method
    variadic_vars<dim> vvars;                                           ///< Bunch of variables for variadic functions
    additional_vars<dim> avars;                                         ///< Bunch of additonal variables
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    md(ui particlenr=0);                                                ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void init(ui particlenr);                                           ///< Copy of the particle number constructor
    ldf dap(ui d,ldf ad);                                               ///< Manipulate particle distances with respect to periodic boundary conditions
    ldf distsq(ui p1,ui p2);                                            ///< Calculate distances between two particles (squared)
    ldf distsq(ldf x1[dim],ldf x2[dim]);                                ///< Calculate distances between two particles (squared)
    ldf distsq(ui p1,ldf x2[dim]);                                      ///< Calculate distances between two particles (squared)
    ldf distsq(ldf x1[dim],ui p2);                                      ///< Calculate distances between two particles (squared)
    ldf dd(ui d,ui p1,ui p2);                                           ///< Calculate difference in particle positions in certain dimension i by particle index
    ldf dd(ui d,ldf x1[dim],ldf x2[dim]);                               ///< Calculate difference in particle positions in certain dimension i by particle index
    ldf dd(ui d,ui p1,ldf x2[dim]);                                     ///< Calculate difference in particle positions in certain dimension i by particle index
    ldf dd(ui d,ldf x1[dim],ui p2);                                     ///< Calculate difference in particle positions in certain dimension i by particle index
    ldf dv(ui d,ui p1,ui p2);                                           ///< Calculate difference in particle velocities in certain dimension i by particle index
    void interactions(ui i,std::vector<std::pair<ui,ui>> &table);       ///< Dump interactions of a certain particle into a table
    void all_interactions(std::vector<std::pair<ui,ui>> &table);        ///< Dump all interaction into a table
    ui add_interaction(ui potential,std::vector<ldf> &parameters);      ///< Add type interaction rule
    ui add_interaction(ui potential,ldf rco,std::vector<ldf> &parameters);///< Add type interaction rule
    bool mod_interaction(ui interaction,ui potential,std::vector<ldf> &parameters);///< Modify type interaction rule
    bool mod_interaction(ui interaction,ui potential,ldf rco,std::vector<ldf> &parameters);///< Modify type interaction rule
    bool rem_interaction(ui interaction);                               ///< Delete type interaction rule
    bool add_typeinteraction(ui type1,ui type2,ui interaction);         ///< Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui interaction);         ///< Modify type interaction rule
    void mad_typeinteraction(ui type1,ui type2,ui interaction);         ///< Force add/mod type interaction rule
    bool add_typeinteraction(ui type1,ui type2,ui potential,std::vector<ldf> &parameters);///< Add type interaction rule
    bool add_typeinteraction(ui type1,ui type2,ui potential,ldf rco,std::vector<ldf> &parameters);///< Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui potential,std::vector<ldf> &parameters);///< Modify type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui potential,ldf rco,std::vector<ldf> &parameters);///< Modify type interaction rule
    void mad_typeinteraction(ui type1,ui type2,ui potential,std::vector<ldf> &parameters);///< Force add/mod type interaction rule
    void mad_typeinteraction(ui type1,ui type2,ui potential,ldf rco,std::vector<ldf> &parameters);///< Force add/mod type interaction rule
    bool rem_typeinteraction(ui type1,ui type2);                        ///< Delete type interaction rule
    ui add_sptype();                                                    ///< Add superparticletype
    bool rem_sptype(ui spt);                                            ///< Delete superparticletype
    bool add_sp_interaction(ui spt,ui p1,ui p2,ui interaction);         ///< Add superparticle interaction rule
    bool mod_sp_interaction(ui spt,ui p1,ui p2,ui interaction);         ///< Modify superparticle interaction rule
    ui mad_sp_interaction(ui spt,ui p1,ui p2,ui interaction);           ///< Force add/mod superparticle interaction rule
    bool add_sp_interaction(ui spt,ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Add superparticle interaction rule
    bool add_sp_interaction(ui spt,ui p1,ui p2,ui potential,ldf rco,std::vector<ldf> &parameters);///< Add superparticle interaction rule
    bool mod_sp_interaction(ui spt,ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Modify superparticle interaction rule
    bool mod_sp_interaction(ui spt,ui p1,ui p2,ui potential,ldf rco,std::vector<ldf> &parameters);///< Modify superparticle interaction rule
    ui mad_sp_interaction(ui spt,ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Force add/mod superparticle interaction rule
    ui mad_sp_interaction(ui spt,ui p1,ui p2,ui potential,ldf rco,std::vector<ldf> &parameters);///< Force add/mod superparticle interaction rule
    bool rem_sp_interaction(ui spt,ui p1,ui p2);                        ///< Delete superparticle interaction rule
    ui add_forcetype(ui force,std::vector<ldf> &parameters);            ///< Add force type
    ui add_forcetype(ui force,std::vector<std::vector<ui>> &plist,std::vector<ldf> &parameters);///< Add force type with plist
    bool mod_forcetype(ui ftype,ui force,std::vector<ldf> &parameters); ///< Modify force type
    bool mod_forcetype(ui ftype,ui force,std::vector<std::vector<ui>> &plist,std::vector<ldf> &parameters);///< Modify force type plist
    bool rem_forcetype(ui ftype);                                       ///< Delete force type
    bool assign_forcetype(ui i,ui ftype);                               ///< Assign force type to particle
    void assign_all_forcetype(ui ftype);                                ///< Assign force type to all particles
    bool unassign_forcetype(ui i,ui ftype);                             ///< Unassign force type to particle
    void unassign_all_forcetype(ui ftype);                              ///< Unassign force type to all particles
    void clear_all_assigned_forcetype();                                ///< Clear all assigned forces
    ui add_hook(ui nohook,std::vector<ldf> &parameters);
    bool mod_hook(ui htype,ui nohook,std::vector<ldf> &parameters);
    bool rm_hook(ui htype);
    bool run_hook(ui htype);
    void run_hooks();
    ldf get_rco(ui i,ui j);                                             ///< Gets the cuttoff radius for a certain pair of particles
    ldf get_rco(ui interaction);                                        ///< Gets the cuttoff radius for a certain interaction
    void set_rco(ldf rco);                                              ///< Sets the cuttoff radius
    void set_rco(ui interaction,ldf rco);                               ///< Sets the cuttoff radius
    void set_ssz(ldf ssz);                                              ///< Sets the skin size radius and its square
    void set_reserve(ldf ssz);                                          ///< Set reserve memory according to skin size
    void set_reserve(ldf ssz,ui M);                                     ///< Set reserve memory according to skin size and some arbitrary number of particles
    void set_type(ui p, ui newtype);                                    ///< Update the type associated with particle p
    void set_index_method(ui method);                                   ///< Set indexmethod
    void index();                                                       ///< Find neighbors
    bool test_index();                                                  ///< Test if we need to run the indexing algorithm
    void thread_index_stick(ui i);                                      ///< Save the particle position at indexing
    ui kdtree_build (ui first, ui last, ui level);                      ///< k-d tree indexing algorithm: tree build function (recursive)
    void kdtree_index (ui first1, ui last1, ui first2, ui last2);       ///< k-d tree indexing algorithm: neighbor finder (recursive)
    void kdtree();                                                      ///< k-d tree indexing algorithm
    void cell();                                                        ///< Cell indexing algorithm
    void thread_cell (ui c);                                            ///< Cell indexer for cell c (thread)
    void bruteforce();                                                  ///< Bruteforce indexing algorithm
    void skinner(ui i,ui j,ldf sszsq);                                  ///< Places interactionneighbor in skin
    void thread_clear_forces(ui i);                                     ///< Clear forces for particle i
    void thread_calc_pot_forces(ui i);                                  ///< Calculate the forces for particle i>j with atomics
    void thread_calc_ext_forces(ui i);                                  ///< Calculate the forces for particle i>j with atomics
    virtual void calc_forces();                                         ///< Calculate the forces between interacting particles
    virtual void recalc_forces();                                       ///< Recalculate the forces between interacting particles for Velocity Verlet
    void update_boundaries();                                           ///< Shifts the periodic boxes appropriately for sheared BC
    void periodicity();                                                 ///< Called after integration to keep the particle within the defined boundaries
    void thread_periodicity(ui i);                                      ///< Apply periodicity to one particle only
    void thread_periodicity(ldf x[dim]);                                ///< Apply periodicity to one particle only
    void thread_seuler(ui i);                                           ///< Symplectic euler integrator (threaded)
    void thread_vverlet_x(ui i);                                        ///< Velocity verlet integrator for position (threaded)
    void thread_vverlet_dx(ui i);                                       ///< Velocity verlet integrator for velocity (threaded)
    void thread_first_order(ui i);                                      ///< First order (Euler) integrator (threaded)
    virtual void integrate();                                           ///< Integrate particle trajectoriess
    void timestep();                                                    ///< Do one timestep
    void timesteps(ui k);                                               ///< Do multiple timesteps
    void import_pos(ldf *x);                                            ///< Load positions from arrays
    template<typename...arg> void import_pos(ldf *x,arg...argv);        ///< Load positions from arrays
    void import_pos(ui i,ldf x);                                        ///< Load position for i from value
    template<typename...arg> void import_pos(ui i,ldf x,arg...argv);    ///< Load position for i from value
    void import_vel(ldf *dx);                                           ///< Load velocity from arrays
    template<typename...arg> void import_vel(ldf *dx,arg...argv);       ///< Load velocity from arrays
    void import_vel(ui i,ldf dx);                                       ///< Load velocity for i from value
    template<typename...arg> void import_vel(ui i,ldf dx,arg...argv);   ///< Load velocity for i from value
    void import_force(ldf *F);                                          ///< Load forces from arrays
    template<typename...arg> void import_force(ldf *F,arg...argv);      ///< Load forces from arrays
    void import_force(ui i,ldf F);                                      ///< Load position for i from value
    template<typename...arg> void import_force(ui i,ldf F,arg...argv);  ///< Load position for i from value
    void export_pos(ldf *x);                                            ///< Save positions from arrays
    template<typename...arg> void export_pos(ldf *x,arg...argv);        ///< Save positions to arrays
    void export_pos(ui i,ldf &x);                                       ///< Save positions from arrays
    template<typename...arg> void export_pos(ui i,ldf &x,arg...argv);   ///< Save positions to arrays
    void export_vel(ldf *dx);                                           ///< Save velocity from arrays
    template<typename...arg> void export_vel(ldf *dx,arg...argv);       ///< Save velocity to arrays
    void export_vel(ui i,ldf &dx);                                      ///< Save positions from arrays
    template<typename...arg> void export_vel(ui i,ldf &dx,arg...argv);  ///< Save positions to arrays
    void export_force(ldf *F);                                          ///< Save forces from arrays
    template<typename...arg> void export_force(ldf *F,arg...argv);      ///< Save forces to arrays
    void export_force(ui i,ldf &F);                                     ///< Save forces from arrays
    template<typename...arg> void export_force(ui i,ldf &F,arg...argv); ///< Save forces to arrays
    ldf direct_readout_x(ui d,ui i);                                    ///< Directly readout a position
    ldf direct_readout_dx(ui d,ui i);                                   ///< Directly readout a velocity
    ldf direct_readout_F(ui d,ui i);                                    ///< Directly readout a forces
    ldf direct_readout(ui d,ui i,uc type);                              ///< Directly readout a position'x'/velocity'v'/forces'F'
    ldf direct_readout(ui i,uc type);                                   ///< Directly readout a position'x'/velocity'v'/forces'F'
    void fix_particle(ui i,bool fix);                                   ///< Fix a particle
    void fix_sp(ui spi,bool fix);                                       ///< Fix a super particles
    ui clone_particle(ui i,ldf x[dim]);                                 ///< Clone a particle and translate
    ui clone_sp(ui spi,ldf x[dim]);                                     ///< Clone a superparticle and translate
    void translate_particle(ui i,ldf x[dim]);                           ///< Translate (or move) a particle
    void translate_sp(ui spi,ldf x[dim]);                               ///< Translate (or move) a super particle
    void drift_particle(ui i,ldf dx[dim]);                              ///< Add velocity to a particle
    void drift_sp(ui spi,ldf dx[dim]);                                  ///< Add velocity to a super particle (all particles the same)
    void heat_particle(ui i,ldf lambda);                                ///< Multiply velocity vector of a particle with a scalar
    void heat_sp(ui spi,ldf lambda);                                    ///< Multiply velocity vectors of a super particle with a scalar (all particles the same)
    void set_position_sp(ui spi,ldf x[dim]);                            ///< Get center of mass of super particle
    void set_velocity_sp(ui spi,ldf dx[dim]);                           ///< Assign velocity to a super particle (all particles the same)
    void get_position_sp(ui spi,ldf x[dim]);                            ///< Get center of mass of super particle
    void get_velocity_sp(ui spi,ldf dx[dim]);                           ///< Get average velocity of a super particle
    ui add_sp(ui sptype);                                               ///< Add a superparticle
    bool rem_sp(ui spi);                                                ///< Remove a superparticle (i.e. the structure, not the particles)
    bool rem_sp_particles(ui spi);                                      ///< Remove all particles in a superparticle
    ui sp_ingest(ui spi,ui i,ui idx=UI_MAX);                            ///< Add a particle to a superparticle
    bool sp_dispose(ui i);                                              ///< Remove a particle from a superparticle
    bool sp_dispose_idx(ui spi,ui idx);                                 ///< Remove a particle from a superparticle
    ui sp_pid(ui spi,ui idx);                                           ///< Reverse lookup of particle id in superparticle
    ui add_particle(ldf mass=1.0,ui ptype=0,bool fixed=false);          ///< Add a particle to the system
    ui add_particle(ldf x[dim],ldf mass=1.0,ui ptype=0,bool fixed=false);///< Add a particle to the system at certain position
    ui add_particle(ldf x[dim],ldf dx[dim],ldf mass=1.0,ui ptype=0,bool fixed=false);///< Add a particle to the system at certain position with certain velocity
    ui add_particle(ldf x[dim],ldf dx[dim],uc bcond[dim],ldf mass=1.0,ui ptype=0,bool fixed=false,bool bconded=true);///< Add a particle to the system at certain position with certain velocity with certain bconds
    void rem_particle(ui i);                                            ///< Remove a particle from the system
    void clear();                                                       ///< Clear all particles and interactions
    void set_damping(ldf coefficient);                                  ///< Enables damping and sets damping coefficient
    bool unset_damping();                                               ///< Disables damping
    void update_skins(ui p1,ui p2);                                     ///< Modify skins after adding/modifying/removing bond
    bool add_bond(ui p1,ui p2,ui interaction);                          ///< Add a bond
    bool mod_bond(ui p1,ui p2,ui interaction);                          ///< Modify a bond
    void mad_bond(ui p1,ui p2,ui interaction);                          ///< Force add/modify bond
    bool add_bond(ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Add a bond
    bool mod_bond(ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Modify a bond
    void mad_bond(ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Force add/modify bond
    bool rem_bond(ui p1,ui p2);                                         ///< Remove a bond from the system
    void assign_unique_types(ui p1, ui p2);                             ///< Assign unique types to particles, modify lookup
    void add_spring(ui p1, ui p2,ldf springconstant,ldf l0);            ///< Add a harmonic bond to the system
    bool add_sp_bond(ui p1,ui p2,ui interaction);                       ///< Add a superparticle bond
    bool mod_sp_bond(ui p1,ui p2,ui interaction);                       ///< Modify a superparticle bond
    void mad_sp_bond(ui p1,ui p2,ui interaction);                       ///< Force add/modify superparticle bond
    bool add_sp_bond(ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Add a superparticle bond
    bool mod_sp_bond(ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Modify a superparticle bond
    void mad_sp_bond(ui p1,ui p2,ui potential,std::vector<ldf> &parameters);///< Force add/modify superparticle bond
    bool rem_sp_bond(ui p1,ui p2);                                      ///< Remove a superparticle bond from the system
    ui clone_sptype(ui sp);                                             ///< Make a new sptype for superparticle sp if it is not unique to sp
    void set_bcond(uc bcond[dim]);                                      ///< Set the global boundary conditions
    void set_pbcond(ui i,uc bcond[dim],bool toggle=true);               ///< Set the boundary conditions for particle i
    void set_spbcond(ui spi,uc bcond[dim],bool toggle=true);            ///< Set the boundary conditions for superparticle spi
    ldf thread_H(ui i);                                                 ///< Measure Hamiltonian for particle i
    virtual ldf thread_T(ui i);                                         ///< Measure kinetic energy for particle i
    virtual ldf thread_V(ui i,bool higher_index_only=false);            ///< Measure potential energy for particle i
    ldf H();                                                            ///< Measure Hamiltonian
    ldf T();                                                            ///< Measure kinetic energy
    ldf V();                                                            ///< Measure potential energy
};

/// Autodiff for Monge patches
template<ui dim> struct duals
{
    ldf x;                                                              ///< Function value
    ldf dx[dim];                                                        ///< First derivatives
    ldf dxdy[dim][dim];                                                 ///< Second derivatives
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    duals();                                                            ///< Constructor
    duals(ldf a);                                                       ///< Constructor
    duals(ldf a,ui i);                                                  ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    duals<dim> operator=(duals<dim> a);                                 ///< Assign operator
    template<class X> duals<dim> operator=(X a);                        ///< Assign foreign type operator
    template<class X> operator X();                                     ///< Cast overload
};

/// This structure defines the Monge patch manifold and its properties
template<ui dim> struct mp
{
    ui patch;                                                           ///< Monge patch type number UI_MAX is custom
    std::vector<ldf> parameters;                                        ///< Monge patch function parameters
    fmpptr<ldf,dim> fmp;                                                ///< Monge patch function
    fmpptr<duals<dim>,dim> dfmp;                                        ///< Derivatives of monge function
    std::vector<duals<dim>> geometryx;                                  ///< Geometric information for particle i in position at position x
    std::vector<duals<dim>> geometryxp;                                 ///< Geometric information for particle i in position at position xp
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mp();                                                               ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void setmp(ui i=MP::FLATSPACE);                                     ///< Picks one of the builtin Monge patches
    void setmp(fmpptr<ldf,dim> f,fmpptr<duals<dim>,dim> df);            ///< Picks a custom Monge patch
    void calc(ui i,ldf x[dim]);                                         ///< Calculate geometric information
    void calc(duals<dim> &z,ldf x[dim]);                                ///< Calculate geometric information on the spot
    ldf f(ldf x[dim]);                                                  ///< Monge patch
    ldf df(ui mu,ldf x[dim]);                                           ///< Monge patch gradient
    ldf ddf(ui mu,ui nu,ldf x[dim]);                                    ///< Monge patch laplacian
    ldf g(ui i,ui mu,ui nu);                                            ///< Monge patch metric tensor
    ldf gp(ui i,ui mu,ui nu);                                           ///< Monge patch metric tensor
    ldf ginv(ui i,ui mu,ui nu);                                         ///< Monge patch metric tensor inverse
    ldf A(ui i,ui sigma,ui mu,ui nu);                                   ///< Monge patch \f$ A_{\sigma \mu \nu} = \Gamma_{\nu \sigma \mu} \f$ where \f$ \Gamma_{\nu \sigma \mu} \f$ are the Christoffel symbols (of the first kind)
};

/// This structure takes care of Monge patch molecular dynamics simulations
template<ui dim> struct mpmd:md<dim>
{
    mp<dim> patch;                                                      ///< Geometric monge patch information
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    using md<dim>::md;
    using md<dim>::N;
    using md<dim>::simbox;
    using md<dim>::particles;
    using md<dim>::network;
    using md<dim>::indexdata;
    using md<dim>::v;
    using md<dim>::f;
    using md<dim>::integrator;
    using md<dim>::avars;
    using md<dim>::thread_clear_forces;
    using md<dim>::thread_calc_ext_forces;
    using md<dim>::periodicity;
    using md<dim>::thread_seuler;
    using md<dim>::thread_vverlet_x;
    using md<dim>::thread_vverlet_dx;
    using md<dim>::distsq;
    using md<dim>::dd;
    using md<dim>::dap;
    using md<dim>::index;
    using md<dim>::test_index;
    using md<dim>::get_rco;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf embedded_distsq(ui p1,ui p2);                                   ///< Calculate distances between two particles (squared)
    ldf embedded_distsq(ldf x1[dim],ldf x2[dim]);                       ///< Calculate distances between two particles (squared)
    ldf embedded_distsq(ui p1,ldf x2[dim]);                             ///< Calculate distances between two particles (squared)
    ldf embedded_distsq(ldf x2[dim],ui p2);                             ///< Calculate distances between two particles (squared)
    ldf embedded_dd_p1(ui d,ui p1,ui p2);                               ///< Calculate particles relative particle in certain dimension i wrt p1
    ldf embedded_dd_p2(ui d,ui p1,ui p2);                               ///< Calculate particles relative particle in certain dimension i wrt p2
    void zuiden_C(ui i,ldf ZC[dim]);                                    ///< Calculates \f$g^{\rho \sigma} C_{\sigma}\f$ for particle i of the van Zuiden integrator
    void zuiden_A(ui i,ldf eps[dim]);                                   ///< Calculates \f$g^{\rho \sigma} A_{\sigma \mu \nu} \epsilon^{\mu} \epsilon^{\nu}\f$ for particle i of the van Zuiden integrator
    void thread_zuiden_wfi(ui i);                                       ///< The van Zuiden integrator without fixed point itterations
    void thread_zuiden_protect(ui i);                                   ///< The van Zuiden integrator with protected fixed point itterations (makes sure you don't get stuck in a loop)
    void thread_zuiden(ui i);                                           ///< The van Zuiden integrator for Riemannian manifolds (fails for pseudo-Riemannian manifolds)
    void thread_history(ui i);                                          ///< Set the history of particle i
    void history();                                                     ///< Set the history of all particles
    void thread_calc_geometry(ui i);                                    ///< Calculate Monge patch derivatives for partice i
    void calc_geometry();                                               ///< Calculate Monge patch derivatives
    void mp_thread_calc_pot_forces(ui i);                               ///< Calculate the forces for particle i>j with atomics
    void integrate() override final;                                    ///< Integrate particle trajectoriess
    void calc_forces() override final;                                  ///< Integrate particle trajectoriess
    void recalc_forces() override final;                                ///< Integrate particle trajectoriess
    ldf thread_T(ui i) override final;                                  ///< Calculate kinetic energy of a particle
    ldf thread_V(ui i,bool higher_index_only=false) override final;     ///< Calculate potential energy
};

#ifndef __libmd_cc__
#ifndef __libmd_src_file__
#include "libmd.cc"
#endif
#endif

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of libmd HEADER file                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

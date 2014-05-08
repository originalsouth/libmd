///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of libmd HEADER file                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef libmd_h
#define libmd_h

#include <cstdio>                                                       //Standard input output (faster than IOstream and also threadsafe) (C)
#include <cstdlib>                                                      //Standard library (C)
#include <cmath>                                                        //Standard math library  (C)
#include <cstring>                                                      //Memcpy and memmove support (C)
#include <vector>                                                       //Vector support (C++)
#include <map>                                                          //Map support (C++)
#include <unordered_map>                                                //Map support (C++)
#include <list>                                                         //List support (C++)
#include <stack>                                                        //Stack support (C++)
#include <set>                                                          //Set support (C++)
#include <unordered_set>                                                //Unordered set support (C++11)
#include <utility>                                                      //Pair support (C++)
#include <limits>                                                       //Limits of types (C++)
#include <thread>                                                       //Thread support (C++11)
#include <mutex>                                                        //Mutex support (C++11)
#include <future>                                                       //Future support (C++11)
#include <algorithm>                                                    //Algorithm support (C++)
#include <functional>

#ifdef FE
#include <fenv.h>                                                        //Floating point exception handling (C)
#endif

#ifdef TIMER
#include <sys/time.h>
long double TicToc()
{
    static timeval start,end;
    start=end,gettimeofday(&end,NULL);
    return (((long double)(end.tv_sec-start.tv_sec))+((long double)(end.tv_usec-start.tv_usec))/1000000.0);
}
#endif

using namespace std;                                                    //Using standard namespace

typedef long double ldf;                                                //long double is now aliased as ldf
typedef unsigned int ui;                                                //unsigned int is now aliased as ui
typedef unsigned char uc;                                               //unsigned char is now aliased as uc

struct INTEGRATOR {enum intergrator:uc {SEULER,VVERLET};};                      ///< Integration options
struct MP_INTEGRATOR {enum mp_integrator:uc {VZ,VZ_P,VZ_WFI,SEULER,VVERLET};};  ///< Monge patch integration options
struct BCOND {enum bcond:uc {NONE,PERIODIC,HARD,BOXSHEAR};};                    ///< Boundary condition options
struct INDEX {enum index:uc {CELL,BRUTE_FORCE,KD_TREE};};                       ///< Indexing options
struct POT {enum pot:ui                                                         ///< Potential options
{
    COULOMB,
    YUKAWA,
    HOOKEAN,
    LJ,
    MORSE,
    FORCEDIPOLE,
    HOOKEANFORCEDIPOLE,
    ANHARMONICSPRING
};};
/// External force options
struct EXTFORCE {enum extforce:ui
{
    DAMPING,
    DISSIPATION
};};
/// Monge patch options
struct MP {enum mp:ui
{
    FLATSPACE,
    GAUSSIANBUMP,
    EGGCARTON,
    MOLLIFIER
};};

//These functions defined outside of libmd
void __libmd__info();                                                   ///< Basic libmd comilation info

/// This structure takes care of multithreading
struct threads
{
    ui nothreads;                                                       ///< Number of threads
    mutex lock;                                                         ///< Thread blocker
    vector<thread> block;                                               ///< Block of threads
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    threads(ui nrthreads=thread::hardware_concurrency());               ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void set(ui nrthreads=thread::hardware_concurrency());              ///< Set the number of threads (default is max)
};

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
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    particle(ldf mass=1.0,ui ptype=0,bool fixed=false);                 ///< Constructor
};

/// This structure contains information about the simulation box
template<ui dim> struct box
{
    ldf L[dim];                                                         ///< Box size
    bool boxShear;                                                      ///< Use sheared box matrix
    ldf vshear[dim][dim];                                               ///< Shear velocity vshear[i][j] is shear velocity in direction i of boundary with normal in direction j. currently vshear[i][i] != 0 results in undefined behaviour.
    ldf Lshear[dim][dim];                                               ///< Box matrix that is updated at each time step. Used to compute distances for shear, in lieu of simbox.L
    ldf LshearInv[dim][dim];                                            ///< Inverse of Lshear[][]
    uc bcond[dim];                                                      ///< Boundary conditions in different dimensions NONE/PERIODIC/HARD(/boxShear)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    box();                                                              ///< Constructor
    void shear_boundary(ui i, ui j, ldf velocity);                      ///< Set up boundary shear velocity in direction i of boundary with normal direction j
    void skew_boundary(ui i, ui j, ldf displacement);                   ///< Skew the simulation box by moving boundary with normal direction j by amount 'displacement' in direction i
    void invert_box();                                                  ///< Invert the Lshear[][] box matrix 
};

/// This structure saves the particle type interactions and calculates the the potentials
struct interactiontype
{
    ui potential;                                                       ///< Type of potential
    vector<ldf> parameters;                                             ///< Parameters of potential
    ldf vco;                                                            ///< Cuttoff potential energy
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interactiontype(ui ppot,vector<ldf> *param,ldf Vco);                ///< Constructor
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
    vector<vector<ui>> particles;                                       ///< Interacting particle list
    vector<ldf> parameters;                                             ///< Parameters for the external force
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    forcetype(ui noexternalforce,vector<vector<ui>> *plist,vector<ldf> *param); ///< Constructor
};

/// This structure introduces "super_particles" i.e. particles that consist of (sub_)particles
struct superparticle
{
    unordered_map<ui,ui> particles;                                     ///< Particles in super particles
    ui sptype;                                                          ///< Super particle type
};

/// This structure caries a lookup device for a specific super particle type
struct superparticletype
{
    map<pair<ui,ui>,ui> splookup;                                       ///< This is the interaction lookup device
};

/// This structure stores all interactions and their types
struct interact
{
    bool update;                                                        ///< Should we update the network
    ldf rco;                                                            ///< R_cutoff radius
    ldf ssz;                                                            ///< Skin radius
    vector<vector<ui>> forces;                                          ///< List of external forces acting on the particles
    vector<forcetype> forcelibrary;                                     ///< Library of external forces
    vector<vector<interactionneighbor>> skins;                          ///< Particle skin by index (array of vector)
    vector<interactiontype> library;                                    ///< This is the interaction library
    unordered_set<ui> free_library_slots;                               ///< Stores free library slots
    map<pair<ui,ui>,ui> lookup;                                         ///< This is the interaction lookup device
    vector<ui> spid;                                                    ///< Super particle identifier array
    vector<superparticle> superparticles;                               ///< Actual super particle array
    vector<superparticletype> sptypes;                                  ///< Super particle type array
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interact();                                                         ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pair<ui,ui> hash(ui type1,ui type2);                                ///< Hash function
    bool probe(ui type1,ui type2);                                      ///< Check if a typeinteraction exists between two types
};

/// Potential functions
template<class X> X COULOMB(X r,vector<ldf> *parameters);
template<class X> X YUKAWA(X r,vector<ldf> *parameters);
template<class X> X HOOKEAN(X r,vector<ldf> *parameters);
template<class X> X MORSE(X r,vector<ldf> *parameters);
template<class X> X FORCEDIPOLE(X r,vector<ldf> *parameters);
template<class X> X HOOKEANFORCEDIPOLE(X r,vector<ldf> *parameters);
template<class X> X ANHARMONICSPRING(X r,vector<ldf> *parameters);

/// External force functions
template<ui dim> void DAMPING(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters,void *sys);
template<ui dim> void DISSIPATION(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters,void *sys);

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

template<class X> using potentialptr=X (*)(X,vector<ldf> *);            ///< Function pointer to potential functions is now called potentialptr

/// This structure takes care of pair potentials (who live outside of the class)
struct pairpotentials
{
    vector<potentialptr<dual>> potentials;                              ///< Pair potential vector
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pairpotentials();                                                   ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(potentialptr<dual> p);                                       ///< Add a potentials
    ldf operator()(ui type,ldf r,vector<ldf> *parameters);              ///< Pair potential executer
    ldf dr(ui type,ldf r,vector<ldf> *parameters);                      ///< Pair potential d/dr executer
};

template<ui dim> using extforceptr=void (*)(ui,vector<ui> *,vector<ldf> *,void *); ///< Function pointer to external force functions is now called extforceptr

/// This structure takes care of additional (external) forces acting on particles
template<ui dim> struct externalforces
{
    vector<extforceptr<dim>> extforces;                                 ///< External forces function container
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    externalforces();                                                   ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(extforceptr<dim> p);                                         ///< Add an external force function
    void operator()(ui type,ui i,vector<ui> *particles,vector<ldf> *parameters,void *sys); ///< Execute external force function
};

/// This structure defines and saves integration metadata
struct integrators
{
    ldf h;                                                              ///< Timestep size
    uc method;                                                          ///< Type of integration (SEULER/VVERLET)
    uc generations;                                                     ///< Maximum generations of timestep
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
        vector<ui> *Cells;                                              ///< List of particles per cell
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

/// This structure stores some cyclic variables for the variadic functions
template<ui dim> struct variadic_vars
{
    vector<ui> vvars;                                                   ///< Container of variables
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    variadic_vars();                                                    ///< Initialize variables (set everyting to zero)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui operator[](ui i);                                                ///< Rotate and return previous for the ith variable
};

/// This structure stores additional variables
template<ui dim> struct additional_vars
{
    ui noftypedamping=numeric_limits<ui>::max();                        ///< This variable stores the number of the ftype that damps/drags the system
    bool export_force_calc=false;                                       ///< This variable tells export_force if the forces have been calculated for this output
    bool reindex=true;                                                  ///< This variable tells if the system needs to be reindexed
};

/// This structure defines the molecular dynamics simulation
template<ui dim> struct md
{
    ui N;                                                               ///< Number of particles
    box<dim> simbox;                                                    ///< Simulation box
    vector<particle<dim>> particles;                                    ///< Particle array
    interact network;                                                   ///< Interaction network
    indexer<dim> indexdata;                                             ///< Data structure for indexing
    pairpotentials v;                                                   ///< Pair potential functor
    externalforces<dim> f;                                              ///< External forces functor
    integrators integrator;                                             ///< Integration method
    threads parallel;                                                   ///< Multithreader
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
    void all_interactions(vector<pair<ui,ui>> &table);                  ///< Dump all interaction into a table
    ui add_interaction(ui potential,vector<ldf> *parameters);           ///< Add type interaction rule
    bool mod_interaction(ui interaction,ui potential,vector<ldf> *parameters);///< Modify type interaction rule
    bool rem_interaction(ui interaction);                               ///< Delete type interaction rule
    bool add_typeinteraction(ui type1,ui type2,ui interaction);         ///< Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui interaction);         ///< Modify type interaction rule
    void mad_typeinteraction(ui type1,ui type2,ui interaction);         ///< Force add/mod type interaction rule
    bool add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);///< Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);///< Modify type interaction rule
    void mad_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);///< Force add/mod type interaction rule
    bool rem_typeinteraction(ui type1,ui type2);                        ///< Delete type interaction rule
    ui add_sptype();                                                    ///< Add superparticletype
    bool rem_sptype(ui spt);                                            ///< Delete superparticletype
    bool add_sp_interaction(ui spt,ui p1,ui p2,ui interaction);         ///< Add superparticle interaction rule
    bool mod_sp_interaction(ui spt,ui p1,ui p2,ui interaction);         ///< Modify superparticle interaction rule
    ui mad_sp_interaction(ui spt,ui p1,ui p2,ui interaction);           ///< Force add/mod superparticle interaction rule
    bool add_sp_interaction(ui spt,ui p1,ui p2,ui potential,vector<ldf> *parameters);///< Add superparticle interaction rule
    bool mod_sp_interaction(ui spt,ui p1,ui p2,ui potential,vector<ldf> *parameters);///< Modify superparticle interaction rule
    ui mad_sp_interaction(ui spt,ui p1,ui p2,ui potential,vector<ldf> *parameters);///< Force add/mod superparticle interaction rule
    bool rem_sp_interaction(ui spt,ui p1,ui p2);                        ///< Delete superparticle interaction rule
    ui add_forcetype(ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters);///< Add force type
    bool mod_forcetype(ui ftype,ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters);///< Modify force type
    bool rem_forcetype(ui ftype);                                       ///< Delete force type
    bool assign_forcetype(ui i,ui ftype);                               ///< Assign force type to particle
    void assign_all_forcetype(ui ftype);                                ///< Assign force type to all particles
    void unassign_forcetype(ui i,ui ftype);                             ///< Unassign force type to particle
    void unassign_all_forcetype(ui ftype);                              ///< Unassign force type to all particles
    void clear_all_assigned_forcetype();                                ///< Clear all assigned forces
    void set_rco(ldf rco);                                              ///< Sets the cuttoff radius and its square
    void set_ssz(ldf ssz);                                              ///< Sets the skin size radius and its square
    void set_reserve(ldf ssz);                                          ///< Set reserve memory according to skin size
    void set_reserve(ldf ssz,ui M);                                     ///< Set reserve memory according to skin size and some arbitrary number of particles
    void set_type(ui p, ui newtype);                                    ///< Update the type associated with particle p
    void set_index_method(ui method);                                   ///< Set indexmethod
    void thread_index(ui i);                                            ///< Find neighbors per cell i
    void index();                                                       ///< Find neighbors
    bool test_index();                                                  ///< Test if we need to run the indexing algorithm
    void thread_index_stick(ui i);                                      ///< Save the particle position at indexing
    ui kdtree_build (ui first, ui last, ui level);                      ///< k-d tree indexing algorithm: tree build function (recursive)
    void kdtree_index (ui first1, ui last1, ui first2, ui last2);       ///< k-d tree indexing algorithm: neighbor finder (recursive)
    void kdtree();                                                      ///< k-d tree indexing algorithm
    void cell();                                                        ///< Cell indexing algorithm
    void thread_cell (ui i);                                            ///< Cell indexer for cell i (thread)
    void bruteforce();                                                  ///< Bruteforce indexing algorithm
    void skinner(ui i,ui j);                                            ///< Places interactionneighbor in skin
    void thread_clear_forces(ui i);                                     ///< Clear forces for particle i
    void thread_calc_forces(ui i);                                      ///< Calculate the forces for particle i>j with atomics
    virtual void calc_forces();                                         ///< Calculate the forces between interacting particles
    virtual void recalc_forces();                                       ///< Recalculate the forces between interacting particles for Velocity Verlet
    void update_boundaries();                                           ///< Shifts the periodic boxes appropriately for sheared BC
    void periodicity();                                                 ///< Called after integration to keep the particle within the defined boundaries
    void thread_periodicity(ui i);                                      ///< Apply periodicity to one particle only
    void thread_periodicity_periodic(ui d,ui i);                        ///< Called by periodicity to keep periodic boundary conditions
    void thread_periodicity_boxshear(ui d,ui i);                        ///< Called by periodicity to keep boxshear boundary conditions
    void thread_periodicity_hard(ui d,ui i);                            ///< Called by periodicity to keep hard boundary conditions
    void thread_seuler(ui i);                                           ///< Symplectic euler integrator (threaded)
    void thread_vverlet_x(ui i);                                        ///< Velocity verlet integrator for position (threaded)
    void thread_vverlet_dx(ui i);                                       ///< Velocity verlet integrator for velocity (threaded)
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
    ui sp_ingest(ui spi,ui i);                                          ///< Add a particle to a superparticle
    bool sp_dispose(ui i);                                              ///< Remove a particle from a superparticle
    ui sp_pid(ui spi,ui idx);                                           ///< Reverse lookup of particle id in superparticle
    ui add_particle(ldf mass=1.0,ui ptype=0,bool fixed=false);          ///< Add a particle to the system
    ui add_particle(ldf x[dim],ldf mass=1.0,ui ptype=0,bool fixed=false);///< Add a particle to the system at certain position
    ui add_particle(ldf x[dim],ldf dx[dim],ldf mass=1.0,ui ptype=0,bool fixed=false);///< Add a particle to the system at certain position with certain velocity
    void rem_particle(ui i);                                            ///< Remove a particle from the system
    void clear();                                                       ///< Clear all particles and interactions
    void set_damping(ldf coefficient);                                  ///< Enables damping and sets damping coefficient
    bool unset_damping();                                               ///< Disables damping
    void uitopptr(vector<particle<dim>*> *x,vector<ui> i);              ///< Convert vector of unsigned integers to particle pointers
    ui pptrtoui(particle<dim> *x);                                      ///< Convert a particle pointer to a particle id
    void update_skins(ui p1,ui p2);                                     ///< Modify skins after adding/modifying/removing bond
    bool add_bond(ui p1,ui p2,ui interaction);                          ///< Add a bond
    bool mod_bond(ui p1,ui p2,ui interaction);                          ///< Modify a bond
    void mad_bond(ui p1,ui p2,ui interaction);                          ///< Force add/modify bond
    bool add_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters);    ///< Add a bond
    bool mod_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters);    ///< Modify a bond
    void mad_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters);    ///< Force add/modify bond
    bool rem_bond(ui p1,ui p2);                                         ///< Remove a bond from the system
    void assign_unique_types(ui p1, ui p2);                             ///< Assign unique types to particles, modify lookup
    void add_spring(ui p1, ui p2,ldf springconstant,ldf l0);            ///< Add a harmonic bond to the system
    bool add_sp_bond(ui p1,ui p2,ui interaction);                       ///< Add a superparticle bond
    bool mod_sp_bond(ui p1,ui p2,ui interaction);                       ///< Modify a superparticle bond
    void mad_sp_bond(ui p1,ui p2,ui interaction);                       ///< Force add/modify superparticle bond
    bool add_sp_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters); ///< Add a superparticle bond
    bool mod_sp_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters); ///< Modify a superparticle bond
    void mad_sp_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters); ///< Force add/modify superparticle bond
    bool rem_sp_bond(ui p1,ui p2);                                      ///< Remove a superparticle bond from the system
    ui clone_sptype(ui sp);                                             ///< Make a new sptype for superparticle sp if it is not unique to sp
    ldf thread_H(ui i);                                                 ///< Measure Hamiltonian for particle i
    virtual ldf thread_T(ui i);                                         ///< Measure kinetic energy for particle i
    virtual ldf thread_V(ui i);                                         ///< Measure potential energy for particle i
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

/// Monge patch function pointer
template<class X,ui dim> using fmpptr=X (*)(X x[dim],vector<ldf> *param);

/// Monge patches (and related)
ldf kdelta(ui i,ui j);
template<class X,ui dim> X FLATSPACE(X x[dim],vector<ldf> *param);
template<class X,ui dim> X GAUSSIANBUMP(X x[dim],vector<ldf> *param);

/// This structure defines the Monge patch manifold and its properties
template<ui dim> struct mp
{
    ui patch;                                                           ///< Monge patch type number numeric_limits<ui>::max() is custom
    vector<ldf> parameters;                                             ///< Monge patch function parameters
    fmpptr<ldf,dim> fmp;                                                ///< Monge patch function
    fmpptr<duals<dim>,dim> dfmp;                                        ///< Derivatives of monge function
    vector<duals<dim>> geometryx;                                       ///< Geometric information for particle i in position at position x
    vector<duals<dim>> geometryxp;                                      ///< Geometric information for particle i in position at position xp
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mp();                                                               ///< Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void setmp(ui i=MP::GAUSSIANBUMP);                                  ///< Picks one of the builtin Monge patches
    void setmp(fmpptr<ldf,dim> f,fmpptr<duals<dim>,dim> df);            ///< Picks a custom Monge patch
    void calc(ui i,ldf x[dim]);                                         ///< Calculate geometric information
    void calc(duals<dim> &z,ldf x[dim]);                                ///< Calculate geometric information on the spot
    ldf f(ldf x[dim]);                                                  ///< Monge patch
    ldf df(ui mu,ldf x[dim]);                                           ///< Monge patch gradient
    ldf ddf(ui mu,ui nu,ldf x[dim]);                                    ///< Monge patch laplacian
    ldf g(ui i,ui mu,ui nu);                                            ///< Monge patch metric tensor
    ldf gp(ui i,ui mu,ui nu);                                           ///< Monge patch metric tensor
    ldf ginv(ui i,ui mu,ui nu);                                         ///< Monge patch metric tensor inverse
    ldf A(ui i,ui sigma,ui mu,ui nu);                                   ///< Monge patch \f$ A_{\sigma \mu \nu} = \Gamma_{\nu \sigma \mu} \f$ where \f$ Gamma_{\nu \sigma \mu} \f$ are the Christoffel symbols (of first kind)
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
    using md<dim>::parallel;
    using md<dim>::periodicity;
    using md<dim>::thread_seuler;
    using md<dim>::thread_vverlet_x;
    using md<dim>::thread_vverlet_dx;
    using md<dim>::distsq;
    using md<dim>::dd;
    using md<dim>::dap;
    using md<dim>::index;
    using md<dim>::test_index;
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
    void mp_thread_calc_forces(ui i);                                   ///< Calculate the forces for particle i>j with atomics
    void integrate() override;                                          ///< Integrate particle trajectoriess
    void calc_forces() override;                                        ///< Integrate particle trajectoriess
    void recalc_forces() override;                                      ///< Integrate particle trajectoriess
    ldf thread_T(ui i) override;                                        ///< Calculate kinetic energy of a particle
    ldf thread_V(ui i) override;                                        ///< Calculate kinetic energy
};

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of libmd HEADER file                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

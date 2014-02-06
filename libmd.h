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
#include <list>                                                         //List support (C++)
#include <set>                                                          //Set support (C++)
#include <utility>                                                      //Pair support (C++)
#include <limits>                                                       //Limits of types (C++)
#include <thread>                                                       //Thread support (C++11)
#include <mutex>                                                        //Mutex support (C++11)
#include <future>                                                       //Future support (C++11)

using namespace std;                                                    //Using standard namespace

typedef long double ldf;                                                //long double is now aliased as ldf
typedef unsigned int ui;                                                //unsigned int is now aliased as ui
typedef unsigned char uc;                                               //unsigned int is now aliased as uc

typedef ldf (*fmpptr)(ldf *,vector<ldf> *);                             //Monge patch function pointer
typedef ldf (*dfmpptr)(ui,ldf *,vector<ldf> *);                         //Monge patch function derivative pointer
typedef ldf (*ddfmpptr)(ui,ui,ldf *,vector<ldf> *);                     //Monge patch function second derivative pointer

enum INTEGRATOR:uc {SEULER,VVERLET};                                    //Integration options
enum MP_INTEGRATOR:uc {MP_VZ,MP_VZ_P,MP_VZ_WFI,MP_SEULER,MP_VVERLET};   //Monge patch integration options
enum BCOND:uc {NONE,PERIODIC,HARD,BOXSHEAR};                            //Boundary condition options
enum INDEX:uc {CELL,BRUTE_FORCE};                                       //Indexing options
enum POT:ui                                                             //Potential options
{
    POT_COULOMB,
    POT_YUKAWA,
    POT_HOOKIAN,
    POT_LJ,
    POT_MORSE,
    POT_FORCEDIPOLE,
    POT_HOOKEANFORCEDIPOLE,
    POT_ANHARMONICSPRING
};
enum EXTFORCE:ui                                                        //External force options
{
    EXTFORCE_DAMPING
};
enum MP:ui                                                              //Monge patch options
{
    MP_FLATSPACE,
    MP_GAUSSIANBUMP
};

//This structure takes care of multithreading
struct threads
{
    ui nothreads;                                                       //Number of threads
    mutex lock;                                                         //Thread blocker
    vector<thread> block;                                               //Block of threads
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    threads(ui nrthreads=thread::hardware_concurrency());               //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void set(ui nrthreads=thread::hardware_concurrency());              //Set the number of threads (default is max)
};

//This structure contains all the information for a single particle
template<ui dim> struct particle                 
{
    ldf m;                                                              //Mass
    ldf x[dim];                                                         //Position
    ldf xp[dim];                                                        //Previous particle position
    ldf xsk[dim];                                                       //Position when building the skinlist
    ldf dx[dim];                                                        //Velocity
    ldf F[dim];                                                         //Forces on particle
    ui type;                                                            //This particle has type number
    bool fix;                                                           //Can this particle move
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    particle(ldf mass=1.0,ui ptype=0,bool fixed=false);                 //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    particle *address();                                                //Return the pointer of the particle
};

//This structure contains information about the simulation box
//TODO: Deformations (talk to Jayson)
//TODO: Wall (talk to Jayson)
template<ui dim> struct box
{
    ldf L[dim];                                                         //Box size
    bool boxShear;                                                      //Use sheared box matrix
    ldf vshear[dim][dim];                                               //Shear velocity vshear[i][j] is shear velocity in direction i of boundary with normal in direction j. currently vshear[i][i] != 0 results in undefined behaviour.
    ldf Lshear[dim][dim];                                               //Box matrix that is updated at each time step. Used to compute distances for shear, in lieu of simbox.L
    ldf LshearInv[dim][dim];                                            //Inverse of Lshear[][]
    uc bcond[dim];                                                      //Boundary conditions in different dimensions NONE/PERIODIC/HARD(/boxShear)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    box();                                                              //Constructor
    void shear_boundary(ui i, ui j, ldf velocity);                      //set up boundary shear velocity in direction i of boundary with normal direction j
    void invert_box();                                                  //invert the Lshear[][] box matrix 
};

//This structure saves the particle type interactions and calculates the the potentials
struct interactiontype
{
    ui potential;                                                       //Type of potential
    vector<ldf> parameters;                                             //Parameters of potential
    ldf vco;                                                            //Cuttoff potential energy
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interactiontype(ui ppot,vector<ldf> *param,ldf Vco);                //Constructor
};

//This struct saves the neighboring particle number and the interaction type library number
struct interactionneighbor
{
    ui neighbor;                                                        //Neighbor number
    ui interaction;                                                     //Interaction number
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interactionneighbor(ui noneighbor,ui nointeraction);                //Constructor
};

struct forcetype
{
    ui externalforce;                                                   //External force type
    vector<vector<ui>> particles;                                       //Interacting particle list
    vector<ldf> parameters;                                             //Parameters for the external force
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    forcetype(ui noexternalforce,vector<vector<ui>> *plist,vector<ldf> *param); //Constructor
};

//This structure stores all interactions and their types
struct interact
{
    bool update;                                                        //Should we update the network
    ldf rco;                                                            //R_cuttoff radius
    ldf rcosq;                                                          //R_cuttoff radius squared
    ldf ssz;                                                            //Skin radius
    ldf sszsq;                                                          //Skin radius squared
    vector<vector<ui>> forces;                                          //List of external forces acting on the particles
    vector<forcetype> forcelibrary;                                     //Library of external forces
    vector<vector<interactionneighbor>> skins;                          //Particle skin by index (array of vector)
    vector<interactiontype> library;                                    //This is the interaction library
    vector<pair<ui,ui>> backdoor;                                       //Inverse lookup device
    map<pair<ui,ui>,ui> lookup;                                         //This is the interaction lookup device
    map<ui,set<ui>> usedtypes;                                          //Map of all used types to points having that type NOTE: no guarantee that this is complete, since user can set particle types without setting this function accordingly!! can change by requiring a set_type() function. TODO
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pair<ui,ui> hash(ui type1,ui type2);                                //Hash function
    bool probe(ui type1,ui type2);                                      //Check if a typeinteraction exists between two types
};

//This structure automatically differentiates first order
struct dual
{
    ldf x;
    ldf dx;
    dual();
    dual(ldf f,ldf fx=1.0);
    dual operator=(dual y);
    void operator+=(dual y);
    void operator-=(dual y);
    template<class X> X operator=(X y);
    template<class X> void operator+=(X y);
    template<class X> void operator-=(X y);
    template<class X> void operator*=(X y);
    template<class X> void operator/=(X y);
    template<class X> bool operator==(X y);
    template<class X> bool operator<=(X y);
    template<class X> bool operator>=(X y);
    template<class X> bool operator<(X y);
    template<class X> bool operator>(X y);
};

typedef dual (*potentialptr)(dual,vector<ldf> *);                       //Function pointer to potential functions is now called potential

//This structure takes care of pair potentials (who live outside of the class)
struct pairpotentials
{
    vector<potentialptr> potentials;                                    //Pair potential vector
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pairpotentials();                                                   //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(potentialptr p);                                             //Add a potentials
    ldf operator()(ui type,ldf r,vector<ldf> *parameters);              //Pair potential executer
    ldf dr(ui type,ldf r,vector<ldf> *parameters);                      //Pair potential d/dr executer
};

template<ui dim> using extforceptr=void (*)(particle<dim> *,vector<particle<dim>*> *,vector<ldf> *);

//This structure takes care of additional (external) forces acting on particles
template<ui dim> struct externalforces
{
    vector<extforceptr<dim>> extforces;                                 //External forces function container
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    externalforces();                                                   //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(extforceptr<dim> p);                                         //Add an external force function
    void operator()(ui type,particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters); //Execute external force function
};

//This structure defines and saves integration metadata
struct integrators
{
    ldf h;                                                              //Timestep size
    uc method;                                                          //Type of integration (SEULER/VVERLET)
    uc generations;                                                     //Maximum generations of timestep
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    integrators();                                                      //Constructor
};

//This structure is specific for the indexer
template<ui dim> struct indexer
{
    uc method;                                                          //Method of indexing
    struct celldatatype
    {
        ui Q[dim];                                                      //Not commented
        ui nCells;                                                      //Total number of cells (= prod(Q))
        ui totNeighbors;                                                //Total number of (potential) neighboring cells to check (= (3^d-1)/2)
        ldf CellSize[dim];                                              //Length of cell in each dimension
        int (*IndexDelta)[dim];                                         //Not commented
        vector<list<ui>> Cells;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        celldatatype();                                                 //Constructor
        ~celldatatype();                                                //Destructor
    };
    celldatatype celldata;                                              //Cell data object
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    indexer();                                                          //Constructor
};

//This structure stores some cyclic variables for the variadic functions
template<ui dim> struct variadic_vars
{
    vector<ui> vvars;                                                   //Container of variables
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    variadic_vars();                                                    //Initialize variables (set everyting to zero)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui operator[](ui i);                                                //Rotate and return previous for the ith variable
};

//This structure stores additional variables
template<ui dim> struct additional_vars
{
    ui noftypedamping;                                                  //This variable stores the number of the ftype that damps/drags the system
    bool export_force_calc;                                             //This variable tells export_force if the forces have been calculated for this output
};

//This structure defines the molecular dynamics simulation
template<ui dim> struct md
{
    ui N;                                                               //Number of particles
    box<dim> simbox;                                                    //Simulation box
    vector<particle<dim>> particles;                                    //Particle array
    interact network;                                                   //Interaction network
    indexer<dim> indexdata;                                             //Data structure for indexing
    pairpotentials v;                                                   //Pair potential functor
    externalforces<dim> f;                                              //External forces functor
    integrators integrator;                                             //Integration method
    threads parallel;                                                   //Multithreader
    variadic_vars<dim> vvars;                                           //Bunch of variables for variadic functions
    additional_vars<dim> avars;                                         //Bunch of additonal variables
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    md();                                                               //Constructor
    md(ui particlenr);                                                  //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf dap(ui i,ldf ad);                                               //Manipulate particle distances with respect to periodic boundary conditions
    ldf distsq(ui p1,ui p2);                                            //Calculate distances between two particles (squared)
    ldf dd(ui i,ui p1,ui p2);                                           //Caculate particles relative particle in certain dimension i
    bool add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Modify type interaction rule
    bool rem_typeinteraction(ui type1,ui type2);                        //Delete type interaction rule
    ui add_forcetype(ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters);    //Add force type
    bool mod_forcetype(ui notype,ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters);    //Modify force type
    bool rem_forcetype(ui notype);                                      //Delete force type
    void assign_forcetype(ui particlenr,ui ftype);                      //Assign force type to particle
    void assign_all_forcetype(ui ftype);                                //Assign force type to all particles
    void unassign_forcetype(ui particlenr,ui ftype);                    //Unassign force type to particle
    void unassign_all_forcetype(ui ftype);                              //Unassign force type to all particles
    void clear_all_assigned_forcetype();                                //Clear all assigned forces
    void set_rco(ldf rco);                                              //Sets the cuttoff radius and its square
    void set_ssz(ldf ssz);                                              //Sets the skin size radius and its square
    void set_type(ui p, ui newtype);                                    //Update the type associated with particle p
    void thread_index(ui i);                                            //Find neighbors per cell i
    void index();                                                       //Find neighbors
    bool test_index();                                                  //Test if we need to run the indexing algorithm
    void thread_index_stick(ui i);                                      //Save the particle position at indexing
    void cell();                                                        //Cell indexing algorithm
    void thread_cell (ui i);                                            //Cell indexer for cell i (thread)
    void bruteforce();                                                  //Bruteforce indexing algorithm
    void thread_clear_forces(ui i);                                     //Clear forces for particle i
    virtual void thread_calc_forces(ui i);                              //Calculate the forces for particle i>j with atomics
    void calc_forces();                                                 //Calculate the forces between interacting particles
    void recalc_forces();                                               //Recalculate the forces between interacting particles for Velocity Verlet
    void update_boundaries();                                           //Shifts the periodic boxes appropriately for sheared BC
    void periodicity();                                                 //Called after integration to keep the particle within the defined boundaries
    void thread_periodicity_periodic(ui d,ui i);                        //Called by periodicity to keep periodic boundary conditions
    void thread_periodicity_boxshear(ui d,ui i);                        //Called by periodicity to keep boxshear boundary conditions
    void thread_periodicity_hard(ui d,ui i);                            //Called by periodicity to keep hard boundary conditions
    void thread_seuler(ui i);                                           //Symplectic euler integrator (threaded)
    void thread_vverlet_x(ui i);                                        //Velocity verlet integrator for position (threaded)
    void thread_vverlet_dx(ui i);                                       //Velocity verlet integrator for velocity (threaded)
    virtual void integrate();                                           //Integrate particle trajectoriess
    void timestep();                                                    //Do one timestep
    void timesteps(ui k);                                               //Do multiple timesteps
    void import_pos(ldf *x);                                            //Load positions from arrays
    template<typename...arg> void import_pos(ldf *x,arg...argv);        //Load positions from arrays
    void import_pos(ui i,ldf x);                                        //Load position for i from value
    template<typename...arg> void import_pos(ui i,ldf x,arg...argv);    //Load position for i from value
    void import_vel(ldf *dx);                                           //Load velocity from arrays
    template<typename...arg> void import_vel(ldf *dx,arg...argv);       //Load velocity from arrays
    void import_vel(ui i,ldf dx);                                       //Load velocity for i from value
    template<typename...arg> void import_vel(ui i,ldf dx,arg...argv);   //Load velocity for i from value
    void import_force(ldf *F);                                          //Load forces from arrays
    template<typename...arg> void import_force(ldf *F,arg...argv);      //Load forces from arrays
    void import_force(ui i,ldf F);                                      //Load position for i from value
    template<typename...arg> void import_force(ui i,ldf F,arg...argv);  //Load position for i from value
    void export_pos(ldf *x);                                            //Save positions from arrays
    template<typename...arg> void export_pos(ldf *x,arg...argv);        //Save positions to arrays
    void export_pos(ui i,ldf &x);                                       //Save positions from arrays
    template<typename...arg> void export_pos(ui i,ldf &x,arg...argv);   //Save positions to arrays
    void export_vel(ldf *dx);                                           //Save velocity from arrays
    template<typename...arg> void export_vel(ldf *dx,arg...argv);       //Save velocity to arrays
    void export_vel(ui i,ldf &dx);                                      //Save positions from arrays
    template<typename...arg> void export_vel(ui i,ldf &dx,arg...argv);  //Save positions to arrays
    void export_force(ldf *F);                                          //Save forces from arrays
    template<typename...arg> void export_force(ldf *F,arg...argv);      //Save forces to arrays
    void export_force(ui i,ldf &F);                                     //Save forces from arrays
    template<typename...arg> void export_force(ui i,ldf &F,arg...argv); //Save forces to arrays
    ldf direct_readout(ui d,ui i,uc type);                              //Directly readout a position'x'/velocity'v'/forces'F'
    ldf direct_readout(ui i,uc type);                                   //Directly readout a position'x'/velocity'v'/forces'F'
    void add_particle(ldf mass=1.0,ui ptype=0,bool fixed=false);        //Add a particle to the system
    void rem_particle(ui particlenr);                                   //Remove a particle from the system
    void clear();                                                       //Clear all particles and interactions
    void set_damping(ldf coefficient);                                  //Enables damping and sets damping coefficient
    void unset_damping();                                               //Disables damping
    void uitopptr(vector<particle<dim>*> *x,vector<ui> i);              //Convert vector of unsigned integers to particle pointers
    void add_bond(ui p1,ui p2,ui itype,vector<ldf> *params);            //Add a bond to the system of arbitrary type
    void add_spring(ui p1, ui p2,ldf springconstant,ldf l0);            //Add a harmonic bond to the system
    bool share_bond(ui p1,ui p2);                                       //Test whether particles p1 and p2 share a bond
    bool rem_bond(ui p1,ui p2);                                         //Remove a bond from the system
    bool mod_bond(ui p1,ui p2,ui itype,vector<ldf> *params);            //Modify a bond in the system
    ldf thread_H(ui i);                                                 //Measure Hamiltonian for particle i
    ldf thread_T(ui i);                                                 //Measure kinetic energy for particle i
    ldf thread_V(ui i);                                                 //Measure potential energy for particle i
    ldf H();                                                            //Measure Hamiltonian
    ldf T();                                                            //Measure kinetic energy
    ldf V();                                                            //Measure potential energy
};

//This structure defines the Monge patch manifold and its properties
template<ui dim> struct mp
{
    vector<ldf> parameters;                                             //Monge patch function parameters
    fmpptr fmp;                                                         //Monge patch function
    dfmpptr dfmp;                                                       //Derivatives of monge function
    ddfmpptr ddfmp;                                                     //Second derivatives of monge function
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mp();                                                               //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void setmp(ui i=1);                                                 //Picks one of the builtin Monge patches
    void setmp(fmpptr f,dfmpptr df,ddfmpptr ddf);                       //Picks a custom Monge patch
    ldf f(ldf x[dim]);                                                  //Monge patch
    ldf df(ui i,ldf x[dim]);                                            //Monge patch derivative
    ldf g(ui i,ui j,ldf x[dim]);                                        //Monge patch metric tensor
    ldf ginv(ui i,ui j,ldf x[dim]);                                     //Monge patch metric tensor inverse
    ldf dg(ui s,ui i,ui j,ldf x[dim]);                                  //Derivatives of metric
};

//This structure takes care of Monge patch molecular dynamics simulations
template<ui dim> struct mpmd:md<dim>
{
    mp<dim> patch;                                                      //Geometric monge patch information
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mpmd();                                                             //Constructor
    mpmd(ui particlenr);                                                //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    using md<dim>::N;
    using md<dim>::simbox;
    using md<dim>::particles;
    using md<dim>::network;
    using md<dim>::indexdata;
    using md<dim>::v;
    using md<dim>::f;
    using md<dim>::integrator;
    using md<dim>::parallel;
    using md<dim>::periodicity;
    using md<dim>::thread_seuler;
    using md<dim>::thread_vverlet_x;
    using md<dim>::thread_vverlet_dx;
    using md<dim>::recalc_forces;
    using md<dim>::distsq;
    using md<dim>::dd;
    using md<dim>::dap;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf embedded_distsq(ui p1,ui p2);                                   //Calculate distances between two particles (squared)
    ldf embedded_dd_p1(ui i,ui p1,ui p2);                               //Calculate particles relative particle in certain dimension i wrt p1
    ldf embedded_dd_p2(ui i,ui p1,ui p2);                               //Calculate particles relative particle in certain dimension i wrt p2
    void zuiden_C(ui i,ldf C[dim]);                                     //Calculates $g{\rho \sigma} C_{\sigma}$ for particle i of the van Zuiden integrator
    void zuiden_A(ui i,ldf eps[dim]);                                   //Calculates $g{\rho \sigma} A_{\sigma \mu \nu} \epsilon^{\mu} \epsilon^{\nu}$ for particle i of the van Zuiden integrator
    void thread_zuiden_wfi(ui i);                                       //The van Zuiden integrator without fixed point itterations
    void thread_zuiden_protect(ui i);                                   //The van Zuiden integrator with protected fixed point itterations (makes sure you don't get stuck in a loop)
    void thread_zuiden(ui i);                                           //The van Zuiden integrator for Riemannian manifolds (fails for pseudo-Riemannian manifolds)
    #if __cplusplus > 199711L
    void thread_calc_forces(ui i) override;                             //Calculate the forces for particle i>j with atomics
    void integrate() override;                                          //Integrate particle trajectoriess
    #else
    #warning "warning: C++11 not found, disabling override, the mpmd is now broken!"
    void thread_calc_forces(ui i);                                      //Calculate the forces for particle i>j with atomics
    void integrate();                                                   //Integrate particle trajectoriess
    #endif
};

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of libmd HEADER file                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

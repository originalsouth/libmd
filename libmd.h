///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of libmd HEADER file                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef libmd_h
#define libmd_h

#include <cstdio>                                                       //Standard input output (faster than IOstream and also threadsafe) (C)
#include <cstdlib>                                                      //Standard library (C)
#include <cstdarg>                                                      //Support for variadic functions (C)
#include <cmath>                                                        //Standard math library  (C)
#include <cstring>                                                      //Memcpy and memmove support (C)
#include <vector>                                                       //Vector support (C++)
#include <map>                                                          //Map support (C++)
#include <list>                                                         //List support (C++)
#include <utility>                                                      //Pair support (C++)
#include <limits>                                                       //Limits of types (C++)
#include <thread>                                                       //Thread support (C++11)
#include <mutex>                                                        //Mutex support (C++11)
#include <future>                                                       //Future support (C++11)

using namespace std;                                                    //Using standard namespace
typedef long double ldf;                                                //long double is now aliased as ldf
typedef unsigned int ui;                                                //unsigned int is now aliased as ui
typedef unsigned char uc;                                               //unsigned int is now aliased as uc
typedef ldf (*potentialptr)(ldf,ldf,vector<ldf> *);                     //Function pointer to potential functions is now called potential
typedef ldf (*fmpptr)(ldf *,vector<ldf> *);                             //Monge patch function pointer
typedef ldf (*dfmpptr)(ui,ldf *,vector<ldf> *);                         //Monge patch function derivative pointer
typedef ldf (*ddfmpptr)(ui,ui,ldf *,vector<ldf> *);                     //Monge patch function second derivative pointer


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
};

//This structure contains information about the simulation box
//TODO: Jayson LEES EDWARDS
//TODO: Deformations (talk to Jayson)
//TODO: Wall (talk to Jayson)
template<ui dim> struct box
{
    ldf L[dim];                                                         //Box size
    uc bcond[dim];                                                      //Boundary conditions in different dimensions NONE/PERIODIC/HARD(/LEESEDWARDS)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    box();                                                              //Constructor
};

//This structure saves the particle type interactions and calculates the the potentials
//TODO: Implement pair potentials and their d/dr
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

//This structure stores all interactions and their types
struct interact
{
    bool update;                                                        //Should we update the network
    ldf rco;                                                            //R_cuttoff radius
    ldf rcosq;                                                          //R_cuttoff radius squared
    ldf sszsq;                                                          //Skin radius squared
    vector<vector<interactionneighbor>> skins;                          //Particle skin by index (array of vector)
    vector<interactiontype> library;                                    //This is the interaction library
    vector<pair<ui,ui>> backdoor;                                       //Inverse lookup device
    map<pair<ui,ui>,ui> lookup;                                         //This is the interaction lookup device
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pair<ui,ui> hash(ui type1,ui type2);                                //Hash function
    bool probe(ui type1,ui type2);                                      //Check if a typeinteraction exists between two types
};

//This structure takes care of pair potentials (who live outside of the class)
struct pairpotentials
{
    vector<potentialptr> potentials;                                    //Pair potential vector
    vector<potentialptr> dpotentialsdr;                                 //Pair potential d/dr vector
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pairpotentials();                                                   //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(potentialptr p,potentialptr dpdr);                           //Add a potential
    ldf operator()(ui type,ldf r,ldf rsq,vector<ldf>* parameters);      //Pair potential executer
    ldf dr(ui type,ldf r,ldf rsq,vector<ldf>* parameters);              //Pair potential d/dr executer
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
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        celldatatype();                                                 //Constructor
        ~celldatatype();                                                //Destructor
    };
    celldatatype celldata;                                              //Cell data object
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    indexer();                                                          //Constructor
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
    integrators integrator;                                             //Integration method
    threads parallel;                                                   //Multithreader
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    md();                                                               //Constructor
    md(ui particlenr);                                                  //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf distsq(ui p1,ui p2);                                            //Calculate distances between two particles (squared)
    ldf dd(ui i,ui p1,ui p2);                                           //Caculate particles relative particle in certain dimension i
    bool add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Modify type interaction rule
    bool rem_typeinteraction(ui type1,ui type2);                        //Delete type interaction rule //TODO:
    void thread_index(ui i);                                            //Find neighbors per cell i (Or whatever Thomas prefers)
    void index();                                                       //Find neighbors
    void cell();                                                        //Cell indexing algorithm
    void bruteforce();                                                  //Bruteforce indexing algorithm
    void thread_clear_forces(ui i);                                     //Clear forces for particle i
    virtual void thread_calc_forces(ui i);                              //Calculate the forces for particle i>j with atomics
    void calc_forces();                                                 //Calculate the forces between interacting particles
    void recalc_forces();                                               //Recalculate the forces between interacting particles for Velocity Verlet
    void thread_periodicity(ui i);                                      //Called after integration to keep the particle within the defined boundaries
    void thread_seuler(ui i);                                           //Symplectic euler integrator (threaded)
    void thread_vverlet_x(ui i);                                        //Velocity verlet integrator for position (threaded)
    void thread_vverlet_dx(ui i);                                       //Velocity verlet integrator for velocity (threaded)
    virtual void integrate();                                           //Integrate particle trajectoriess
    void timestep();                                                    //Do one timestep
    void timesteps(ui k);                                               //Do multiple timesteps
    void import_pos(...);                                               //Load positions from arrays
    void import_vel(...);                                               //Load velocity from arrays
    void import_force(...);                                             //Load forces from arrays
    void export_pos(...);                                               //Load positions from arrays
    void export_vel(...);                                               //Load velocity from arrays
    void export_force(...);                                             //Load forces from arrays
    void add_particle(ldf mass=1.0,ui ptype=0,bool fixed=false);        //Add a particle to the system
    void rem_particle(ui particlenr);                                   //Remove a particle from the system
    void clear();                                                       //Clear all particles and interactions
    void add_bond();                                                    //Add a bond to the system //TODO: Jayson
    void add_bonds();                                                   //Add multiple bond to the system //TODO: Jayson
    void rem_bond();                                                    //Remove a bond to the system //TODO: Jayson
    void rem_bonds();                                                   //Remove multiple bond to the system //TODO: Jayson
    void mod_bond();                                                    //Modify a bond to the system //TODO: Jayson
    void mod_bonds();                                                   //Modify multiple bond to the system //TODO: Jayson
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
    using md<dim>::integrator;
    using md<dim>::parallel;
    using md<dim>::thread_periodicity;
    using md<dim>::thread_seuler;
    using md<dim>::thread_vverlet_x;
    using md<dim>::thread_vverlet_dx;
    using md<dim>::recalc_forces;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf embedded_distsq(ui p1,ui p2);                                   //Calculate distances between two particles (squared)
    ldf embedded_dd_p1(ui i,ui p1,ui p2);                               //Calculate particles relative particle in certain dimension i wrt p1
    ldf embedded_dd_p2(ui i,ui p1,ui p2);                               //Calculate particles relative particle in certain dimension i wrt p2
    void zuiden_C(ui i,ldf C[dim]);                                     //Calculates $g{\rho \sigma} C_{\sigma}$ for particle i of the van Zuiden integrator
    void zuiden_A(ui i,ldf eps[dim]);                                   //Calculates $g{\rho \sigma} A_{\sigma \mu \nu} \epsilon^{\mu} \epsilon^{\nu}$ for particle i of the van Zuiden integrator
    void thread_zuiden_wfi(ui i);                                       //The van Zuiden integrator without fixed point itterations
    void thread_zuiden_protect(ui i);                                   //The van Zuiden integrator with protected fixed point itterations (makes sure you don't get stuck in a loop)
    void thread_zuiden(ui i);                                           //The van Zuiden integrator for Riemannian manifolds (fails for pseudo-Riemannian manifolds)
    void thread_calc_forces(ui i) override final;                       //Calculate the forces for particle i>j with atomics
    void integrate() override final;                                    //Integrate particle trajectoriess
};

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of libmd HEADER file                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

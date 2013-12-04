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
#include <utility>                                                      //Pair support (C++)
#include <limits>                                                       //Limits of types (C++)
#include <thread>                                                       //Thread support (C++11)
#include <mutex>                                                        //Mutex support (C++11)

using namespace std;                                                    //Using standard namespace
typedef long double ldf;                                                //long double is now aliased as ldf
typedef unsigned int ui;                                                //unsigned int is now aliased as ui
typedef unsigned char uc;                                               //unsigned int is now aliased as uc
typedef ldf (*potentialptr)(ldf,ldf,vector<ldf> *);                     //Function pointer to potential functions is now called potential
typedef ldf (*metricptr)(ldf *);                                        //Function pointer to metric element that defines the curvature

//Potential declarations
ldf COULOMB(ldf r,ldf rsq,vector<ldf> *parameters);
ldf YUKAWA(ldf r,ldf rsq,vector<ldf> *parameters);
ldf HOOKIAN(ldf r,ldf rsq,vector<ldf> *parameters);
ldf LJ(ldf r,ldf rsq,vector<ldf> *parameters);
ldf MORSE(ldf r,ldf rsq,vector<ldf> *parameters);

//TODO: Automatic differentation?
//Potential derivative declaretions
ldf dCOULOMBdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dYUKAWAdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dHOOKIANdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dLJdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dMORSEdr(ldf r,ldf rsq,vector<ldf> *parameters);

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

//This structure defines the molecular dynamics simulation
template<ui dim> struct md
{
    ui N;                                                               //Number of particles
    box<dim> simbox;                                                    //Simulation box
    vector<particle<dim>> particles;                                    //Particle array
    interact network;                                                   //Interaction network
    pairpotentials v;                                                   //Pair potential functor
    integrators integrator;                                             //Integration method
    threads parallel;                                                   //Multithreader
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    md();                                                               //Constructor
    md(ui particlenr);                                                  //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf distsq(ui p1,ui p2);                                            //Calculate distances between two particles (squared)
    ldf dd(ui i,ui p1,ui p2);                                           //Caculate particles relative particle in certain dimension i
    void add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Add type interaction rule
    void mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Modify type interaction rule
    void rem_typeinteraction(ui type1,ui type2);                        //Delete type interaction rule //TODO:
    void thread_index(ui i);                                            //Find neighbors per cell i (Or whatever Thomas prefers)
    void index();                                                       //Find neighbors
    void thread_clear_forces(ui i);                                     //Clear forces for particle i
    void thread_calc_forces(ui i);                                      //Calculate the forces for particle i>j with atomics
    void calc_forces();                                                 //Calculate the forces between interacting particles
    void recalc_forces();                                               //Recalculate the forces between interacting particles for Velocity Verlet
    void thread_periodicity(ui i);                                      //Called after integration to keep the particle within the defined boundaries
    void thread_seuler(ui i);                                           //Symplectic euler integrator (threaded)
    void thread_vverlet_x(ui i);                                        //Velocity verlet integrator for position (threaded)
    void thread_vverlet_dx(ui i);                                       //Velocity verlet integrator for velocity (threaded)
    void integrate();                                                   //Integrate particle trajectoriess
    void timestep();                                                    //Do one timestep
    void timesteps(ui k);                                               //Do multiple timesteps
    void import_pos(...);                                               //Load positions from arrays
    void import_vel(...);                                               //Load velocity from arrays
    void import_force(...);                                             //Load forces from arrays
    void export_pos(...);                                               //Load positions from arrays
    void export_vel(...);                                               //Load velocity from arrays
    void export_force(...);                                             //Load forces from arrays
    void add_particle();                                                //Add a particle to the system //TODO:
    void rem_particle();                                                //Remove a particle from the system //TODO:
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

template<ui dim> struct geometry
{
    bool curvature;                                                     //Is there any curvature
    //add stuff here
};

//TODO:
//This structure takes care of Monge patch molecular dynamics
template<ui dim> struct mpmd
{
    ui N;                                                               //Number of particles
    ui nothreads;                                                       //Number of threads
    box<dim+1> simbox;                                                  //Simulation box
    geometry<dim> manifold;                                             //Geometric information
    vector<particle<dim>> particles;                                    //Particle array
    interact network;                                                   //Interaction network
    pairpotentials v;                                                   //Pair potential functor
    integrators integrator;                                             //Integration method
    threads parallel;                                                   //Multithreader
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mpmd();                                                             //Constructor
    mpmd(ui particlenr);                                                //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf distsq(ui p1,ui p2);                                            //Calculate distances between two particles (squared)
    ldf dd(ui i,ui p1,ui p2);                                           //Caculate particles relative particle in certain dimension i
    void add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Add type interaction rule
    void mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Modify type interaction rule
    void rem_typeinteraction(ui type1,ui type2);                        //Delete type interaction rule //TODO:
    void thread_index(ui i);                                            //Find neighbors per cell i (Or whatever Thomas prefers)
    void index();                                                       //Find neighbors
    void thread_clear_forces(ui i);                                     //Clear forces for particle i
    void thread_calc_forces(ui i);                                      //Calculate the forces for particle i>j with atomics
    void calc_forces();                                                 //Calculate the forces between interacting particles
    void recalc_forces();                                               //Recalculate the forces between interacting particles for Velocity Verlet
    void thread_periodicity(ui i);                                      //Called after integration to keep the particle within the defined boundaries
    void thread_seuler(ui i);                                           //Symplectic euler integrator (threaded)
    void thread_vverlet_x(ui i);                                        //Velocity verlet integrator for position (threaded)
    void thread_vverlet_dx(ui i);                                       //Velocity verlet integrator for velocity (threaded)
    void integrate();                                                   //Integrate particle trajectoriess
    void timestep();                                                    //Do one timestep
    void timesteps(ui k);                                               //Do multiple timesteps
    void import_pos(...);                                               //Load positions from arrays
    void import_vel(...);                                               //Load velocity from arrays
    void import_force(...);                                             //Load forces from arrays
    void export_pos(...);                                               //Load positions from arrays
    void export_vel(...);                                               //Load velocity from arrays
    void export_force(...);                                             //Load forces from arrays
    void add_particle();                                                //Add a particle to the system //TODO:
    void rem_particle();                                                //Remove a particle from the system //TODO:
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

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of libmd HEADER file                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

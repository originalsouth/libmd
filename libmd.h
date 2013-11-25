///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of libmd HEADER file                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef libmd_h
#define libmd_h

#include <cstdio>                                                       //Standard input output (faster than IOstream and also threadsafe) (C)
#include <cstdlib>                                                      //Standard library (C)
#include <cstdarg>                                                      //Support for variadic functions (C)
#include <cmath>                                                        //Standard math library  (C)
#include <vector>                                                       //Vector support (C++)
#include <map>                                                          //Map support (C++)
#include <utility>                                                      //Pair support (C++)
#include <limits>                                                       //Limits of types (C++)
#include <thread>                                                       //Thread support (C++11)
#include <future>                                                       //Future support (C++11)
#include <atomic>                                                       //Atomic support (C++11)

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

//Potential derivative declaretions
ldf dCOULOMBdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dYUKAWAdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dHOOKIANdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dLJdr(ldf r,ldf rsq,vector<ldf> *parameters);
ldf dMORSEdr(ldf r,ldf rsq,vector<ldf> *parameters);

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
    particle();                                                         //Constructor
    particle(ldf mass,ui ptype,bool fixed);                             //Constructor
};

//This structure contains information about the simulation box
//TODO: Jayson LEES EDWARDS
//TODO: Deformations (talk to Jayson)
//TODO: Wall (talk to Jayson)
template<ui dim> struct box
{
    ldf L[dim];                                                         //Box size
    uc bcond[dim];                                                      //Boundary conditions in different dimensions NONE/PERIODIC/HARD(/LEESEDWARDS)
};

//This structure saves the particle type interactions and calculates the the potentials
//TODO: Implement pair potentials and their d/dr
struct interactiontype
{
    ui potential;                                                       //Type of potential
    vector<ldf> parameters;                                             //Parameters of potential
    ldf vco;                                                            //Cuttoff potential energy
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interactiontype();                                                  //Constructor
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
    ui nothreads;                                                       //Number of threads
    box<dim> simbox;                                                    //Simulation box
    vector<particle<dim>> particles;                                    //Particle array
    interact network;                                                   //Interaction network
    pairpotentials v;                                                   //Pair potential functor
    integrators integrator;                                             //Integration method
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    md();                                                               //Constructor
    md(ui particlenr);                                                  //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ldf distsq(ui p1,ui p2);                                            //Calculate distances between two particles (squared)
    ldf dd(ui i,ui p1,ui p2);                                           //Caculate particles relative particle in certain dimension i
    void add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Add type interaction rule //FIXME:
    void mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);   //Modify type interaction rule //TODO:
    void rem_typeinteraction(ui type1,ui type2);                        //Delete type interaction rule //FIXME:
    void thread_index(ui i);                                            //Find neighbors per cell i (Or whatever Thomas prefers)
    void index();                                                       //Find neighbors
    void thread_clear_forces(ui i);                                     //Clear forces for particle i
    void thread_calc_forces(ui i);                                      //Calculate the forces for particle i>j with atomics
    void calc_forces();                                                 //Calculate the forces between interacting particles
    void recalc_forces();                                               //Recalculate the forces between interacting particles for Velocity Verlet
    void thread_integrate(ui i,ui gen);                                 //Integrate trajectory position for particle i
    void integrate();                                                   //Integrate particle trajectoriess
    void timestep();                                                    //Do one timestep
    void timesteps(ui k);                                               //Do multiple timesteps
    void import_pos(...);                                               //Load positions from arrays //TODO:
    void import_vel(...);                                               //Load velocity from arrays //TODO:
    void import_force(...);                                             //Load forces from arrays //TODO:
    void export_pos(...);                                               //Load positions from arrays //TODO:
    void export_vel(...);                                               //Load velocity from arrays //TODO:
    void export_force(...);                                             //Load forces from arrays //TODO:
    void add_particle();                                                //Add a particle to the system //TODO:
    void rem_particle();                                                //Remove a particle from the system //TODO:
    ldf thread_H(ui i);                                                 //Measure Hamiltonian for particle i //TODO:
    ldf thread_T(ui i);                                                 //Measure kinetic energy for particle i //TODO:
    ldf thread_V(ui i);                                                 //Measure potential energy for particle i //TODO:
    ldf H();                                                            //Measure Hamiltonian //TODO:
    ldf T();                                                            //Measure kinetic energy //TODO:
    ldf V();                                                            //Measure potential energy //TODO:
};

template<ui dim> struct geometry
{
    bool curvature;                                                     //Is there any curvature
    //add stuff here
};

template<ui dim> struct cmd
{

};

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of libmd HEADER file                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

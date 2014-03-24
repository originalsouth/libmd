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
#include <stack>                                                        //Stack support (C++)
#include <set>                                                          //Set support (C++)
#include <unordered_set>                                                //Unordered set support (C++11)
#include <utility>                                                      //Pair support (C++)
#include <limits>                                                       //Limits of types (C++)
#include <thread>                                                       //Thread support (C++11)
#include <mutex>                                                        //Mutex support (C++11)
#include <future>                                                       //Future support (C++11)
#include <algorithm>                                                    //Algorithm support (C++)

#ifdef FE
#include <cfenv>                                                        //Floating point exception handling (C/C++11)
#endif

using namespace std;                                                    //Using standard namespace

typedef long double ldf;                                                //long double is now aliased as ldf
typedef unsigned int ui;                                                //unsigned int is now aliased as ui
typedef unsigned char uc;                                               //unsigned char is now aliased as uc

enum INTEGRATOR:uc {SEULER,VVERLET};                                    //Integration options
enum MP_INTEGRATOR:uc {MP_VZ,MP_VZ_P,MP_VZ_WFI,MP_SEULER,MP_VVERLET};   //Monge patch integration options
enum BCOND:uc {NONE,PERIODIC,HARD,BOXSHEAR};                            //Boundary condition options
enum INDEX:uc {CELL,BRUTE_FORCE,KD_TREE};                               //Indexing options
enum POT:ui                                                             //Potential options
{
    POT_COULOMB,
    POT_YUKAWA,
    POT_HOOKEAN,
    POT_LJ,
    POT_MORSE,
    POT_FORCEDIPOLE,
    POT_HOOKEANFORCEDIPOLE,
    POT_ANHARMONICSPRING
};
///External force options
enum EXTFORCE:ui
{
    EXTFORCE_DAMPING,
    EXTFORCE_DISSIPATION
};
///Monge patch options
enum MP:ui
{
    MP_FLATSPACE,
    MP_GAUSSIANBUMP,
    MP_EGGCARTON,
    MP_MOLLIFIER
};

//These functions defined outside of libmd
void __libmd__info();                                                   //Basic libmd comilation info

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

///This structure saves the external force functions and calculates them
struct forcetype
{
    ui externalforce;                                                   //External force type
    vector<vector<ui>> particles;                                       //Interacting particle list
    vector<ldf> parameters;                                             //Parameters for the external force
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    forcetype(ui noexternalforce,vector<vector<ui>> *plist,vector<ldf> *param); //Constructor
};

//This structure introduces "super_particles" i.e. particles that built from sub_particles
struct superparticle
{
    map<ui,ui> particles;                                               //Particles in super particles
    ui sptype;                                                          //Super particle type
};

//This structure caries a lookup device for a specific super particle type
struct superparticletype
{
    map<pair<ui,ui>,ui> splookup;                                       //This is the interaction lookup device
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
    unordered_set<ui> free_library_slots;                               //Stores free library slots
    map<pair<ui,ui>,ui> lookup;                                         //This is the interaction lookup device
    vector<ui> spid;                                                    //Super particle identifier array
    vector<superparticle> superparticles;                               //Actual super particle array
    vector<superparticletype> sptypes;                                  //Super particle type array
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    interact();                                                         //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pair<ui,ui> hash(ui type1,ui type2);                                //Hash function
    bool probe(ui type1,ui type2);                                      //Check if a typeinteraction exists between two types
};

//Potential functions
template<class X> X COULOMB(X r,vector<ldf> *parameters);
template<class X> X YUKAWA(X r,vector<ldf> *parameters);
template<class X> X HOOKEAN(X r,vector<ldf> *parameters);
template<class X> X MORSE(X r,vector<ldf> *parameters);
template<class X> X FORCEDIPOLE(X r,vector<ldf> *parameters);
template<class X> X HOOKEANFORCEDIPOLE(X r,vector<ldf> *parameters);
template<class X> X ANHARMONICSPRING(X r,vector<ldf> *parameters);

//External force functions
template<ui dim> void DAMPING(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters);
template<ui dim> void DISSIPATION(particle<dim> *p,vector<particle<dim>*> *particles,vector<ldf> *parameters);

//This structure automatically differentiates first order
struct dual
{
    ldf x;                                                              //Function value
    ldf dx;                                                             //Function derivative value
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dual();                                                             //Constructor
    dual(ldf f,ldf fx=0.0);                                             //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dual operator=(dual y);                                             //Assign operator
    void operator+=(dual y);                                            //Add-assign operator
    void operator-=(dual y);                                            //Subtract-assign operator
    template<class X> operator X();                                     //Cast overload
    template<class X> X operator=(X y);                                 //Assign foreign type operator
    template<class X> void operator+=(X y);                             //Add-assign foreign type operator
    template<class X> void operator-=(X y);                             //Subtract-assign foreign type operator
    template<class X> void operator*=(X y);                             //Multiply-assign foreign type operator
    template<class X> void operator/=(X y);                             //Devide-assign foreign type operator
    template<class X> bool operator==(X y);                             //Test equality to foreign type
    template<class X> bool operator<=(X y);                             //Test if smaller or equal than foreign type
    template<class X> bool operator>=(X y);                             //Test if greater or equal than foreign type
    template<class X> bool operator<(X y);                              //Test if smaller than foreign type
    template<class X> bool operator>(X y);                              //Test if greater than foreign type
};

template<class X> using potentialptr=X (*)(X,vector<ldf> *);            //Function pointer to potential functions is now called potentialptr

//This structure takes care of pair potentials (who live outside of the class)
struct pairpotentials
{
    vector<potentialptr<dual>> potentials;                              //Pair potential vector
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pairpotentials();                                                   //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ui add(potentialptr<dual> p);                                       //Add a potentials
    ldf operator()(ui type,ldf r,vector<ldf> *parameters);              //Pair potential executer
    ldf dr(ui type,ldf r,vector<ldf> *parameters);                      //Pair potential d/dr executer
};

template<ui dim> using extforceptr=void (*)(particle<dim> *,vector<particle<dim>*> *,vector<ldf> *); //Function pointer to external force functions is now called extforceptr

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
    struct kdtreedatatype
    {
        ui (*Idx);                                                      //Indices of particles, ordered by tree-structure
        ui DivideByDim[30];                                             //Dimension to divide by at each recursion level (assuming N <= 2^30)
        ldf (*Pmin)[dim],(*Pmax)[dim];                                  //Minimum and maximum value of each coordinate, for every subtree
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        kdtreedatatype();                                               //Constructor
        ~kdtreedatatype();                                              //Destructor
    };
    kdtreedatatype kdtreedata;
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
    md(ui particlenr=0);                                                //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void init(ui particlenr);                                           //Copy of the particle number constructor
    ldf dap(ui i,ldf ad);                                               //Manipulate particle distances with respect to periodic boundary conditions
    ldf distsq(ui p1,ui p2);                                            //Calculate distances between two particles (squared)
    ldf dd(ui i,ui p1,ui p2);                                           //Calculate difference in particle positions in certain dimension i
    void all_interactions(vector<pair<ui,ui>> &table);                  //Dump all interaction into a table
    ui add_interaction(ui potential,vector<ldf> *parameters);           //Add type interaction rule
    bool mod_interaction(ui interaction,ui potential,vector<ldf> *parameters);//Modify type interaction rule
    bool rem_interaction(ui interaction);                               //Delete type interaction rule
    bool add_typeinteraction(ui type1,ui type2,ui interaction);         //Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui interaction);         //Modify type interaction rule
    void mad_typeinteraction(ui type1,ui type2,ui interaction);         //Force add/mod type interaction rule
    bool add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);//Add type interaction rule
    bool mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);//Modify type interaction rule
    void mad_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters);//Force add/mod type interaction rule
    bool rem_typeinteraction(ui type1,ui type2);                        //Delete type interaction rule
    ui add_sp_interaction(ui spt,ui p1,ui p2,ui interaction);           //Add type interaction rule
    bool mod_sp_interaction(ui spt,ui p1,ui p2,ui interaction);         //Modify type interaction rule
    bool rem_sp_interaction(ui spt,ui p1,ui p2);                        //Delete type interaction rule
    bool rem_sp_interaction(ui spt);                                    //Delete type interaction rule
    ui add_forcetype(ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters);//Add force type
    bool mod_forcetype(ui notype,ui force,vector<vector<ui>> *noparticles,vector<ldf> *parameters);//Modify force type
    bool rem_forcetype(ui notype);                                      //Delete force type
    void assign_forcetype(ui particlenr,ui ftype);                      //Assign force type to particle
    void assign_all_forcetype(ui ftype);                                //Assign force type to all particles
    void unassign_forcetype(ui particlenr,ui ftype);                    //Unassign force type to particle
    void unassign_all_forcetype(ui ftype);                              //Unassign force type to all particles
    void clear_all_assigned_forcetype();                                //Clear all assigned forces
    void set_rco(ldf rco);                                              //Sets the cuttoff radius and its square
    void set_ssz(ldf ssz);                                              //Sets the skin size radius and its square
    void set_reserve(ldf ssz);                                          //Set reserve memory according to skin size
    void set_reserve(ldf ssz,ui M);                                     //Set reserve memory according to skin size and some arbitrary number of particles
    void set_type(ui p, ui newtype);                                    //Update the type associated with particle p
    void set_index_method(ui method);                                   //Set indexmethod
    void thread_index(ui i);                                            //Find neighbors per cell i
    void index();                                                       //Find neighbors
    bool test_index();                                                  //Test if we need to run the indexing algorithm
    void thread_index_stick(ui i);                                      //Save the particle position at indexing
    ui kdtree_build (ui first, ui last, ui level);                      //k-d tree indexing algorithm: tree build function (recursive)
    void kdtree_index (ui first1, ui last1, ui first2, ui last2);       //k-d tree indexing algorithm: neighbor finder (recursive)
    void kdtree();                                                      //k-d tree indexing algorithm
    void cell();                                                        //Cell indexing algorithm
    void thread_cell (ui i);                                            //Cell indexer for cell i (thread)
    void bruteforce();                                                  //Bruteforce indexing algorithm
    void skinner(ui i,ui j);                                            //Places interactionneighbor in skin
    void thread_clear_forces(ui i);                                     //Clear forces for particle i
    virtual void thread_calc_forces(ui i);                              //Calculate the forces for particle i>j with atomics
    void calc_forces();                                                 //Calculate the forces between interacting particles
    void recalc_forces();                                               //Recalculate the forces between interacting particles for Velocity Verlet
    void update_boundaries();                                           //Shifts the periodic boxes appropriately for sheared BC
    void periodicity();                                                 //Called after integration to keep the particle within the defined boundaries
    void thread_periodicity(ui i);                                      //Apply periodicity to one particle only
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
    void fix_particle(ui i,bool fix);                                   //Fix a particle
    void fix_particles(ui spi,bool fix);                                //Fix a super particles
    ui clone_particle(ui i,ldf x[dim]);                                 //Fix a particle
    ui clone_particles(ui spi,ldf x[dim]);                              //Fix a particle
    void translate_particle(ui i,ldf x[dim]);                           //Translate (or move) a particle
    void translate_particles(ui spi,ldf x[dim]);                        //Translate (or move) a super particle
    void drift_particle(ui i,ldf dx[dim]);                              //Add velocity to a particle
    void drift_particles(ui spi,ldf dx[dim]);                           //Add velocity to a super particle (all particles the same)
    void set_position_particles(ui spi,ldf x[dim]);                     //Get center of mass of super particle
    void set_velocity_particles(ui spi,ldf dx[dim]);                    //Assign velocity to a super particle (all particles the same)
    void get_position_particles(ui spi,ldf x[dim]);                     //Get center of mass of super particle
    void get_velocity_particles(ui spi,ldf dx[dim]);                    //Get average velocity of a super particle
    ui sp_ingest(ui spi,ui i);                                          //Add a particle to a super particle
    ui sp_ingest(ui spi,ui sptype,ui i);                                //Add a particle to a super particle
    void sp_dispose(ui spi);                                            //Remove particle from a super particle
    void sp_p_dispose(ui i);                                            //Remove particle from its super particle
    ui add_particle(ldf mass=1.0,ui ptype=0,bool fixed=false);          //Add a particle to the system
    ui add_particle(ldf x[dim],ldf mass=1.0,ui ptype=0,bool fixed=false);//Add a particle to the system at certain position
    ui add_particle(ldf x[dim],ldf dx[dim],ldf mass=1.0,ui ptype=0,bool fixed=false);//Add a particle to the system at certain position with certain velocity
    void rem_particle(ui i);                                            //Remove a particle from the system
    void rem_particles(ui spi);                                         //Remove a super particle
    void clear();                                                       //Clear all particles and interactions
    void set_damping(ldf coefficient);                                  //Enables damping and sets damping coefficient
    void unset_damping();                                               //Disables damping
    void uitopptr(vector<particle<dim>*> *x,vector<ui> i);              //Convert vector of unsigned integers to particle pointers
    bool add_bond(ui p1,ui p2,ui interaction);                          //Add a bond
    bool mod_bond(ui p1,ui p2,ui interaction);                          //Modify a bond
    void mad_bond(ui p1,ui p2,ui interaction);                          //Force add/modify bond
    bool add_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters);    //Add a bond
    bool mod_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters);    //Modify a bond
    void mad_bond(ui p1,ui p2,ui potential,vector<ldf> *parameters);    //Force add/modify bond
    bool rem_bond(ui p1,ui p2,bool force=false);                        //Remove a bond from the system
    void add_spring(ui p1, ui p2,ldf springconstant,ldf l0);            //Add a harmonic bond to the system
    ldf thread_H(ui i);                                                 //Measure Hamiltonian for particle i
    virtual ldf thread_T(ui i);                                         //Measure kinetic energy for particle i
    virtual ldf thread_V(ui i);                                         //Measure potential energy for particle i
    ldf H();                                                            //Measure Hamiltonian
    ldf T();                                                            //Measure kinetic energy
    ldf V();                                                            //Measure potential energy
};

//Autodiff for Monge patches
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

//Monge patch function pointer
template<class X,ui dim> using fmpptr=X (*)(X x[dim],vector<ldf> *param);

//Monge patches (and related)
ldf kdelta(ui i,ui j);
template<class X,ui dim> X FLATSPACE(X x[dim],vector<ldf> *param);
template<class X,ui dim> X GAUSSIANBUMP(X x[dim],vector<ldf> *param);

//This structure defines the Monge patch manifold and its properties
template<ui dim> struct mp
{
    vector<ldf> parameters;                                             //Monge patch function parameters
    fmpptr<ldf,dim> fmp;                                                //Monge patch function
    fmpptr<duals<dim>,dim> dfmp;                                        //Derivatives of monge function
    vector<duals<dim>> geometryx;
    vector<duals<dim>> geometryxp;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mp();                                                               //Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void setmp(ui i=1);                                                 //Picks one of the builtin Monge patches
    void setmp(fmpptr<ldf,dim> f,fmpptr<duals<dim>,dim> df);            //Picks a custom Monge patch
    void calc(ui i,ldf x[dim]);
    ldf f(ldf x[dim]);                                                  //Monge patch
    ldf g(ui i,ui mu,ui nu);                                            //Monge patch metric tensor
    ldf gp(ui i,ui mu,ui nu);                                           //Monge patch metric tensor
    ldf ginv(ui i,ui mu,ui nu);                                         //Monge patch metric tensor inverse
    ldf G(ui i,ui sigma,ui mu,ui nu);                                   //Monge patch Christoffel symbols (of first kind)
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
    void zuiden_C(ui i,ldf ZC[dim]);                                    //Calculates $g{\rho \sigma} C_{\sigma}$ for particle i of the van Zuiden integrator
    void zuiden_A(ui i,ldf eps[dim]);                                   //Calculates $g{\rho \sigma} A_{\sigma \mu \nu} \epsilon^{\mu} \epsilon^{\nu}$ for particle i of the van Zuiden integrator
    void thread_zuiden_wfi(ui i);                                       //The van Zuiden integrator without fixed point itterations
    void thread_zuiden_protect(ui i);                                   //The van Zuiden integrator with protected fixed point itterations (makes sure you don't get stuck in a loop)
    void thread_zuiden(ui i);                                           //The van Zuiden integrator for Riemannian manifolds (fails for pseudo-Riemannian manifolds)
    void thread_history(ui i);                                          //Set the history of particle i
    void history();                                                     //Set the history of all particles
    void thread_calc_geometry(ui i);                                    //Calculate Monge patch derivatives for partice i
    void calc_geometry();                                               //Calculate Monge patch derivatives
    void thread_calc_forces(ui i) override;                             //Calculate the forces for particle i>j with atomics
    void integrate() override;                                          //Integrate particle trajectoriess
    ldf thread_T(ui i) override;                                        //Calculate kinetic energy of a particle
    ldf thread_V(ui i) override;                                        //Calculate kinetic energy
};

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of libmd HEADER file                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

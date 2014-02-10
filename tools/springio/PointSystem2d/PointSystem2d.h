#ifndef _POINTSYSTEM_H_
#define _POINTSYSTEM_H_

#include <vector>
#include <deque>
#include <set>
#include <map>
#include <list>
#include <string>
#include <stdio.h>
#include <math.h>
#include "Point2d.h"

using namespace std;
typedef vector<Point2d> PointArray;
typedef list<Bond> BondArray; // don't care about indices of bonds, just want to traverse them!
typedef deque<Dihedral> DihedralArray; // each dihedral consists of four indices
typedef vector<Facet> FacetArray;
typedef vector<int> intArray;
typedef vector<double> doubleArray;
typedef vector<intArray > intArrayArray;
//~ typedef vector<double> dblArray;
typedef set<Bond* > BondPointerArray;
typedef set<Facet* > FacetPointerArray;
typedef vector<BondPointerArray > BondPointerArrayArray;
typedef vector<FacetPointerArray > FacetPointerArrayArray;
//~ typedef vector<
typedef vector<PointNeighbours> PointNbrArray;
typedef pair<int, int> PtPair;
typedef map<PtPair, BondArray::iterator> BondPointerIndex;


class PointSystem2d {
    public:
    // parameters
        int N;
        double lx; // box periodic bc in x dirn
        double ly; // box periodic bc in y dirn
        double shearx;
        
        // data objects
        PointArray points;
        BondArray bonds;
        BondArray force_dipoles;
        FacetArray facets;
        PointNbrArray nbrinfo;
        BondPointerIndex bondpointers;
        
        PointSystem2d();
        bool morse;

        void initialize(string, string, double, double, double, double, double); // if using surface evolver curvatures, no need for dihedrals/adjacent facets!
        void updatePoints(string);
        void updateShear(double);
        void addPoint(Point2d);
        void removeLastPt();
        
        PtPair ordered_key(int, int);
        void addBond(int,int,double,double);
        void addBondsFrom(string,double,double);
        void removeBond(int,int);
        void removeBond(Bond);
        void clearBonds();
        
        void bp_update_last(int, int);
        void make_nbrs(int, int);
        
        double conn();
        int randomBondIndex();
        Bond bondFromIdx(int);
        
        bool random_cut(double);
        bool uniform_cut(double);
        
        void random_force_dipoles(double, double); // populate bond network with force dips
        
        // helper functions
        double bondlength(Bond);
        double sphr; // sphere confinement
        double sphpot; // sphere confinement potential
        double pressure;
        
        double espring(Bond);
        
        double estretch();
        double eforcedipole();
        Point2d fstretch(int, int);
        double fmag(int, int);
        double fmag(Bond);
        double energy();
        
        // for conjugate gradient
        void d_estretch(PointArray*);
        void d_eforcedipole(PointArray*);
        void d_energy(PointArray*);                                   
        
        // sphere clamp
        double esphere();
        void d_esphere(PointArray*);
        
        // volume and area
        double surf();
        double vol();
        double elementsurf(Facet);
        double elementvol(Facet);
        double epressure();
        void d_epressure(PointArray*);
        
        // morse potential bonds
        double emorse();
        void d_emorse(PointArray*);
        
        //refpts
        PointArray refpoints;
        BondArray refbonds;
        intArray pinpts;
        void makeref(int, double, double, double); // make reference Point2d and tether
        void makeref_self(int, double); // tether to self (immobilize particle)
        double erefpts(); // tether some points to reference points
        void d_erefpts(PointArray*);
        void clearRefs();


        vector<double> stress_tensor();
        void shake(double);
        
        // i/o
        int write_points(string);
        void write_bonds(string);
        void write_data(string);
        void print_energies(string);
};

#endif

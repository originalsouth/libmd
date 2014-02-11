//      PointSystem2d.cpp
//      
//      Copyright 2011 Jayson Paulose <jpaulose@jpaulose-laptop>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

// Modified 2012/04/05 for use to minimize just a system of springs on sphere

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include "PointSystem2d.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

#define ZEROCORR 0 // correct for extra 1 in data file (if present)
#define ERROR_BONDCUT "error: target connectivity is below current connectivity!"
#define MIN_CONN 3

boost::mt19937 rng;

PointSystem2d::PointSystem2d() {
    points = PointArray();
    bonds = BondArray();
    force_dipoles = BondArray();
    facets = FacetArray();
    refpoints = PointArray();
    refbonds = BondArray();
    pinpts = intArray();
    nbrinfo = PointNbrArray();
    bondpointers = BondPointerIndex();
    N=0;
    morse = false;
    lx=1000;
    ly=1000;
    shearx=0;
    //~ double lambda, gamma, v0;
    //~ bool spontCurv = false;
}


// use this constructor for vertex-based curvature energies (SE or cotangent)
void PointSystem2d::initialize(string ptfile, string bfile, 
                    double Sc=1, double Bl0=-1, double lx0=10., double ly0=10., double shearx0=0.) {
    FILE* inputM;
    double x, y;
    int p1, p2;
    
    lx = lx0;
    ly = ly0;
    shearx = shearx0;
    
    // initialize points
    inputM = fopen(ptfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%lf %lf\n", &x, &y);
        addPoint(Point2d(x,y,lx, ly,shearx));
    } 
    N = points.size();

    // initialize bonds
    inputM = fopen(bfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%d %d\n", &p1, &p2);
        p1 = p1 - ZEROCORR; p2 = p2 - ZEROCORR;
        addBond(p1,p2,Sc,Bl0);
    }
}


// update all points with coords from a particular inputfile
void PointSystem2d::updatePoints(string ptfile) {
    FILE* inputM;
    double x, y;
    int idx = 0;
    
    // initialize points
    inputM = fopen(ptfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%lf %lf\n", &x, &y);
        points[idx].set(x,y,lx,ly,shearx);
        idx++;
    } 
}

void PointSystem2d::updateShear(double shear) {
    shearx = shear;
    for (PointArray::iterator it = points.begin(); it != points.end(); it++)
        it->setshear(shear);
}

void PointSystem2d::addPoint(Point2d pt) {
    points.push_back(pt);
    nbrinfo.push_back(PointNeighbours());
    N++;
}

void PointSystem2d::removeLastPt() {
    points.pop_back();
    N--;
}

void PointSystem2d::clearBonds() {
    bonds.clear();
}

void PointSystem2d::clearRefs() {
    refpoints.clear();
    refbonds.clear();
}

void PointSystem2d::bp_update_last(int p1, int p2) {
    // update bondpointerindex with pointer to last bond
    bondpointers.insert(make_pair(ordered_key(p1,p2), --bonds.end()));
}


void PointSystem2d::make_nbrs(int p1, int p2) {
    nbrinfo[p1].add_nbr(p2);
    nbrinfo[p2].add_nbr(p1);
}


void PointSystem2d::addBond(int p1, int p2, double k=1., double l0 = -1.) {
    if (l0 < 0.) {
        Bond nb = Bond(p1,p2,k,0);
        nb.setlength(bondlength(nb));
        bonds.push_back(nb);
    }
    else bonds.push_back(Bond(p1,p2,k,l0));
    
    // update the bond pointer index
    bp_update_last(p1,p2);
    
    // set as mutual neighbors
    make_nbrs(p1,p2);
}

void PointSystem2d::addBondsFrom(string bfile, double k=1., double l0 = -1.) {
    // initialize bonds
    FILE* inputM = fopen(bfile.c_str(), "r");
    int p1, p2;
    while (!(feof(inputM))) {
        fscanf(inputM, "%d %d\n", &p1, &p2);
        p1 = p1 - ZEROCORR; p2 = p2 - ZEROCORR;
        if (l0 < 1) {
            Bond nb = Bond(p1,p2,k,0);
            nb.setlength(bondlength(nb));
            bonds.push_back(nb);
        }
        else bonds.push_back(Bond(p1,p2,k,l0));
    }
}


PtPair PointSystem2d::ordered_key(int p1, int p2) {
    if (p1 < p2) return make_pair(p1, p2); 
    else return make_pair(p2, p1);
}


void PointSystem2d::removeBond(int p1, int p2) {
    // remove all traces of bond including from the point adjacency lists
    PtPair key_to_remove = ordered_key(p1,p2);
    
    if (bondpointers.count(key_to_remove) > 0) {
        BondArray::iterator el_to_remove = bondpointers[key_to_remove];
        bonds.erase(el_to_remove);
         // remove from bondpointers list
        bondpointers.erase(key_to_remove);
        // remove neighbour reference
        nbrinfo[p1].remove_nbr(p2);
        nbrinfo[p2].remove_nbr(p1);
    }
}

void PointSystem2d::removeBond(Bond b) {
    removeBond(b.p1(), b.p2());
}

double PointSystem2d::conn() {
    return bonds.size()*2./N;
}

double PointSystem2d::bondlength(Bond b) {
    Point2d p1 = points[b.p1()];
    Point2d p2 = points[b.p2()];
    return p1.distanceTo(p2);
}

// bond cutting algos
int randomIndex(int container_size) {
    boost::uniform_int<> dist(0, container_size-1);
    return dist(rng);
}

int PointSystem2d::randomBondIndex() {
    //~ boost::uniform_int<> dist(0, bonds.size()-1);
    //~ return dist(rng);
    return randomIndex(bonds.size());
}

Bond PointSystem2d::bondFromIdx(int idx) {
    // consider boost::container::flat_set if this is slow
    BondArray::iterator it = bonds.begin();
    advance(it, idx);
    return *it;
}


bool PointSystem2d::random_cut(double conn_target) {
    /* Cut bonds at random, subject to no-local-floppiness constraint
     * z >=3
     */
    printf(" old connectivity: %1.4f no. of bonds: %lu\n", conn(), bonds.size());
    if (conn() < conn_target || conn_target < 0) {
        cerr << ERROR_BONDCUT << endl;
        return false;
    }
    
    while (conn() > conn_target) {
        Bond bond_to_cut = bondFromIdx(randomBondIndex());
        if (nbrinfo[bond_to_cut.p1()].degree() > MIN_CONN &&
            nbrinfo[bond_to_cut.p2()].degree() > MIN_CONN) 
                removeBond(bond_to_cut);
    }
    printf(" new connectivity: %1.4f no. of bonds: %lu\n", conn(), bonds.size());
    return (conn() < conn_target);
}


bool PointSystem2d::uniform_cut(double conn_target) {
    /* cut highest-coordinated bonds first */
    
    printf(" old connectivity: %1.4f no. of bonds: %lu\n", conn(), bonds.size());
    if (conn() < conn_target || conn_target < 0) {
        cerr << ERROR_BONDCUT << endl;
        return false;
    }
    
    int maxDegree = 0;
    for (int i = 0; i < N; i++) {
        int degree = nbrinfo[i].degree();
        if (degree > maxDegree) maxDegree = degree;
    }
    
    for (int c1 = maxDegree; c1 > MIN_CONN; c1--) {
        for (int c2 = c1; c2 > MIN_CONN; c2--) {
            // collect all pairs with c1 and c2
            cout << " c1 "<< c1 <<" c2 "<<c2<<endl;
            typedef list<PtPair > ConnSet;
            ConnSet c12bonds = ConnSet();
            for (BondArray::iterator it = bonds.begin(); it != bonds.end(); it++) {
                if ((nbrinfo[it->p1()].degree() == c1 && nbrinfo[it->p2()].degree() == c2) ||
                    (nbrinfo[it->p2()].degree() == c1 && nbrinfo[it->p1()].degree() == c2)) {
                    c12bonds.push_back(ordered_key(it->p1(), it->p2()));
                }
            }
            
            while (c12bonds.size() > 0) {
                // pick a bond randomly, remove it
                int randpos = randomIndex(c12bonds.size());
                ConnSet::iterator it2 = c12bonds.begin();
                advance(it2, randpos);
                
                int pa = it2->first; int pb = it2->second;
                removeBond(pa, pb);
                if (conn() <= conn_target) {
                    printf(" new connectivity: %1.4f no. of bonds: %lu\n", conn(), bonds.size());
                    return true;
                }
                
                // remove from the set of pairs all pairs that involved the points in the random bond
                for (ConnSet::iterator it3 = c12bonds.begin(); it3 != c12bonds.end(); ) {
                    if (it3->first == pa || it3->second == pa) c12bonds.erase(it3++);
                    else if (it3->first == pb || it3->second == pb) c12bonds.erase(it3++);
                    else it3++;
                }
            }
        }
    }

    if (conn() > conn_target) return false;
    else return true;
}

void PointSystem2d::random_force_dipoles(double f, double p) {
    // introduce force dipoles to a fraction p of the bonds with force f
    for (int i = 0; i <= bonds.size()*p; i++) {
        Bond bond_to_squeeze = bondFromIdx(randomBondIndex());
        force_dipoles.push_back(Bond(bond_to_squeeze.p1(), bond_to_squeeze.p2(), f)); 
    }
}


//~ // volume and area
//~ double PointSystem2d::elementsurf(Facet f) {
    //~ Point2d p1 = points[f.p1()];
    //~ Point2d p2 = points[f.p2()];
    //~ Point2d p3 = points[f.p3()];
    //~ 
    //~ Point2d v12 = p2 - p1; //really vectors
    //~ Point2d v13 = p3 - p1;
    //~ 
    //~ Point2d vcross = v12 % v13;
    //~ 
    //~ return sqrt(vcross*vcross)/2;
//~ }
//~ 
//~ 
//~ double PointSystem2d::elementvol(Facet f) {
    //~ Point2d p1 = points[f.p1()];
    //~ Point2d p2 = points[f.p2()];
    //~ Point2d p3 = points[f.p3()];
    //~ 
    //~ return (p1*(p2%p3))/6;
//~ }
//~ 
//~ double PointSystem2d::surf() {
    //~ double sa = 0;
    //~ for (FacetArray::iterator it = facets.begin(); it != facets.end(); it++) {
        //~ sa += elementsurf(*it);
    //~ }
    //~ return sa;
//~ }
//~ 
//~ double PointSystem2d::vol() {
    //~ double v = 0;
    //~ for (FacetArray::iterator it = facets.begin(); it != facets.end(); it++) {
        //~ v += elementvol(*it);
    //~ }
    //~ return v;
//~ }
//~ 
//~ double PointSystem2d::epressure() {
    //~ return pressure*vol();
//~ }
//~ 
//~ void PointSystem2d::d_epressure(PointArray* grdvec) {
    //~ for (FacetArray::iterator it = facets.begin(); it != facets.end(); it++) {
        //~ Point2d p1 = points[it->p1()];
        //~ Point2d p2 = points[it->p2()];
        //~ Point2d p3 = points[it->p3()];
        //~ 
        //~ (*grdvec)[it->p1()] = (*grdvec)[it->p1()]+(p2 % p3)*pressure/6;
        //~ (*grdvec)[it->p2()] = (*grdvec)[it->p2()]+(p3 % p1)*pressure/6;
        //~ (*grdvec)[it->p3()] = (*grdvec)[it->p3()]+(p1 % p2)*pressure/6;
    //~ }
//~ }


// harmonic spring
double PointSystem2d::espring(Bond b) {
    double length = bondlength(b);
    return 0.5*b.k()*(length-b.l0())*(length-b.l0());
}


double PointSystem2d::estretch() {
    double e_stretch = 0;
    for (BondArray::iterator it = bonds.begin(); it != bonds.end(); it++) {
        e_stretch += espring(*it);
    }
    return e_stretch;
}

double PointSystem2d::eforcedipole() {
    double e_fd = 0;
    for (BondArray::iterator it = force_dipoles.begin(); it != force_dipoles.end(); it++) {
        e_fd += it->k()*bondlength(*it); // use k() to store the force magnitude
    }
    return e_fd;
}

Point2d PointSystem2d::fstretch(int p1, int p2) {
    // return the vector force on p1 due to p2, 0 if no bond
    PtPair bondkey = ordered_key(p1,p2);
    if (bondpointers.count(bondkey) == 0) return 0.;
    
    BondArray::iterator it = bondpointers[bondkey];
    Point2d bd = points[p2]-points[p1];
    double r = sqrt(bd*bd);
    
    return bd*(it->k()*(r - it->l0())/r);
}

double PointSystem2d::fmag(int p1, int p2) {
    // return the scalar force on p1 due to p2, 0 if no bond
    PtPair bondkey = ordered_key(p1,p2);
    if (bondpointers.count(bondkey) == 0) return 0.;
    
    BondArray::iterator it = bondpointers[bondkey];
    Point2d bd = points[p2]-points[p1];
    double r = sqrt(bd*bd);
    
    return it->k()*(r - it->l0());
}

double PointSystem2d::fmag(Bond b) {
    // return the scalar force on p1 due to p2, 0 if no bond
    double r = bondlength(b);
    return b.k()*(r - b.l0());
}


void PointSystem2d::d_estretch(PointArray* grdvec) {
    for (BondArray::iterator it = bonds.begin(); it != bonds.end(); it++) {
        int p1 = it->p1();
        int p2 = it->p2();
        Point2d bd = points[p1]-points[p2];
        double r = sqrt(bd*bd);
        (*grdvec)[p1] = (*grdvec)[p1]+bd*(it->k()*(r - it->l0())/r);
        (*grdvec)[p2] = (*grdvec)[p2]-bd*(it->k()*(r - it->l0())/r);
    }
}

void PointSystem2d::d_eforcedipole(PointArray* grdvec) {
    for (BondArray::iterator it = force_dipoles.begin(); it != force_dipoles.end(); it++) {
        int p1 = it->p1();
        int p2 = it->p2();
        Point2d bd = points[p1]-points[p2];
        double r = sqrt(bd*bd);
        (*grdvec)[p1] = (*grdvec)[p1]+bd*(it->k()/r);
        (*grdvec)[p2] = (*grdvec)[p2]-bd*(it->k()/r);
    }
}

//~ double PointSystem2d::emorse() {
    //~ double e_stretch = 0;
    //~ double ea;
    //~ for (BondArray::iterator it = bonds.begin(); it != bonds.end(); it++) {
        //~ ea = exp(-it->k()*(bondlength(*it)-it->l0()));
        //~ e_stretch += ea*ea - 2*ea;
    //~ }
    //~ return e_stretch;
//~ }
//~ 
//~ void PointSystem2d::d_emorse(PointArray* grdvec) {
    //~ for (BondArray::iterator it = bonds.begin(); it != bonds.end(); it++) {
        //~ int p1 = it->p1();
        //~ int p2 = it->p2();
        //~ Point2d bd = points[p1]-points[p2];
        //~ double r = sqrt(bd*bd);
        //~ double ea = exp(-it->k()*(r-it->l0()));
        //~ (*grdvec)[p1] = (*grdvec)[p1]+bd*(2*it->k()*(-ea*ea + ea)/r);
        //~ (*grdvec)[p2] = (*grdvec)[p2]-bd*(2*it->k()*(-ea*ea + ea)/r);
    //~ }
//~ }
//~ 
//~ double PointSystem2d::esphere() {
    //~ double en = 0;
    //~ double delta;
    //~ for (int i =0; i < N; i++) {
        //~ delta = sqrt(points[i]*points[i])-sphr;
        //~ en += delta*delta/2.;
    //~ }
    //~ return en*sphpot;
//~ }
//~ 
//~ void PointSystem2d::d_esphere(PointArray* grdvec) {
    //~ for (int i = 0; i < N; i++) {
        //~ (*grdvec)[i] = (*grdvec)[i] + points[i]*((1.-sphr/sqrt(points[i]*points[i]))*sphpot);
    //~ }
//~ }
//~ 
double PointSystem2d::erefpts() {
    double e_ref = 0;
    for (BondArray::iterator it = refbonds.begin(); it != refbonds.end(); it++) {
        int pmesh = it->p1();
        int pref = it->p2();
        Point2d bd = points[pmesh]-refpoints[pref];
        
        if (it->k() > 0) { // x-y restrict only
            bd.set(bd.x(), bd.y(), lx, ly, shearx); 

            double r = sqrt(bd*bd);
            //~ // e_ref += 0.5*(r - it->l0())*(r - it->l0())*it->k();
            e_ref += 0.5*r*r*it->k();
        }
        
        //~ if (it->k() < 0) { // z restrict only
            //~ bd.set(0, 0, bd.z()); 
//~ 
            //~ double r = sqrt(bd*bd);
            // e_ref += 0.5*(r - it->l0())*(r - it->l0())*it->k();
            //~ e_ref += -0.5*r*r*it->k(); // correct for setting k < 0
        //~ }
    }
    return e_ref;
}

void PointSystem2d::d_erefpts(PointArray* grdvec) {
    for (BondArray::iterator it = refbonds.begin(); it != refbonds.end(); it++) {
        int pmesh = it->p1();
        int pref = it->p2();
        Point2d bd = points[pmesh]-refpoints[pref];
        
        if (it->k() > 0) {
            // x-y restrict only
            bd.set(bd.x(), bd.y(), lx, ly, shearx);

            //~ double r = sqrt(bd*bd);
            //~ // (*grdvec)[pmesh] = (*grdvec)[pmesh]+bd*(it->k()*(r - it->l0())/r); // only Point2d belonging to mesh has gradient updated
            (*grdvec)[pmesh] = (*grdvec)[pmesh]+bd*it->k(); // only Point2d belonging to mesh has gradient updated
            //~ // (*grdvec)[p2] = (*grdvec)[p2]-bd*(it->k()*(r - it->l0())/r);
        }
        //~ if (it->k() < 0) {
            //~ // x-y restrict only
            //~ bd.set(0, 0, bd.z());
//~ 
            //~ double r = sqrt(bd*bd);
            // (*grdvec)[pmesh] = (*grdvec)[pmesh]+bd*(it->k()*(r - it->l0())/r); // only Point2d belonging to mesh has gradient updated
            //~ (*grdvec)[pmesh] = (*grdvec)[pmesh] - bd*it->k(); // compensate for k < 0
            // (*grdvec)[p2] = (*grdvec)[p2]-bd*(it->k()*(r - it->l0())/r);
        //~ }
    }
}
//~ 
//~ // tie Point2d idx to a reference Point2d (x,y,z) with a spring of constant k
void PointSystem2d::makeref(int idx, double x, double y,  double k) {
    refpoints.push_back(Point2d(x,y,lx,ly,shearx));
    int pref = refpoints.size()-1;
    double l0 = points[idx].distanceTo(refpoints[pref]);
    refbonds.push_back(Bond(idx, pref, k, l0));
}

void PointSystem2d::makeref_self(int idx, double k) {
    makeref(idx, points[idx].x(), points[idx].y(), k);
}

double PointSystem2d::energy() {
    double en = 0;
    //~ if (morse) en = emorse();
    //~ else en = estretch();
    en = estretch();
    //~ en += ebend(); 
    //~ if (gamma > 0) en += esurf();
    //~ if (lambda > 0) en += evol();
    //~ if (pressure != 0) en += epressure();
    //~ en += esphere();
    en += erefpts();
    en += eforcedipole();
    return en;
}

void PointSystem2d::d_energy(PointArray* grd) {
    //~ if (morse) d_emorse(grd);
    d_estretch(grd);
    //~ else d_estretch(grd);
    //~ d_ebend(grd);
    //~ if (gamma > 0) d_esurf(grd);
    //~ if (lambda > 0) d_evol(grd);
    //~ if (pressure != 0) d_epressure(grd);
    //~ d_esphere(grd);
    d_erefpts(grd);
    d_eforcedipole(grd);
}


vector<double> PointSystem2d::stress_tensor() {
    static const double arr[] = {0,0,0,0};
    vector<double> res (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    Point2d fij;
    Point2d rji;
    for (BondArray::iterator it = bonds.begin(); it != bonds.end(); it++) {
        fij = fstretch(it->p1(), it->p2());
        rji = points[it->p2()]-points[it->p1()];
        res[0] += fij.x()*rji.x();
        res[1] += fij.x()*rji.y();
        res[2] += fij.y()*rji.x();
        res[3] += fij.y()*rji.y();
    }
    return res;
}

// apply random shake of size <= amplitude in all directions
void PointSystem2d::shake(double amplitude) {
    for (int idx = 0; idx < N; idx++) {
        points[idx].set(points[idx] + Point2d((double) rand()/RAND_MAX-.5, (double) rand()/RAND_MAX-.5)*amplitude);
    }
}


// i/o functions
int PointSystem2d::write_points(string filename) {
    FILE* op = fopen(filename.c_str(),"w");
    for (PointArray::iterator it = points.begin(); it != points.end(); it++) {
        fprintf(op, "%2.8f %2.8f\n", it->x(), it->y());
    }
    fclose(op);
    return 0;
}

void PointSystem2d::write_bonds(string filename) {
    FILE* op = fopen(filename.c_str(),"w");
    for (BondArray::iterator it = bonds.begin(); it != bonds.end(); it++) {
        fprintf(op, "%d %d\n", it->p1()+ZEROCORR, it->p2()+ZEROCORR);
    }
    fclose(op);
}

void PointSystem2d::write_data(string filename) {
    // write points, bonds, facets with .pts, .bds, .fct extensions
    write_points(filename + ".pts");
    write_bonds(filename + ".bds");
    //~ write_facets(filename + ".fct");
}

void PointSystem2d::print_energies(string prefix = "") {
    printf(" %s\te_stretch\t%2.8e\n",prefix.c_str(), estretch());
    //~ printf(" %s\te_morse\t%2.8e\n",prefix.c_str(), emorse());
    //~ printf(" %s\te_clamp\t%2.8e\n",prefix.c_str(), esphere());
    //~ printf(" %s\tvol\t%2.8f\n",prefix.c_str(), vol());
    //~ printf(" %s\te_pres\t%2.8e\n",prefix.c_str(), epressure());
    //~ printf(" %s\te_ref\t%2.8e\n",prefix.c_str(), erefpts());
    printf(" %s\te_tot\t%2.8e\n",prefix.c_str(), energy());
}

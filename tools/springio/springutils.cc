////////////////////////////////////
// Utilities for 2D spring networks     //
// Uses PointSystem2d class; dependent on Boost libraries
////////////////////////////////////
#include "PointSystem2d/PointSystem2d.h"
#include "PointSystem2d/PointSystem2d.cpp"
#include "PointSystem2d/Point2d.h"
#include "springio.cc"


void read_points(string ptfile, PointSystem2d &pts, ldf boxsize) {
    /* Read two-dimensional point data from ptfile into 'PointSystem2d' structure.
     * Each row of ptfile must contain two entries, corresponding to 'x' and 'y' coordinates */
	pts.lx = boxsize;
	pts.ly = boxsize;
	
    FILE* inputM;
    double xin, yin;
    
    inputM = fopen(ptfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%lf %lf\n", &xin, &yin);
        pts.addPoint(Point2d(xin,yin,boxsize, boxsize,0));
    }
}

void read_bonds(string bfile, PointSystem2d &pts) {
    /* Read in connectivity data from bfile into PointSystem2d structure.
     * each row of bfile contains five entries: idx1 idx2 bondtype springconstant restlength */
    ui p1in, p2in, dummy;
    ldf kin, l0in;
    FILE* inputM = fopen(bfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%d %d %d %Lf %Lf\n", &p1in, &p2in, &dummy, &kin, &l0in);
        // spring with k and r0
        pts.addBond(p1in-INDEXSHIFT, p2in-INDEXSHIFT, kin, l0in);
    }
}


void ps2md(PointSystem2d &pts, md<2> &sys) {
    /* transfer data from PointSystem2d structure to md structure */
	ldf x[pts.N];
	ldf y[pts.N];
	
	vector<double> xv(0);
    vector<double> yv(0);
	
	// copy over points
	for (int i = 0; i < pts.N; i++) {
		xv.push_back(pts.points[i].x()); yv.push_back(pts.points[i].y());
	}
	copy(xv.begin(), xv.end(), x);
    copy(yv.begin(), yv.end(), y);
	
	sys.import_pos(x,y);
	
	// set interactions
	for (int i = 0; i < pts.N; i++) sys.set_type(i,i); // not essential but makes add_spring() work less.
	for (BondArray::iterator it = pts.bonds.begin(); it != pts.bonds.end(); it++) {
        sys.add_spring(it->p1(), it->p2(),it->k(),it->l0());
	}
}

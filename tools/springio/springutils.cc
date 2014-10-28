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
        dummy = fscanf(inputM, "%lf %lf\n", &xin, &yin);
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
        dummy = fscanf(inputM, "%d %d %d " F_LDF " " F_LDF "\n", &p1in, &p2in, &dummy, &kin, &l0in);
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

vector<double> stress_tensor(md<2> &sys) {
    static const double arr[] = {0.,0.,0.,0.};
    vector<double> res (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    ldf rcosq = pow(sys.network.rco,2);
    for (ui i = 0; i < sys.N; i++) {
        // stress contribution due to central forces
        for(ui j=sys.network.skins[i].size()-1;j<UI_MAX;j--) if(i>sys.network.skins[i][j].neighbor)
        {   
            const ldf rsq=sys.distsq(i,sys.network.skins[i][j].neighbor);
            if(!sys.network.update or rsq<rcosq)
            {   
                vector<double> fij(2,0.);
                vector<double> rji(2);
                const ldf r=sqrt(rsq);
                const ldf dVdr=sys.v.dr(sys.network.library[sys.network.skins[i][j].interaction].potential,r,&sys.network.library[sys.network.skins[i][j].interaction].parameters);
                for(ui d=0;d<2;d++) {
                    rji[d] = sys.dd(d,i,sys.network.skins[i][j].neighbor);
                    fij[d] = rji[d]*dVdr/r;
                }
                res[0] += fij[0]*rji[0];
                res[1] += fij[0]*rji[1];
                res[2] += fij[1]*rji[0];
                res[3] += fij[1]*rji[1];
            }
        }
        // stress contribution due to dissipative springs. NOTE: does not consider all pair forces, only the DISSIPATION pair force
        if(sys.network.forcelibrary.size() and sys.network.forces[i].size()) for(ui k=sys.network.forces[i].size()-1;k<UI_MAX;k--)
        {
            ui ftype=sys.network.forces[i][k];
            if(sys.network.forcelibrary[ftype].externalforce == EXTFORCE::DISSIPATION and sys.network.forcelibrary[ftype].particles.size() and sys.network.forcelibrary[ftype].particles[i].size())
            {
                vector<ui> plist = sys.network.forcelibrary[ftype].particles[i];
                ldf b = sys.network.forcelibrary[ftype].parameters[0];
                for (vector<ui>::iterator it = plist.begin(); it != plist.end(); it++) {
                    ui j = *it;
                    vector<double> fij(2,0.);
                    vector<double> rji(2);
                    for(ui d=0;d<2;d++) {
                        rji[d] = sys.dd(d,i,j);
                        fij[d] = b*sys.dv(d,i,j);
                    }
                    res[0] += fij[0]*rji[0];
                    res[1] += fij[0]*rji[1];
                    res[2] += fij[1]*rji[0];
                    res[3] += fij[1]*rji[1];
                }
            }
        }
    }
    return res;
}


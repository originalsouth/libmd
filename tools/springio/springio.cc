////////////////////////////////////
// I/O for 2D spring networks     //
// Convenience functions for      // 
// interfacing with files such as //
// Stephan's MD input/output      //
////////////////////////////////////
#include "PointSystem2d/PointSystem2d.h"
#include "PointSystem2d/Point2d.h"


using namespace std;

#define INDEXSHIFT 0 // set to 1 for Stephan's conn.txt files, for which point indexing starts from 1

ui number_of_lines(string ptfile) {
    int nl = 0;
    string line;
    ifstream myfile(ptfile);

    while (getline(myfile, line))
        ++nl;
    return nl;
}

void read_points_ulrich(string ptfile, PointSystem2d &pts, ldf boxsize) {
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

void read_points_ulrich(string ptfile, ldf *x, ldf* y) {
    FILE* inputM;
    double xin, yin;
    
    vector<double> xv(0);
    vector<double> yv(0);
    inputM = fopen(ptfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%lf %lf\n", &xin, &yin);
        xv.push_back(xin); yv.push_back(yin);
    }
    
    copy(xv.begin(), xv.end(), x);
    copy(yv.begin(), yv.end(), y);
}

void read_bonds_ulrich(string bfile, PointSystem2d &pts) {
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
	ldf x[pts.N];
	ldf y[pts.N];
	
	vector<double> xv(0);
    vector<double> yv(0);
	
	// copy over points
	for (ui i = 0; i < pts.N; i++) {
		xv.push_back(pts.points[i].x()); yv.push_back(pts.points[i].y());
	}
	copy(xv.begin(), xv.end(), x);
    copy(yv.begin(), yv.end(), y);
	
	sys.import_pos(x,y);
	
	// set interactions
	for (ui i = 0; i < pts.N; i++) sys.particles[i].type = i;
	for (BondArray::iterator it = pts.bonds.begin(); it != pts.bonds.end(); it++) {
		vector<ldf> a = {it->k(), it->l0()};
		sys.add_typeinteraction(it->p1(), it->p2(), 2, &a);
	}
}

template<ui dim> void read_bonds_ulrich(string bfile, md<dim> &sys) {
    ui p1in, p2in, dummy;
    ldf kin, l0in;
    FILE* inputM = fopen(bfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%d %d %d %Lf %Lf\n", &p1in, &p2in, &dummy, &kin, &l0in);
        // spring with k and r0
        sys.add_spring(p1in-INDEXSHIFT, p2in-INDEXSHIFT,kin,l0in);
    }
}

template<ui dim> void write_points_x(string filename, md<dim> &sys) {
    /* write N*dim array of point positions */
    FILE* op = fopen(filename.c_str(),"w");
    for (int i = 0; i < sys.N; i++) {
        for (int d = 0; d < dim; d++) {
            fprintf(op, "%2.8Lf ", sys.particles[i].x[d]);  
        }
        fprintf (op, "\n");
    }
    fclose(op);
}

template<ui dim> void write_points_v(string filename, md<dim> &sys) {
    /* write N*dim array of point positions */
    FILE* op = fopen(filename.c_str(),"w");
    for (int i = 0; i < sys.N; i++) {
        for (int d = 0; d < dim; d++) {
            fprintf(op, "%2.8Lf ", sys.particles[i].dx[d]);
        }
        fprintf (op, "\n");
    }
    fclose(op);
}

template<ui dim> void write_points_f(string filename, md<dim> &sys) {
    /* write N*dim array of point positions */
    FILE* op = fopen(filename.c_str(),"w");
    for (int i = 0; i < sys.N; i++) {
        for (int d = 0; d < dim; d++) {
            fprintf(op, "%2.8Lf ", sys.particles[i].F[d]);
        }
        fprintf (op, "\n");
    }
    fclose(op);
}

template<ui dim> void write_bonds(string filename, md<dim> &sys) {
	// TODO: add interface functions to md class so that user does not touch the inner workings
    FILE* op = fopen(filename.c_str(),"w");
    for (ui i = 0; i < sys.N; i++) {
        for(ui j=sys.network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>sys.network.skins[i][j].neighbor) {
            fprintf(op, "%d %d\n", i, sys.network.skins[i][j].neighbor);
        }
    }
    fclose(op);
}

template<ui dim> void write_data(string prefix, md<dim> &sys) {
	write_points_x(prefix+".pts",sys);
	write_points_v(prefix+".vel",sys);
	write_points_f(prefix+".f",sys);
	write_bonds(prefix+".bds",sys);
}


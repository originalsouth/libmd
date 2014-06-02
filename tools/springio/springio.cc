////////////////////////////////////
// I/O for 2D spring networks     //
// Convenience functions for      // 
// interfacing with files such as //
// Stephan's MD input/output      //
////////////////////////////////////
#include <iostream>
#include <fstream>

using namespace std;

#define INDEXSHIFT 0 // set to 1 for Stephan's conn.txt files, for which point indexing starts from 1
int dummy;  // swallows integer returned by fscanf

ui number_of_lines(string ptfile) {
    int nl = 0;
    string line;
    ifstream myfile(ptfile);

    while (getline(myfile, line))
        ++nl;
    return nl;
}

void read_points(string ptfile, ldf *x, ldf* y) {
    /* Read two-dimensional point data from ptfile into 'x' and 'y' arrays 
     * Each row of ptfile must contain two entries, corresponding to 'x' and 'y' coordinates */
    FILE* inputM;
    double xin, yin;
    
    vector<double> xv(0);
    vector<double> yv(0);
    inputM = fopen(ptfile.c_str(), "r");
    while (!(feof(inputM))) {
        dummy = fscanf(inputM, "%lf %lf\n", &xin, &yin);
        xv.push_back(xin); yv.push_back(yin);
    }
    
    copy(xv.begin(), xv.end(), x);
    copy(yv.begin(), yv.end(), y);
}

template<ui dim> void read_bonds(string bfile, md<dim> &sys) {
    /* Read in connectivity data from bfile into md structure, creating harmonic bonds.
     * each row of bfile contains five entries: idx1 idx2 bondtype springconstant restlength */
    ui p1in, p2in, dummy;
    ldf kin, l0in;
    FILE* inputM = fopen(bfile.c_str(), "r");
    while (!(feof(inputM))) {
        dummy = fscanf(inputM, "%d %d %d %Lf %Lf\n", &p1in, &p2in, &dummy, &kin, &l0in);
        // spring with k and r0
        sys.add_spring(p1in-INDEXSHIFT, p2in-INDEXSHIFT,kin,l0in);
    }
}

template<ui dim> void read_bonds(string bfile, md<dim> &sys,vector<vector<ui>> &nbrlist, ldf kfactor=1.) {
    /* Read neighbors into a list of neighbor indices, in addition to adding springs to md structure */
    ui p1in, p2in, dummy;
    ldf kin, l0in;
    FILE* inputM = fopen(bfile.c_str(), "r");
    while (!(feof(inputM))) {
        dummy = fscanf(inputM, "%d %d %d %Lf %Lf\n", &p1in, &p2in, &dummy, &kin, &l0in);
        // spring with k and r0
        sys.add_spring(p1in-INDEXSHIFT, p2in-INDEXSHIFT,kin*kfactor,l0in);
        // update nbrlist
        nbrlist[p1in-INDEXSHIFT].push_back(p2in-INDEXSHIFT);
        nbrlist[p2in-INDEXSHIFT].push_back(p1in-INDEXSHIFT);
    }
}

template<ui dim> void write_points_x(string filename, md<dim> &sys) {
    /* write N*dim array of point positions */
    FILE* op = fopen(filename.c_str(),"w");
    for (unsigned int i = 0; i < sys.N; i++) {
        for (unsigned int d = 0; d < dim; d++) {
            fprintf(op, "%2.8Lf ", sys.particles[i].x[d]);  
        }
        fprintf (op, "\n");
    }
    fclose(op);
}

template<ui dim> void write_points_v(string filename, md<dim> &sys) {
    /* write N*dim array of point velocities */
    FILE* op = fopen(filename.c_str(),"w");
    for (unsigned int i = 0; i < sys.N; i++) {
        for (unsigned int d = 0; d < dim; d++) {
            fprintf(op, "%2.8Lf ", sys.particles[i].dx[d]);
        }
        fprintf (op, "\n");
    }
    fclose(op);
}

template<ui dim> void write_points_f(string filename, md<dim> &sys) {
    /* write N*dim array of forces */
    FILE* op = fopen(filename.c_str(),"w");
    for (unsigned int i = 0; i < sys.N; i++) {
        for (unsigned int d = 0; d < dim; d++) {
            fprintf(op, "%2.8Lf ", sys.particles[i].F[d]);
        }
        fprintf (op, "\n");
    }
    fclose(op);
}

template<ui dim> void write_bonds(string filename, md<dim> &sys) {
	/* write connectivity data from the skins structures. */
    FILE* op = fopen(filename.c_str(),"w");
    for (unsigned int i = 0; i < sys.N; i++) {
        for(unsigned int j=sys.network.skins[i].size()-1;j<UI_MAX;j--) if(i>sys.network.skins[i][j].neighbor) {
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


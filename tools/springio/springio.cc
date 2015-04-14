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

string dataprecision = "%2.8"; // precision for string datafile output
string floatformat = string(F_LDF).erase(0,1); // strip leading '%' from F_LDF format string
string formatstring_ = dataprecision + floatformat + " ";

const char* formatstring = formatstring_.c_str(); // format string  for datafile output

ui number_of_lines(string ptfile) {
    int nl = 0;
    string line;
    ifstream myfile(ptfile);

    while (getline(myfile, line))
        ++nl;
    return nl;
}

int read_points(string ptfile, ldf *x, ldf* y) {
    /* Read two-dimensional point data from ptfile into 'x' and 'y' arrays 
     * Each row of ptfile must contain two entries, corresponding to 'x' and 'y' coordinates */
    FILE* inputM;
    ldf xin, yin;
    
    vector<ldf> xv(0);
    vector<ldf> yv(0);
    inputM = fopen(ptfile.c_str(), "r");
    if(inputM!=NULL) {
    while (!(feof(inputM))) {
        dummy = fscanf(inputM, F_LDF " " F_LDF "\n", &xin, &yin);
        xv.push_back(xin); yv.push_back(yin);
    }
    
    copy(xv.begin(), xv.end(), x);
    copy(yv.begin(), yv.end(), y);
    fclose(inputM);
    return 0;
    } else {
        cout << "Couldn't open " << ptfile << endl;
        return 1;
    }
}


int read_points(string ptfile, ldf *x, ldf* y, ldf* z) {
    /* Read three-dimensional point data from ptfile into 'x', 'y', 'z' arrays 
     * Each row of ptfile must contain three entries */
    FILE* inputM;
    ldf xin, yin, zin;
    
    vector<ldf> xv(0);
    vector<ldf> yv(0);
    vector<ldf> zv(0);
    inputM = fopen(ptfile.c_str(), "r");
    if(inputM!=NULL) {
    while (!(feof(inputM))) {
        dummy = fscanf(inputM, F_LDF " " F_LDF " " F_LDF "\n", &xin, &yin, &zin);
        xv.push_back(xin); yv.push_back(yin); zv.push_back(zin);
    }
    
    copy(xv.begin(), xv.end(), x);
    copy(yv.begin(), yv.end(), y);
    copy(zv.begin(), zv.end(), z);
    fclose(inputM);
    return 0;
    } else {
        cout << "Couldn't open " << ptfile << endl;
        return 1;
    }
}



template<ui dim> int read_bonds(string bfile, md<dim> &sys) {
    /* Read in connectivity data from bfile into md structure, creating harmonic bonds.
     * each row of bfile contains five entries: idx1 idx2 bondtype springconstant restlength */
    ui p1in, p2in, dummy;
    ldf kin, l0in;
    FILE* inputM = fopen(bfile.c_str(), "r");
    if(inputM!=NULL) {
    while (!(feof(inputM))) {
        dummy = fscanf(inputM, "%d %d %d " F_LDF " " F_LDF "\n", &p1in, &p2in, &dummy, &kin, &l0in);
        // spring with k and r0
        sys.add_spring(p1in-INDEXSHIFT, p2in-INDEXSHIFT,kin,l0in);
    }
    fclose(inputM);
    return 0;
    } else {
        cout << "Couldn't open " << bfile << endl;
        return 1;
    }
}

template<ui dim> int read_bonds(string bfile, md<dim> &sys,vector<vector<ui>> &nbrlist, ldf kfactor=1.) {
    /* Read neighbors into a list of neighbor indices, in addition to adding springs to md structure */
    ui p1in, p2in, dummy;
    ldf kin, l0in;
    FILE* inputM = fopen(bfile.c_str(), "r");
    if(inputM!=NULL) {
    while (!(feof(inputM))) {
        dummy = fscanf(inputM, "%d %d %d " F_LDF " " F_LDF "\n", &p1in, &p2in, &dummy, &kin, &l0in);
        // spring with k and r0
        sys.add_spring(p1in-INDEXSHIFT, p2in-INDEXSHIFT,kin*kfactor,l0in);
        // update nbrlist
        nbrlist[p1in-INDEXSHIFT].push_back(p2in-INDEXSHIFT);
        nbrlist[p2in-INDEXSHIFT].push_back(p1in-INDEXSHIFT);
    }
    fclose(inputM);
    return 0;
    } else {
        cout << "Couldn't open " << bfile << endl;
        return 1;
    }
}

template<ui dim> void write_points_x(string filename, md<dim> &sys) {
    /* write N*dim array of point positions */
    FILE* op = fopen(filename.c_str(),"w");
    for (unsigned int i = 0; i < sys.N; i++) {
        for (unsigned int d = 0; d < dim; d++) {
            fprintf(op, formatstring, sys.particles[i].x[d]);  
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
            fprintf(op, formatstring, sys.particles[i].dx[d]);
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
            fprintf(op, formatstring, sys.particles[i].F[d]);
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

template<ui dim> void write_energies(string filename, md<dim> &sys) {
    /* write energies (H, T, V) to filename : use scientific notation */
    FILE* op = fopen(filename.c_str(),"w");
    fprintf(op, formatstring, sys.H());
    fprintf(op, formatstring, sys.T());
    fprintf(op, formatstring, sys.V());
    fprintf (op, "\n");
    fclose(op);
}

template<ui dim> void append_energies(string filename, md<dim> &sys) {
    /* append(!) energies (H, T, V) to filename : use scientific notation */
    FILE* op = fopen(filename.c_str(),"a");
    fprintf(op, formatstring, sys.H());
    fprintf(op, formatstring, sys.T());
    fprintf(op, formatstring, sys.V());
    fprintf (op, "\n");
    fclose(op);
}

template<ui dim> void write_data(string prefix, md<dim> &sys) {
	write_points_x(prefix+".pts",sys);
	write_points_v(prefix+".vel",sys);
	write_points_f(prefix+".f",sys);
	write_bonds(prefix+".bds",sys);
}


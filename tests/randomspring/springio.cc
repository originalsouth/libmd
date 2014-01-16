////////////////////////////////////
// I/O for 2D spring networks     //
// Convenience functions for      // 
// interfacing with files such as //
// Stephan's MD input/output      //
////////////////////////////////////

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

template<ui dim> void read_bonds_ulrich(string bfile, md<dim> &sys) {
    ui p1in, p2in, dummy;
    ldf kin, l0in;
    FILE* inputM = fopen(bfile.c_str(), "r");
    while (!(feof(inputM))) {
        fscanf(inputM, "%d %d %d %Lf %Lf\n", &p1in, &p2in, &dummy, &kin, &l0in);
        // spring with k and r0
        vector<ldf> a={kin,l0in};
        sys.add_typeinteraction(p1in-INDEXSHIFT, p2in-INDEXSHIFT, 2, &a);
    }
}

template<ui dim> void write_points(string filename, md<dim> &sys) {
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

template<ui dim> void write_bonds(string filename, md<dim> &sys) {
    //~ FILE* op = fopen(filename.c_str(),"w");
    //~ for (int i = 0; i < sys.N) {
        //~ for (int d = 0; d < dim; d++) {
            //~ fprintf(op, "%2.8f ", sys.particles[i].x[d]);
        //~ }
        //~ fprintf (op, "\n");
    //~ }
    //~ fclose(op);
}

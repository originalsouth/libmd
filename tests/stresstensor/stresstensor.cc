/***************************************************************************************************
 *   stresstensor.cc
 * 
 * Tests the stress tensor measurement from tools/springio, for harmonic springs,
 * and matches to the theoretical expression for the shear and bulk moduli.
 * 
 * Usage: stresstensor <shear strain> <bulk strain>
 * Output: stress tensor components as a function of time; should converge to the expected stresses
 * for a random harmonic spring network at long times.
 * Program runs until terminated by user.
 * 
 * This program does *not* test the dissipative component of the stress tensor.
 * 
 ***************************************************************************************************/

#include "../../libmd.cc"
#include "../../tools/springio/springutils.cc"

using namespace std;

// Initial configuration of points and bonds; random six-fold coordinated network with springs of spring constant 1
string ptfile = "smd.pts";
string bfile = "smd.bds";
ldf relaxedboxsize = 45.5;  // box size for which the initial configuration is relaxed under periodic BC
ldf springconst = 1.; // spring constant (column 4 in smd.bds)

int main(int argc, char* argv[])
{   
    // Parse command-line arguments
    if (argc < 3) { cout <<" Usage: stresstensor <shear strain> <bulk strain>" << endl; return EXIT_FAILURE; }
    ldf gamma = atof(argv[1]); // shear strain applied to boundaries perpendicular to x
    ldf eta = atof(argv[2]);   // bulk strain applied as a uniform dilation of the spring network
    
    // initialize points
    // read points from ptfile into arrays
    ui systemsize = number_of_lines(ptfile);
    ldf x[systemsize];
    ldf y[systemsize];
    read_points(ptfile,x,y);
    
    // make md system
    md<2> sys(systemsize);
    sys.parallel.set(1);
    sys.set_rco(15.);
    sys.set_ssz(15.);
    
    
    ldf boxsize = relaxedboxsize*(1+eta);      // Induce uniform dilation by increasing box size
    sys.simbox.L[0]=boxsize; 
    sys.simbox.L[1]=boxsize;
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.integrator.method=INTEGRATOR::SEULER;
    sys.integrator.h = 0.001;
    
    // import the point positions read from ptfile
    sys.import_pos(x,y);

    
    // initialize bonds directly from bfile
    vector<vector<ui>> springnbrs(sys.N);
    read_bonds(bfile,sys,springnbrs);

    // index once for a static spring network
    sys.indexdata.method = INDEX::CELL;
    sys.index();
    sys.network.update=false;
    
    // add a uniform drag to all particles; this lets the springs relax to a static solution for the applied strains
    sys.set_damping(.1);
    
    // shear the box until the required shear is obtained. Do this over some time steps, so that particles have a chance to relax
    sys.simbox.shear_boundary(1,0,gamma*relaxedboxsize/(100*sys.integrator.h)); // This sets the shear velocity so that the box has sheared by gamma*L in 100 time steps
    sys.timesteps(100);
    sys.simbox.shear_boundary(1,0,0.);  // turn of dynamic shear; box is now frozen in static shear condition
    
    // theoretical values of sxx, sxy and syy for a 6-coordinated spring network
    ldf Y = 2.*springconst/sqrt(3.);  // 2D Young's modulus
    ldf nu = 1./3.;                   // Poisson ratio
    ldf sxy = gamma*Y/(2*(1+nu));     // sigma_xy
    ldf sxx = Y*eta/(1-nu);           // also equal to sigma_yy

    for(ui h=0;h<UI_MAX;h++)
    {
        sys.timesteps(1000);
        cout << h <<" ";
        vector<double> sij = stress_tensor(sys);
        cout << "\tsij: " 
             << sij[0]/(boxsize*boxsize) << " "<< sij[1]/(boxsize*boxsize) << " "
             << sij[2]/(boxsize*boxsize) << " "<< sij[3]/(boxsize*boxsize) 
             << " \tsim/theory: " 
             << sij[0]/(boxsize*boxsize)/sxx << " "<< sij[1]/(boxsize*boxsize)/sxy << " "<< sij[3]/(boxsize*boxsize)/sxx<<endl; 
        
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

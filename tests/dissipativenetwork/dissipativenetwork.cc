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
#include "../../tools/BaX/BaX.h"
#include "../../tools/springio/springutils.cc"

using namespace std;

// Initial configuration of points and bonds; random six-fold coordinated network with springs of spring constant 1
string ptfile = "smd.pts";
string bfile = "smd.bds";
ldf boxsize = 45.5;  // box size for which the initial configuration is relaxed under periodic BC
ldf springconst = 1.; // spring constant (column 4 in smd.bds)
ldf Y = 2.*springconst/sqrt(3.);  // 2D Young's modulus
ldf nu = 1./3.;                   // Poisson ratio



int main(int argc, char* argv[])
{   
    bool pngout=true;
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    bmp.fillup(BLACK);

    
    // Parse command-line arguments
    if (argc < 2) { cout <<" Usage: dissipativenetwork <shear rate>" << endl; return EXIT_FAILURE; }
    ldf gammadot = atof(argv[1]); // shear rate, applied to boundaries perpendicular to x
    
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
    
    // add a dissipative spring force type, using the neighbor list from the input file
    vector<ldf> dissipativecoeff = {1.};
    ui dissipationForceIndex = sys.add_forcetype(EXTFORCE_DISSIPATION,&springnbrs,&dissipativecoeff);
    sys.assign_all_forcetype(dissipationForceIndex);
    
    // shear the box according to user-input shear rate
    sys.simbox.shear_boundary(1,0,gammadot*boxsize); 
    
    // initialize particle velocities to follow the shear rate
    ldf vx[systemsize];
    ldf vy[systemsize];
    for (ui i = 0; i < systemsize; i++) {
        vx[i]=0.;
        vy[i]=gammadot*x[i];
    }
    sys.import_vel(vx,vy);
    

    for(ui h=0;h<numeric_limits<ui>::max();h++)
    {
        sys.timesteps(1000);
        
        if (pngout) {
            for(ui i=0;i<systemsize;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
            bmp.save_png_seq(const_cast<char *>("sim"));
        }
        
        
        // theoretical values of sxx, sxy and syy for a 6-coordinated spring network
        ldf gamma = sys.simbox.Lshear[1][0]/boxsize; // instantaneous shear due to accumulated shear rate
        ldf sxy = gamma*Y/(2*(1+nu));     // sigma_xy
        ldf sxx = 0;           // also equal to sigma_yy
        
        cout << h <<" ";
        vector<double> sij = stress_tensor(sys);
        cout << "\tshear: " << gamma << " "
             << "\tdissipative sij: " 
             << sij[0]/(boxsize*boxsize)-sxx << " "<< sij[1]/(boxsize*boxsize)-sxy << " "
             << sij[2]/(boxsize*boxsize)-sxy << " "<< sij[3]/(boxsize*boxsize)-sxy 
             << " \tsim/theory: " 
             << sij[0]/(boxsize*boxsize)/sxx << " "<< sij[1]/(boxsize*boxsize)/sxy << " "<< sij[3]/(boxsize*boxsize)/sxx<<endl; 
        
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

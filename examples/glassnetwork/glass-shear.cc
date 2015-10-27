///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LIBMD__LONG_DOUBLE__
#define LIBMD__LONG_DOUBLE__
#endif

#include "../../libmd.h"
#include "../../tools/BaX/BaX.h"
#include "../../tools/springio/springio.cc"

using namespace std;

string ptfile = "glass.pts";
string bfile = "glass.bds";

int main()
{
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    bmp.fillup(BLACK);

    // initialize points
    // need correct-sized arrays
    ui systemsize = number_of_lines(ptfile);
    cout << systemsize << endl;
    ldf x[systemsize];
    ldf y[systemsize];

    ldf* dx = new ldf[systemsize];
    ldf* dy = new ldf[systemsize];

    for (ui i = 0; i < systemsize; i++) { dx[i] = 0.; dy[i] = 0.; }

    // make md system
    md<2> sys(systemsize);
    //~ sys.set_damping(1.0); //This should damp....
    sys.set_rco(15.);
    sys.set_ssz(15.);

    ldf lx = 76.68;
    ldf ly = 59.0282;
    sys.simbox.L[0]=lx;
    sys.simbox.L[1]=ly;
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;

    // initialize points
    read_points(ptfile,x,y);
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);


    // initialize bonds
    read_bonds(bfile,sys);

    // choose indexing algorithm
    sys.index();
    sys.network.update=false;

    // impose shear
    sys.simbox.shear_boundary(1,0,0.1);

    for(ui h=0;h<200;h++)
    {
        write_data("sim"+std::to_string(h), sys);

        sys.timesteps(1000);
        cout << h << endl;
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

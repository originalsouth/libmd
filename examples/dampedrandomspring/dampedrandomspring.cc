///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"
#include "../../tools/springio/springio.cc"

using namespace std;

string ptfile = "smd.pts";
string bfile = "smd.bds";

bool pngout=false;

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
    
    read_points(ptfile,x,y);
    
    // make md system
    md<2> sys(systemsize);
    sys.set_rco(15.);
    sys.set_ssz(15.);
    
    ldf boxsize = 50.5;
    sys.simbox.L[0]=boxsize; // make a larger box to force relaxation. sytem size for no forces: 45.5
    sys.simbox.L[1]=boxsize;
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.integrator.method=INTEGRATOR::VVERLET;
    sys.import_pos(x,y);

    
    // initialize bonds
    vector<vector<ui>> springnbrs(sys.N);
    read_bonds(bfile,sys,springnbrs);

    // index once for a static spring network
    sys.indexdata.method = INDEX::CELL;
    sys.index();
    sys.network.update=false;
    
    // add a dissipative spring force type, using the neighbor list from the input file
    vector<ldf> dissipativecoeff = {0.01};
    ui dissipationForceIndex = sys.add_forcetype(EXTFORCE::DISSIPATION,&springnbrs,&dissipativecoeff);
    sys.assign_all_forcetype(dissipationForceIndex);
    
    // timestep
    sys.integrator.h = 0.0001;

    for(ui h=0;h<10;h++)
    {
        if (pngout) {
            for(ui i=0;i<systemsize;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
            bmp.save_png_seq(const_cast<char *>("sim"));
        }
        
        //~ write_points_x("sim"+std::to_string(h)+".pts", sys);
        //~ write_bonds("sim"+std::to_string(h)+".bds", sys);
        write_data("sim"+std::to_string(h), sys);
        
        sys.timesteps(100000);
        cout << h << endl;
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

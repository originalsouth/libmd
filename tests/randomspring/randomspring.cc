///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../BaX/BaX.h"
#include "springio.cc"

using namespace std;


//~ string ptfile = "triangle.pts";
//~ string bfile = "triangle.bds";
string ptfile = "smd.pts";
string bfile = "smd.bds";

int main()
{   
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN,BLUE};
    bmp.fillup(BLACK);
    
    // initialize points
    // need correct-sized arrays
    ui systemsize = number_of_lines(ptfile);
    cout << systemsize << endl;
    ldf x[systemsize];
    ldf y[systemsize];
    
    read_points_ulrich(ptfile,x,y);
    
    // make md system
    md<2> sys(systemsize);
    sys.parallel.set(1);
    sys.network.rcosq=400.0;
    sys.network.rco=20.0;
    sys.network.sszsq=500.0;
    sys.simbox.L[0]=50.5; // make a larger box to force relaxation. sytem size for no forces: 45.5
    sys.simbox.L[1]=50.5;
    sys.simbox.bcond[0]=1;
    sys.simbox.bcond[1]=1;
    sys.integrator.method=1;
    sys.import_pos(&x,&y);


    // worst case: each particle a unique type, each bond a unique interaction
    for (int i = 0; i < systemsize; i++) sys.particles[i].type = i;
    
    // initialize bonds
    read_bonds_ulrich(bfile,sys);

    // seg faults with thomas's cell method (seg faults if the following line is commented). why?
    sys.indexdata.method = 1;
    sys.index();
    sys.network.update=false;

    for(ui h=0;h<100;h++)
    {
        for(ui i=0;i<systemsize;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
        bmp.save_png_seq(const_cast<char *>("sim"));
        sys.timesteps(100);
        cout << h << endl;
    }
    for(ui i=0;i<systemsize;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,BLUE);
    bmp.save_png_seq(const_cast<char *>("sim"));
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

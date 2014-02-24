///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

ldf x[2]={-0.5,0.5};
ldf y[2]={0.0,0.0};
ldf dx[2]={0.0,0.0};
ldf dy[2]={-0.5,0.5};

int main()
{
    __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    bmp.fillup(BLACK);
    md<2> sys(2);
    sys.parallel.set(2);
    sys.set_rco(1.0);
    sys.set_ssz(2.0);
    sys.simbox.L[0]=10.0;
    sys.simbox.L[1]=10.0;
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.integrator.method=INTEGRATOR::VVERLET;
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    vector<ldf> a={-1.0};
    sys.add_typeinteraction(0,0,POT::POT_COULOMB,&a);
    sys.index();
    sys.network.update=false;
    for(ui h=0;h<500;h++)
    {
        for(ui i=0;i<2;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,pix[i]);
        bmp.save_png_seq(const_cast<char *>("sim"));
        sys.timesteps(10);
    }
    for(ui i=0;i<2;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,pix[i]);
    bmp.save_png_seq(const_cast<char *>("sim"));
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

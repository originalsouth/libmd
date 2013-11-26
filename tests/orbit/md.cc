///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../BaX/BaX.h"

using namespace std;

ldf x[2]={-0.5,0.5};
ldf y[2]={0.0,-0.0};
ldf dx[2]={0.0,0.0};
ldf dy[2]={-0.5,0.5};

int main()
{
    unsigned int W=2000,H=2000;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    bmp.fillup(BLACK);
    md<2> sys(2);
    sys.network.rcosq=25.0;
    sys.network.rco=5.0;
    sys.network.sszsq=30.0;
    sys.simbox.L[0]=10.0;
    sys.simbox.L[1]=10.0;
    sys.simbox.bcond[0]=1;
    sys.simbox.bcond[1]=1;
    sys.integrator.method=1;
    for(ui i=0;i<2;i++)
    {
        sys.particles[i].x[0]=x[i];
        sys.particles[i].x[1]=y[i];
        sys.particles[i].dx[0]=dx[i];
        sys.particles[i].dx[1]=dy[i];  
    }
    vector<ldf> a={-1.0};
    sys.add_typeinteraction(0,0,0,&a);
    sys.index();
    sys.network.update=false;
    for(ui h=0;h<1000;h++)
    {
        for(ui i=0;i<2;i++) bmp.set(2,W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,pix[i]);
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

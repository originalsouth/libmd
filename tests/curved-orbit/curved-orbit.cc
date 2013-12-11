///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../BaX/BaX.h"

using namespace std;

ldf x[2]={-0.5,0.5};
ldf y[2]={0.0,0.0};
ldf dx[2]={0.0,0.0};
ldf dy[2]={-0.5,0.5};

int main()
{
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    bmp.fillup(BLACK);
    mpmd<2> sys(2);
    sys.patch.setmp(1);
    sys.parallel.set(2);
    sys.network.rcosq=100.0;
    sys.network.rco=10.0;
    sys.network.sszsq=120.0;
    sys.simbox.L[0]=10.0;
    sys.simbox.L[1]=10.0;
    sys.simbox.bcond[0]=1;
    sys.simbox.bcond[1]=1;
    sys.integrator.method=0;
    sys.import_pos(&x,&y);
    sys.import_vel(&dx,&dy);
    sys.particles[0].xp[0]=sys.particles[0].x[0]-sys.particles[0].dx[0]*sys.integrator.h;
    sys.particles[0].xp[1]=sys.particles[0].x[1]-sys.particles[0].dx[1]*sys.integrator.h;
    sys.particles[1].xp[0]=sys.particles[1].x[0]-sys.particles[1].dx[0]*sys.integrator.h;
    sys.particles[1].xp[1]=sys.particles[1].x[1]-sys.particles[1].dx[1]*sys.integrator.h;
    vector<ldf> a={-1.0};
    sys.add_typeinteraction(0,0,0,&a);
    sys.index();
    sys.network.update=false;
    for(ui h=0;h<2000;h++)
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

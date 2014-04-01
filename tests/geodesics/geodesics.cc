///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

ui N=60;

int main()
{
    __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={CYAN,MAGENTA,BLUE,RED,YELLOW,GREEN};
    bmp.fillup(BLACK);
    mpmd<2> sys(N);
    sys.set_rco(1.0);
    sys.set_ssz(1.0);
    sys.simbox.L[0]=4.0;
    sys.simbox.L[1]=4.0;
    sys.simbox.bcond[0]=BCOND::NONE;
    sys.simbox.bcond[1]=BCOND::NONE;
    sys.patch.setmp(MP::GAUSSIANBUMP);
    sys.integrator.method=MP_INTEGRATOR::VZ;
    ldf x[N],y[N],dx[N],dy[N];
    for(ui i=0;i<N;i++) x[i]=-sys.simbox.L[0]/2.0,y[i]=-sys.simbox.L[1]/2.0+i*sys.simbox.L[1]/N;
    for(ui i=0;i<N;i++) dx[i]=0.04,dy[i]=0.0;
    sys.integrator.h=1e-2;
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    sys.history();
    sys.network.update=false;
    FILE *energy;
    energy=fopen("energy.ls","w");
    for(ui h=0;h<1000;h++)
    {
        for(ui i=0;i<N;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H-(H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2),pix[i%6]);
        fprintf(energy,"%u;%Lf;%Lf;%Lf\n",h,sys.V(),sys.T(),sys.H());
        sys.timesteps(10);
    }
    for(ui i=0;i<N;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H-(H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2),pix[i%6]);
    bmp.save_png(const_cast<char *>("geodesics"));
    fprintf(energy,"%u;%Lf;%Lf;%Lf\n",1000,sys.V(),sys.T(),sys.H());
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define PASS_WARNING
#define PASS_ERROR

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

const ui N=156;
ldf x[N+100],y[N+100],zero[N+100]={0.0};

int main()
{
    __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    bmp.fillup(BLACK);
    md<2> sys(N);
    srand48(666);
    sys.simbox.bcond[0]=BCOND::HARD;
    sys.simbox.bcond[1]=BCOND::HARD;
    sys.set_rco(1.0);
    sys.set_ssz(2.0);
    for(ui i=0;i<N;i++) x[i]=-drand48()*sys.simbox.L[0]/2.0;
    for(ui i=0;i<N;i++) y[i]=drand48()*sys.simbox.L[1]-sys.simbox.L[1]/2.0;
    sys.import_pos(x,y);
    sys.import_vel(zero,zero);
    vector<ldf> a={1.0,100.0};
    vector<ldf> b={100.0,sys.simbox.L[1]/100.0};
    sys.add_typeinteraction(0,0,POT::POT_YUKAWA,&a);
    sys.add_typeinteraction(0,1,POT::POT_YUKAWA,&a);
    sys.add_typeinteraction(1,1,POT::POT_HOOKIAN,&b);
    ldf spx[]={0.0,-sys.simbox.L[1]/2.0};
    for(ui i=0;i<100;i++) sys.add_sp_interaction(0,i,i+1,2);
    for(ui i=0;i<100;i++) sys.sp_ingest(0,0,sys.add_particle(spx)),sys.particles[N+i].type=1,spx[1]+=sys.simbox.L[1]/100.0;
    sys.export_pos(x,y);
    for(ui i=0;i<N;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[0]);
    for(ui i=N;i<N+100;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[1]);
    bmp.save_png_seq(const_cast<char *>("sim"));
    for(ui k=0;k<1000;k++)
    {
        sys.timesteps(1000);
        sys.export_pos(x,y);
        bmp.fillup(BLACK);
        for(ui i=0;i<N;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[0]);
        for(ui i=N;i<N+100;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[1]);
        bmp.save_png_seq(const_cast<char *>("sim"));
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

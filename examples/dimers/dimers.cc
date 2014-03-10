///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define PASS_WARNING
#define PASS_ERROR
#define DEBUG_LEVEL 0

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

ldf ix[2]={0.08,-0.08};
ldf iy[2]={0.08,-0.08};
ldf zero[2]={0.0,0.0};

int main()
{
    __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    mpmd<2> sys(2);
    sys.set_rco(0.75);
    sys.set_ssz(1.0);
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.import_pos(ix);
    sys.import_pos(iy);
    sys.import_vel(zero);
    sys.import_vel(zero);
    sys.history();
    vector<ldf> a={1.0,100.0};
    vector<ldf> b={100.0,0.16};
    sys.add_typeinteraction(0,0,POT::POT_YUKAWA,&a);
    sys.add_typeinteraction(0,1,POT::POT_YUKAWA,&a);
    sys.add_typeinteraction(1,1,POT::POT_HOOKEAN,&b);
    sys.sp_ingest(0,0);
    sys.sp_ingest(0,1);
    fprintf(stderr,"Mario");
    for(ui i=0;i<1;i++) 
    {
        ldf delx[2]={0.0,2.0};
        sys.clone_particles(0,delx);
    }
    fprintf(stderr,"Luigi");
    ldf x[sys.N],y[sys.N];
    sys.export_pos(x,y);
    bmp.fillup(BLACK);
    for(ui i=0;i<sys.N;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[0]);
    bmp.save_png_seq(const_cast<char *>("sim"));
    for(ui k=0;k<100000;k++)
    {
        sys.timesteps(1000);
        sys.export_pos(x,y);
        bmp.fillup(BLACK);
        for(ui i=0;i<sys.N;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[0]);
        bmp.save_png_seq(const_cast<char *>("sim"));
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

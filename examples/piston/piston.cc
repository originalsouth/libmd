///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

const ui N=17466;
ldf x[N],y[N],zero[N]={0.0};

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
    for(ui i=0;i<N;i++) x[i]=-drand48()*sys.simbox.L[0]/2.0;
    for(ui i=0;i<N;i++) y[i]=drand48()*sys.simbox.L[1]-sys.simbox.L[1]/2.0;
    sys.import_pos(x,y);
    sys.import_vel(zero,zero);
    for(ui i=0;i<N;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,pix[0]);
    bmp.save_png_seq(const_cast<char *>("sim"));
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

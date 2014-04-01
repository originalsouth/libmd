///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

ui N=64,M=100;

int main()
{
    __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    md<2> sys(N+M);
    sys.set_rco(1.0);
    sys.set_ssz(3.0);
    sys.simbox.bcond[0]=BCOND::HARD;
    sys.simbox.bcond[1]=BCOND::HARD;
    ldf x[N+M],y[N+M],dx[M+N],dy[M+N];
    for(ui i=0;i<N;i++) x[i]=(-(1.0+i)/((ldf)(N)))*sys.simbox.L[0]/2.0;
    for(ui i=0;i<N;i++) y[i]=2.0*x[i]+sys.simbox.L[1]/2.0;
    for(ui i=0;i<N;i++) dx[i]=0.1;
    for(ui i=0;i<N;i++) dy[i]=0.0;
    for(ui i=N;i<N+M;i++) x[i]=0.1;
    for(ui i=N;i<N+M;i++) y[i]=-sys.simbox.L[1]/2.0+(((ldf)i-N)*sys.simbox.L[1]/((ldf)M));
    for(ui i=N;i<N+M;i++) dx[i]=0.0;
    for(ui i=N;i<N+M;i++) dy[i]=0.0;
    for(ui i=N;i<N+M;i++) sys.set_type(i,1);
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    vector<ldf> a={1.0,100.0};
    vector<ldf> b={10.0,sys.simbox.L[1]/((ldf)M)};
    sys.add_typeinteraction(0,0,POT::POT_YUKAWA,&a);
    sys.add_typeinteraction(0,1,POT::POT_YUKAWA,&a);
    sys.add_typeinteraction(1,1,POT::POT_HOOKEAN,&b);
    sys.add_sptype();
    for(ui i=N;i<N+M;i++) sys.sp_ingest(0,0,i);
    for(ui i=N;i<N+M;i++) sys.add_sp_interaction(0,i-N,i-N+1,2);
    sys.export_pos(x,y);
    bmp.fillup(BLACK);
    for(ui i=0;i<N;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[0]);
    for(ui i=N;i<N+M;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[1]);
    bmp.save_png_seq(const_cast<char *>("sim"));
    for(ui k=0;k<100;k++)
    {
        sys.timesteps(1000);
        sys.export_pos(x,y);
        bmp.fillup(BLACK);
        for(ui i=0;i<N;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[0]);
        for(ui i=N;i<N+M;i++) bmp.solidkykel(2.0,W*x[i]/sys.simbox.L[0]+W/2.0,H*y[i]/sys.simbox.L[1]+H/2,pix[1]);
        bmp.save_png_seq(const_cast<char *>("sim"));
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define PASS_ERRORS
#define PASS_WARNINGS
#define DEBUG_LEVEL 3

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

ldf urand()
{

    static unsigned long int rseed=666UL;
    return (rseed=16807UL*rseed%2147483647UL)/2147483648.0;
}

ldf zero[2]={};

int main()
{
    __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    bmp.fillup(BLACK);
    mpmd<2> sys(0);
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.set_index_method(INDEX::BRUTE_FORCE);
    sys.set_rco(1.0);
    sys.set_ssz(2.0);
    vector<ldf> a={1.0,100.0};
    vector<ldf> b={10.0,sys.simbox.L[1]/16.0};
    sys.add_interaction(POT::POT_YUKAWA,&a);
    sys.add_interaction(POT::POT_HOOKEAN,&b);
    sys.add_typeinteraction(0,0,0);
    sys.add_typeinteraction(0,1,1);
    sys.add_typeinteraction(1,1,0);
    ldf spx[2]={};
    for(ui i=0;i<2;i++) sys.add_sp_interaction(0,i,i+1,1);
    for(ui i=0;i<2;i++) sys.sp_ingest(0,0,sys.add_particle(spx)),sys.set_type(i,1),spx[1]+=sys.simbox.L[1]/16.0;
    for(ui i=0;i<3;i++)
    {
        ldf spxx[2]={sys.simbox.L[0]*urand()/2.0,sys.simbox.L[0]*urand()/2.0};
        sys.clone_particles(i,spxx);
    }
    sys.import_vel(zero,zero);
    sys.history();
    for(ui i=0;i<sys.N;i++)
    {
        ldf x=sys.direct_readout(0,i,'x');
        ldf y=sys.direct_readout(1,i,'x');
        printf("%u;%Lf;%Lf\n",i,x,y);
        bmp.solidkykel(2.0,W*x/sys.simbox.L[0]+W/2.0,H-H*y/sys.simbox.L[1]+H/2,pix[1]);
    }
    bmp.save_png_seq(const_cast<char *>("sim"));
    for(ui k=0;k<1000;k++)
    {
        sys.index();
        for(ui i=0;i<sys.N;i++)
        {
            ldf x=sys.direct_readout(0,i,'x');
            ldf y=sys.direct_readout(1,i,'x');
            printf("%u;%Lf;%Lf\n",i,x,y);
        }
        sys.calc_forces();
        for(ui i=0;i<sys.N;i++)
        {
            ldf x=sys.direct_readout(0,i,'x');
            ldf y=sys.direct_readout(1,i,'x');
            printf("%u;%Lf;%Lf\n",i,x,y);
        }
        sys.integrate();
        for(ui i=0;i<sys.N;i++)
        {
            ldf x=sys.direct_readout(0,i,'x');
            ldf y=sys.direct_readout(1,i,'x');
            printf("%u;%Lf;%Lf\n",i,x,y);
        }
        bmp.fillup(BLACK);
        for(ui i=0;i<sys.N;i++)
        {
            ldf x=sys.direct_readout(0,i,'x');
            ldf y=sys.direct_readout(1,i,'x');
            bmp.solidkykel(2.0,W*x/sys.simbox.L[0]+W/2.0,H-H*y/sys.simbox.L[1]+H/2,pix[1]);
        }
        bmp.save_png_seq(const_cast<char *>("sim"));
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

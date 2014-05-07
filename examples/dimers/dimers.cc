///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

ldf urand()
{

    static unsigned long int rseed=666UL;
    return (rseed=16807UL*rseed%2147483647UL)/2147483648.0;
}

ldf zero[500]={};

int main()
{
        __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H),bg(W,H);
    color pix[]={RED,GREEN},cc;
    mpmd<2> sys(0);
    for(ui i=0;i<W;i++) for(ui j=0;j<H;j++)
    {
        ldf fx[]={(((ldf)i)/W-0.5)*sys.simbox.L[0],(((ldf)j)/H-0.5)*sys.simbox.L[1]};
        bg.upset(i,j,cc.scale(sys.patch.f(fx),BLACK,WHITE));
    }
    bmp.import(&bg);
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.set_index_method(INDEX::KD_TREE);
    sys.set_rco(1.0);
    sys.set_ssz(2.0);
    vector<ldf> a={1.0,100.0};
    vector<ldf> b={1000.0,sys.simbox.L[1]/160.0};
    sys.add_interaction(POT::YUKAWA,&a);
    sys.add_interaction(POT::HOOKEAN,&b);
    sys.add_typeinteraction(0,0,0);
    sys.add_typeinteraction(0,1,1);
    sys.add_typeinteraction(1,1,0);
    ldf spx[2]={};
    sys.add_superparticle(sys.add_sptype());
    sys.add_sp_interaction(0,0,1,1);
    for(ui i=0;i<2;i++) sys.sp_ingest(0,sys.add_particle(spx)),sys.set_type(i,1),spx[1]+=sys.simbox.L[1]/160.0;
    sys.import_vel(zero,zero);
    sys.history();
    for(ui i=0;i<249;i++)
    {
        ldf spxx[2]={sys.simbox.L[0]*urand()/2.0,sys.simbox.L[0]*urand()/2.0};
        sys.clone_sp(i,spxx);
    }
    sys.import_vel(zero,zero);
    sys.history();
    for(ui i=0;i<sys.N;i++)
    {
        ldf x=sys.direct_readout(0,i,'x');
        ldf y=sys.direct_readout(1,i,'x');
        bmp.solidkykel(2.0,W*x/sys.simbox.L[0]+W/2.0,H-H*y/sys.simbox.L[1]+H/2,pix[1]);
    }
    bmp.save_png_seq(const_cast<char *>("sim"));
    FILE *energy=fopen("energy.ls","w");
    fprintf(energy,"%u;%Lf;%Lf;%Lf\n",0,sys.V(),sys.T(),sys.H());
    for(ui k=0;k<1000;k++)
    {
        sys.timesteps(100);
        bmp.import(&bg);
        for(ui i=0;i<sys.N;i++)
        {
            ldf x=sys.direct_readout(0,i,'x');
            ldf y=sys.direct_readout(1,i,'x');
            bmp.solidkykel(2.0,W*x/sys.simbox.L[0]+W/2.0,H-H*y/sys.simbox.L[1]+H/2,pix[1]);
        }
        fprintf(energy,"%u;%Lf;%Lf;%Lf\n",k+1,sys.V(),sys.T(),sys.H());
        bmp.save_png_seq(const_cast<char *>("sim"));
    }
    fclose(energy);
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

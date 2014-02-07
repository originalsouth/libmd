///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

ldf x[4]={-1.5,-0.5,0.5,1.5};
ldf y[4]={0.0,0.0,0.0,0.0};
ldf dx[4]={0.0,0.0,0.0,0.0};
ldf dy[4]={0.0,-0.5,0.5,0.0};


template<ui dim> void print_interactions(md<dim> &sys) {
    for (ui i = 0; i < 4; i++) printf("point %d has type %d.\n", i, sys.particles[i].type);
    for (map<pair<ui,ui>,ui>::iterator it = sys.network.lookup.begin(); it != sys.network.lookup.end(); it++) {
        interactiontype itype = sys.network.library[it->second];
        printf("type1: %d type2: %d interaction: %d\n", (it->first).first, (it->first).second, itype.potential);
    }
    printf("\n");
}

int main()
{
    __libmd__info();
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={GREEN,RED,BLUE,GREEN};
    bmp.fillup(BLACK);
    md<2> sys(4);
    sys.set_rco(10.0);
    sys.set_ssz(15.0);
    sys.simbox.L[0]=10.0;
    sys.simbox.L[1]=10.0;
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    sys.integrator.method=INTEGRATOR::VVERLET;
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    vector<ldf> a={1.0};
    sys.add_typeinteraction(0,0,POT::POT_COULOMB,&a);
    print_interactions(sys);
    sys.add_spring(0,1,1,1);
    print_interactions(sys);
    sys.add_spring(1,2,1,1);
    print_interactions(sys);
    sys.add_spring(2,3,1,1);
    print_interactions(sys);
    sys.rem_bond(2,3);
    print_interactions(sys);
    sys.index();
    sys.network.update=false;
    
    
    for(ui h=0;h<1000;h++)
    {
        for(ui i=0;i<4;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,pix[i]);
        bmp.save_png_seq(const_cast<char *>("sim"));
        sys.timesteps(100);
    }
    for(ui i=0;i<2;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,pix[i]);
    bmp.save_png_seq(const_cast<char *>("sim"));
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

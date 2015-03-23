///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LIBMD__LONG_DOUBLE__
#define LIBMD__LONG_DOUBLE__
#endif

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

#define CHAINSIZE 100

using namespace std;

bool pngout = true;

template<ui dim> void print_interactions(md<dim> &sys) {
    for (ui i = 0; i < 4; i++) printf("point %d has type %d.\n", i, sys.particles[i].type);
    for (map<pair<ui,ui>,ui>::iterator it = sys.network.lookup.begin(); it != sys.network.lookup.end(); it++) {
        interactiontype itype = sys.network.library[it->second];
        printf("type1: %d type2: %d interaction: %d\n", (it->first).first, (it->first).second, itype.potential);
    }
    printf("\n");
}


template<ui dim> void printmatrix (ldf A[dim][dim])
{ for (ui i = 0; i < dim; i++)
    for (ui j = 0; j < dim; j++)
      fprintf(stdout,"% 9.9Lf" F_UC "", A[i][j], j<dim-1?' ':'\n');
  fprintf(stdout,"\n");
}

template<ui dim> void print_network(md<dim> &sys) {
    for (ui i=0; i < sys.N; i++) {
         for(ui j=sys.network.skins[i].size()-1;j<UI_MAX;j--) printf("%d's neighbor: %d\tdistsq: %1.4Lf\n",i,sys.network.skins[i][j].neighbor,sys.distsq(i,sys.network.skins[i][j].neighbor));
     }
}


ldf x[CHAINSIZE];
ldf y[CHAINSIZE]={0.0};
ldf dx[CHAINSIZE]={0.0};
ldf dy[CHAINSIZE]={0.0};

int main()
{
    unsigned int W=1000,H=50;
    md<2> sys(CHAINSIZE);
    sys.set_rco(1.1);
    sys.set_ssz(1.1);
    sys.simbox.L[0]=CHAINSIZE*1.0;
    sys.simbox.L[1]=5.0;

    // line up a chain of ptcs
    for (int i=0; i < CHAINSIZE; i++) x[i]=i*1.0-sys.simbox.L[0]/2.+0.5;
    // pulse to first ptc
    dx[0] = 0.5;

    sys.import_pos(x,y);
    sys.import_vel(dx,dy);

    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    vector<ldf> a={1.0,1.0,2.5}; // hertzian
    sys.add_typeinteraction(0,0,POT::ANHARMONICSPRING,a);


    sys.indexdata.method=INDEX::CELL;

    sys.index();

    print_network(sys);
    sys.network.update=false;

    sys.set_rco(20.);
    sys.set_ssz(20.);


    sys.integrator.method=INTEGRATOR::SEULER;



    //~ for (ui i = 0; i < 5; i++) {
        //~ sys.particles[i].xp[0]=sys.particles[i].x[0]-sys.particles[i].dx[0]*sys.integrator.h;
        //~ sys.particles[i].xp[1]=sys.particles[i].x[1]-sys.particles[i].dx[1]*sys.integrator.h;
    //~ }

    for(ui h=0;h<400;h++)
    {
        if (pngout) {
            bitmap bmp(W,H);
            bmp.fillup(BLACK);
            for(ui i=0;i<CHAINSIZE;i++) bmp.set(3,W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
            bmp.save_png(const_cast<char *>(("sim"+std::to_string(h)).c_str()));
        }

        sys.timesteps(100);
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

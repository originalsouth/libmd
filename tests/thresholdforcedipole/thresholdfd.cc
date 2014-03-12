///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

bool pngout = false;

template<ui dim> void print_interactions(md<dim> &sys) {
    for (ui i = 0; i < 4; i++) printf("point %d has type %d.\n", i, sys.particles[i].type);
    for (map<pair<ui,ui>,ui>::iterator it = sys.network.lookup.begin(); it != sys.network.lookup.end(); it++) {
        interactiontype itype = sys.network.library[it->second];
        printf("type1: %d type2: %d interaction: %d\n", (it->first).first, (it->first).second, itype.potential);
    }
    printf("\n");
}

ldf x[2]={-.5,.5};
ldf y[2]={2.0,2.0};
ldf vx = 0.1;
ldf dx[2]={vx,-vx};
ldf dy[2]={0.0,0.0};

template<ui dim> void printmatrix (ldf A[dim][dim])
{ for (ui i = 0; i < dim; i++)
    for (ui j = 0; j < dim; j++)
      fprintf(stdout,"% 9.9Lf%c", A[i][j], j<dim-1?' ':'\n');
  fprintf(stdout,"\n");
}

template<ui dim> void print_network(md<dim> &sys) {
    for (ui i=0; i < sys.N; i++) {
         for(ui j=sys.network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) printf("%d's neighbor: %d\tdistsq: %1.4Lf\n",i,sys.network.skins[i][j].neighbor,sys.distsq(i,sys.network.skins[i][j].neighbor));
     }
}

int main()
{
    unsigned int W=500,H=500;
    md<2> sys(2);
    sys.set_rco(2.5);
    sys.set_ssz(2.5);
    sys.simbox.L[0]=5.0;
    sys.simbox.L[1]=5.0;
    
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;

    vector<ldf> hfd={1.0,1.0,0.1,.05};
    sys.add_bond(0,1,POT::POT_HOOKEANFORCEDIPOLE,&hfd);
    
    sys.indexdata.method=INDEX::CELL;
    
    sys.index();
    
    print_network(sys);
    sys.network.update=false;
    
    
    
    sys.integrator.method=INTEGRATOR::SEULER;
    
    
    for(ui h=0;h<400;h++)
    {
        for (ui i = 0; i < 2; i++) fprintf(stdout,"pos %1.8Lf ",sys.particles[i].x[0]);
        fprintf(stdout,"\n");
        for (ui i = 0; i < 2; i++) fprintf(stdout,"vel %1.8Lf ",sys.particles[i].dx[0]);
        fprintf(stdout,"\n");
        for (ui i = 0; i < 2; i++) fprintf(stdout,"force %1.8Lf ",sys.particles[i].F[0]);
        fprintf(stdout,"\n");
        fprintf(stdout,"\n");
        
        if (pngout) {
            bitmap bmp(W,H);
            bmp.fillup(BLACK);
            for(ui i=0;i<2;i++) bmp.set(3,W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
            bmp.save_png(const_cast<char *>(("sim"+std::to_string(h)).c_str()));
        }

        sys.timesteps(100);
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
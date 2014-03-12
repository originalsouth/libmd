///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../BaX/BaX.h"
#include "springio.cc"

using namespace std;

template<ui dim> void printmatrix (ldf A[dim][dim])
{ for (ui i = 0; i < dim; i++)
    for (ui j = 0; j < dim; j++)
      fprintf(stdout,"% 9.9Lf%c", A[i][j], j<dim-1?' ':'\n');
  fprintf(stdout,"\n");
}

//~ string ptfile = "triangle.pts";
//~ string bfile = "triangle.bds";
string ptfile = "smd.pts";
string bfile = "smd.bds";

int main()
{   
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN,BLUE};
    bmp.fillup(BLACK);
    
    // initialize points
    // need correct-sized arrays
    ui systemsize = number_of_lines(ptfile);
    cout << systemsize << endl;
    ldf x[systemsize];
    ldf y[systemsize];
    
    read_points_ulrich(ptfile,x,y);
    
    // make md system
    md<2> sys(systemsize);
    sys.parallel.set(1);
    sys.network.rcosq=400.0;
    sys.network.rco=20.0;
    sys.network.sszsq=500.0;
    
    ldf boxsize = 45.5;
    sys.simbox.L[0]=boxsize; // make a larger box to force relaxation. sytem size for no forces: 45.5
    sys.simbox.L[1]=boxsize;
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    
    // try shearing
    sys.simbox.shear_boundary(1,0,-0.1);
    
    
    
    
    sys.integrator.method=1;
    sys.import_pos(x,y);

    
    // initialize bonds
    read_bonds_ulrich(bfile,sys);

    // choose indexing algorithm
    sys.indexdata.method = 1;
    sys.index();
    
    // test deletion of a particle
    sys.rem_particle(5);
    sys.rem_particle(10);
    sys.rem_particle(15);
    sys.rem_particle(20);
    sys.network.update=false;

    for(ui h=0;h<10;h++)
    {
        for(ui i=0;i<sys.N;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
        bmp.save_png_seq(const_cast<char *>("sim"));
        
        write_points("sim"+std::to_string(h)+".pts", sys);
        write_bonds("sim"+std::to_string(h)+".bds", sys);
        
        sys.timesteps(1000);
        cout << h << endl;
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
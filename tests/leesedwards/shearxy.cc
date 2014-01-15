///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../BaX/BaX.h"

using namespace std;

ldf x[5]={-2., -1., 0., 1., 2.};
ldf y[5]={0.0,0.0,0.0,0.0,0.0};
ldf dx[5]={0.0,0.0,0.0,0.0,0.0};
//~ ldf dx[5]={0.05,-.05,0.0,0.0,0.0};
//~ ldf dy[5]={0.01,-.01,0.0,0.0,0.0};
ldf dy[5]={0.0,0.0,0.0,0.0,0.0};

template<ui dim> void printmatrix (ldf A[dim][dim])
{ for (ui i = 0; i < dim; i++)
    for (ui j = 0; j < dim; j++)
      fprintf(stderr,"% 9.9Lf%c", A[i][j], j<dim-1?' ':'\n');
  fprintf(stderr,"\n");
}

int main()
{
    unsigned int W=500,H=500;
    bitmap bmp(W,H);
    color pix[]={RED,GREEN};
    bmp.fillup(BLACK);
    md<2> sys(5);
    sys.parallel.set(2);
    sys.network.rcosq=1.21;
    sys.network.rco=1.1;
    sys.network.sszsq=1.21;
    sys.simbox.L[0]=5.0;
    sys.simbox.L[1]=5.0;
    
    // testing lees-edwards: shear boundaries perpendicular to x in y-direction
    // first index using ordinary PBC
    sys.simbox.bcond[0]=BCOND::PERIODIC;
    sys.simbox.bcond[1]=BCOND::PERIODIC;
    vector<ldf> a={1.0,1.0};
    sys.add_typeinteraction(0,0,2,&a);
    sys.index();
    sys.network.update=false;
    
    // now set shear Lees-Edwards of x boundary in y direction
    sys.simbox.shear_boundary(1,0,-0.01);
    
    sys.integrator.method=1;
    
        
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    for (ui i = 0; i < 5; i++) {
        sys.particles[i].xp[0]=sys.particles[i].x[0]-sys.particles[i].dx[0]*sys.integrator.h;
        sys.particles[i].xp[1]=sys.particles[i].x[1]-sys.particles[i].dx[1]*sys.integrator.h;
    }
    
    for(ui h=0;h<200;h++)
    {
        for (ui i = 0; i < 5; i++) fprintf(stderr,"%1.8Lf ",sys.particles[i].x[0]);
        fprintf(stderr,"\n");
        for (ui i = 0; i < 5; i++) fprintf(stderr,"%1.8Lf ",sys.particles[i].x[1]);
        fprintf(stderr,"\n");
        
        fprintf(stderr,"\n");
        //~ fprintf(stderr,"\n");
        for(ui i=0;i<5;i++) bmp.set(3,W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
        bmp.save_png_seq(const_cast<char *>("sim"));
        sys.timesteps(5000);
        //~ cout << sys.simbox.xshear[0] << endl;
    }
    for(ui i=0;i<5;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
    bmp.save_png_seq(const_cast<char *>("sim"));
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LIBMD__LONG_DOUBLE__
#define LIBMD__LONG_DOUBLE__
#endif

#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"

using namespace std;

bool pngout = false; // change to 'true' if also want PNG output

ldf x[2]={0.,1.};
ldf y[2]={0.0,0.0};


int main(int argc, char* argv[])
{
    // Read the particle velocity from command-line
    if (argc < 3) { cout << "Syntax: dissipativespring <dampingratio> <v0>" << endl; exit(0); }
    ldf xi = atof(argv[1]);
    ldf vx = atof(argv[2]);
    ldf dx[2]={0.,vx};
    ldf dy[2]={0.,0.};

    unsigned int W=500,H=500;
    md<2> sys(2);
    sys.set_rco(1.1);
    sys.set_ssz(1.1);
    sys.simbox.L[0]=5.0;
    sys.simbox.L[1]=5.0;

    // pin particle zero as base of spring
    sys.fix_particle(0,true);

    sys.import_pos(x,y);
    sys.import_vel(dx,dy);

    vector<ldf> a={1.0,1.0};
    sys.add_spring(0,1,1.,1.);

    // add dissipative component
    vector<vector<ui>> springnbrs(sys.N);
    springnbrs[0].push_back(1);
    springnbrs[1].push_back(0);

    // add a dissipative spring force type, using the neighbor list from the input file
    vector<ldf> dissipativecoeff = {2*xi};
    ui dissipationForceIndex = sys.add_forcetype(EXTFORCE::DISSIPATION,&springnbrs,&dissipativecoeff);
    sys.assign_all_forcetype(dissipationForceIndex);

    sys.indexdata.method=INDEX::CELL;

    sys.index();

    sys.network.update=false;

    sys.set_rco(2.5);
    sys.set_ssz(2.5);


    sys.integrator.method=INTEGRATOR::VVERLET;

    ldf t = 0.;
    ui elapsed_steps = 0;

    for(ui h=0;h<400;h++)
    {
        t = elapsed_steps*sys.integrator.h;
        for (ui i = 0; i < 2; i++) fprintf(stdout,"t %1.8Lf pos %1.8Lf %1.8Lf ",t,sys.particles[i].x[0],sys.particles[i].x[1]);
        fprintf(stdout,"\n");
        for (ui i = 0; i < 2; i++) fprintf(stdout,"t %1.8Lf vel %1.8Lf %1.8Lf ",t,sys.particles[i].dx[0],sys.particles[i].dx[1]);
        fprintf(stdout,"\n");
        for (ui i = 0; i < 2; i++) fprintf(stdout,"t %1.8Lf force %1.8Lf %1.8Lf ",t,sys.particles[i].F[0],sys.particles[i].F[1]);
        fprintf(stdout,"\n");
        fprintf(stdout,"\n");

        if (pngout) {
            bitmap bmp(W,H);
            bmp.fillup(BLACK);
            for(ui i=0;i<2;i++) bmp.set(3,W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
            bmp.save_png(const_cast<char *>(("sim"+std::to_string(h)).c_str()));
        }

        sys.timesteps(100);
        elapsed_steps += 100;
    }
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**************************************************************************************
 * Test of a rectangular skewed box with hard boundary conditions
 * 
 * Use the boxshear structure to simulate a simulation box in the form of a rhombus
 * with side length 5.0 and internal angle \pi/3. A single particle is placed in
 * the centre of the box and given a velocity parallel to one side, which should
 * lead to a periodic orbit as it is reflected around the box. At the end, the particle
 * position is given after a time has elapsed corresponding to 
 * 
 * Set pngout to 'true' to get snapshots of the simulation, and dataout to 'true' to get
 * particle position and velocity at intermediate steps. You may want to reduce the number
 * of cycles with pngout=true, or the image files will take up a lot of disk space.
 * ************************************************************************************/

#include <iostream>
#include "../../libmd.cc"
#include "../../tools/BaX/BaX.h"


using namespace std;

bool pngout = false;
bool dataout = false;


ldf L = 5.0; // Box size (side of rhombus)
ldf vx = 0.1; // Initial particle velocity (only in x direction)

int main()
{
    unsigned int W=500,H=round(500*sqrt(3.)/2.);
    md<2> sys(1);
    sys.simbox.L[0]=L;
    sys.simbox.L[1]=L*sqrt(3.)/2.;
    sys.set_rco(2.5);
    sys.set_ssz(2.5);
    
    ldf x[1] = {0.}; ldf y[1] = {0.};ldf dx[1] = {vx}; ldf dy[1] = {0.};
    sys.import_pos(x,y);
    sys.import_vel(dx,dy);
    
    sys.indexdata.method=INDEX::CELL;
    
    sys.index();
    sys.integrator.method=INTEGRATOR::VVERLET;
    
    sys.simbox.skew_boundary(0,1,2.5);     // shear box statically by length 2.5 along x direction
    
    sys.simbox.bcond[0]=BCOND::HARD;
    sys.simbox.bcond[1]=BCOND::HARD;
    
    ui n_cycles = 100;
    ui loops_per_cycle = 150; // number of main loops below that brings the particle through a complete loop, calculated from the integrator timestep and box size.
    
    for(ui h=0;h<n_cycles*loops_per_cycle;h++)
    {   
        sys.timesteps(1000);

        if (dataout)
        {    
            for (ui i = 0; i < sys.N; i++) fprintf(stdout,"pos %1.8Lf %1.8Lf ",sys.particles[i].x[0],sys.particles[i].x[1]);
            ldf sx=0,sy=0;
            for (ui i=0; i<2; i++)  { sx += sys.simbox.LshearInv[0][i]*sys.particles[0].x[i]; sy += sys.simbox.LshearInv[1][i]*sys.particles[0].x[i]; }
            for (ui i = 0; i < sys.N; i++) fprintf(stdout,"spos %1.8Lf %1.8Lf\n",sx,sy);
        }

        if (pngout) {
            bitmap bmp(W,H);
            bmp.fillup(BLACK);
            for(ui i=0;i<sys.N;i++) bmp.solidkykel(3,W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
            bmp.save_png(const_cast<char *>(("sim"+std::to_string(h)).c_str()));
        }
    }
    printf("Particle position: x %1.8Le\ty %1.8Le\n(should be zero within numerical precision)\n",sys.particles[0].x[0], sys.particles[0].x[1]);
    
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

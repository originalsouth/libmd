///////////////////////////////////////////////////////////////////////////////////
//  Example that builds a spring network from input files and command-line args  //
//  Code contributed by Bryan Gin-ge Chen                                        //
///////////////////////////////////////////////////////////////////////////////////


#ifndef LIBMD__LONG_DOUBLE__
#define LIBMD__LONG_DOUBLE__
#endif

#include "../../libmd.h"
#include "../../tools/BaX/BaX.h"
#include "../../tools/springio/springio.cc"

#define DEFAULT_BOXSIZE 100.

using namespace std;

int main(int argc,char *argv[])
{   

    char ptfile[200];
    char bfile[200];
    char vfile[200];
    char outFile[200]="sim_out"; //,outFile2[200],outFileloop[200];
    /* set global parameters from command line */
    int count=0,assigned;
    unsigned int W,H, writesteps,ministeps;
    ldf boxsizex=DEFAULT_BOXSIZE,boxsizey=DEFAULT_BOXSIZE,xskew=0;
    bool fail=false,posassigned=false,bdsassigned=false,velassigned=false,makebmps=false,energyonly=false;

    if(argc>1) {
        do {
            count++;
            if(argv[count][0]=='-'){  // check if we have an option string -x
                switch(argv[count][1]){  // switch on the option

                    case 'h' : // "Help" -h 
                        printf("-s Lengthx Lengthy Xskew -t (total writes) (steps between writes) -p (positions file) -b (bonds file) -v (vel file) -o (output prefix) [-B bitmapW bitmapH -E (toggle for energy-only mode)]\n");
                        return EXIT_SUCCESS;
                        break;

                    case 's' : // "Size" -s LENGTHX LENGTHY XSKEW
                        if(count+3 <argc){
                            count++;
                            assigned=sscanf(argv[count],F_LDF, &boxsizex);
                            if(assigned!=1) { 
                                printf("assigning X length failed in slot %d\n",count); 
                                fail=true;
                            } 
                            count++;
                            assigned=sscanf(argv[count],F_LDF, &boxsizey);
                            if(assigned!=1) {
                                printf("assigning Y length failed in slot %d\n",count);
                                fail=true;
                            } 
                            count++;
                            assigned=sscanf(argv[count],F_LDF, &xskew);
                            if(assigned!=1) {
                                printf("assigning X skew failed in slot %d\n",count);
                                fail=true;
                            } 
                        } else {
                            count=argc;
                            printf("too few arguments or bad input to -s (box sizex) (boxsize y) (box x skew) in slot %d\n",count);
                            fail=true;
                        }
                        break;

                    case 't' : // "timesteps" -t (total writes) (steps between writes)
                        if(count+2<argc){
                            count++;
                            assigned=sscanf(argv[count],"%d", &writesteps);
                            if(assigned!=1) {
                                printf("assigning total number of writes failed in slot %d\n",count);
                                fail=true;
                            } 
                            count++;
                            assigned=sscanf(argv[count],"%d", &ministeps);
                            if(assigned!=1) {
                                printf("assigning number of timesteps between writes failed in slot %d\n",count);
                                fail=true;
                            } 
                        } else {
                            count=argc;
                            printf("too few arguments or bad input to -t (total file writes) (steps betwen writes) in slot %d\n",count);
                            fail=true;
                        }
                        break;

                    case 'p' : // "positions file" -p posFilename
                        if(count+1<argc){
                            count++;
                            assigned=sscanf(argv[count],"%s", ptfile);
                            if(assigned!=1) {
                                printf("assigning position filename failed in slot %d\n",count);
                                fail=true;
                            } else {
                                // inputfiletrig=1;  // note! this loads immediately!
                                posassigned=true;
                                printf("loading positions from %s\n",ptfile);
                                // input_surface(ptfile,x);
                            }
                        } else {
                            count=argc;
                            printf("too few arguments or bad input to -p (positions file name) in slot %d\n",count); 
                            fail=true;
                        }
                        break;                  

                    case 'v' : // "velocities file" -v velFilename
                        if(count+1<argc){
                            count++;
                            assigned=sscanf(argv[count],"%s", vfile);
                            if(assigned!=1) {
                                printf("assigning velocity filename failed in slot %d\n",count);
                                fail=true;
                            } else {
                                // inputfiletrig=1;  // note! this loads immediately!
                                velassigned=true;
                                printf("loading velocities from %s\n",vfile);
                                // input_surface(ptfile,x);
                            }
                        } else {
                            count=argc;
                            printf("too few arguments or bad input to -v (velocities file name) in slot %d\n",count); 
                            fail=true;
                        }
                        break;                  

                    case 'b' : // "bonds file" -p bondsFilename
                        if(count+1<argc){
                            count++;
                            assigned=sscanf(argv[count],"%s", bfile);
                            if(assigned!=1) {
                                printf("assigning bond filename failed in slot %d\n",count);
                                fail=true;
                            } else {
                                // inputfiletrig=1;  // note! this loads immediately!
                                bdsassigned=true;
                                printf("loading bonds from %s\n",bfile);
                                // input_surface(ptfile,x);
                            }
                        } else {
                            count=argc;
                            printf("too few arguments or bad input to -b (bonds filename) in slot %d\n",count); 
                            fail=true;
                        }
                        break;

                    case 'o' : // "Output prefix" -o outputprefix
                        if(count+1<argc){
                            count++;
                            assigned=sscanf(argv[count],"%s", outFile);
                            if(assigned!=1) {
                                printf("assigning output filename failed in slot %d\n",count);
                                fail=true;
                            }
                        } else { 
                            count=argc;
                            printf("too few arguments or bad input to -o (output filename prefix) in slot %d\n",count);
                            fail=true;
                        } 
                        break;

                    case 'B' :  // "bitmap" -B Width Height
                        if(count+2<argc){
                            count++;
                            assigned=sscanf(argv[count],"%d", &W);
                            if(assigned!=1) {
                                printf("assigning bitmap width failed in slot %d\n",count);
                                fail=true;
                            }
                            count++;
                            assigned=sscanf(argv[count],"%d", &H);
                            if(assigned!=1) {
                                printf("assigning bitmap height failed in slot %d\n",count);
                                fail=true;
                            }
                            if(!fail) {
                                makebmps=true;
                                cout << W << " by " << H << " bitmaps will be drawn." << endl; 
                            }
                        } else { 
                            count=argc;
                            printf("too few arguments or bad input to -B (bitmap width) (bitmap height) in slot %d\n",count);
                            fail=true;
                        } 
                        break;

                    case 'E' :  // "energy only" -E
                        energyonly=true; 
                        printf("Energy only mode activated.\n");
                        break;

                    default :
                        count++;
                        printf("bad option in slot %d\n",count);
                        fail=true;
                        break;
                }

            } else {  // not an option string?  move forward a slot and try again
                count++;
                printf("did not get an option as expected in slot %d\n",count);
                fail=true;
            } 

        } while(count+1<argc);
    } 

    if(fail) { 
        printf("Couldn't parse command line parameters; run with -h to see command line parameters.\n");
        return EXIT_FAILURE;
    }

    if(posassigned && velassigned && bdsassigned) {
        printf("X length of box = " F_LDF ", Y length of box = " F_LDF ", X skew = " F_LDF"\n", boxsizex,boxsizey,xskew);
        printf("Writing H, T, V in %s.energies file, %d lines, with %d timesteps between each line.\n", outFile, writesteps, ministeps);
        if(!energyonly) { 
            printf("Will also write positions and bonds in %d files, with %d timesteps between writes:\n", writesteps,ministeps);    
            printf("Positions in %s###.pts and bonds in %s###.bds.\n", outFile,outFile);
        }

        // initialize points
        // need correct-sized arrays
        ui systemsize = number_of_lines(ptfile);
        cout << "number of particles: " << systemsize << endl;
        ldf x[systemsize];
        ldf y[systemsize];
        ldf vx[systemsize];
        ldf vy[systemsize];

        if(!read_points(ptfile,x,y)) cout << "points successfully read!" << endl;
        else fail=true;

        if(!read_points(vfile,vx,vy)) cout << "velocities successfully read!"<< endl;
        else fail=true;

        // make md system
        md<2> sys(systemsize);
        //~ sys.set_damping(1.0); //This should damp....
        sys.set_rco(15.);
        sys.set_ssz(15.);

        sys.simbox.L[0]=boxsizex; // make a larger box to force relaxation. sytem size for no forces: 45.5
        sys.simbox.L[1]=boxsizey;
        sys.simbox.skew_boundary(0,1,xskew);
        sys.simbox.bcond[0]=BCOND::BOXSHEAR;
        sys.simbox.bcond[1]=BCOND::BOXSHEAR;
        // try shearing
        //~ sys.simbox.shear_boundary(1,0,-0.1);

        sys.integrator.method=1;
        sys.import_pos(x,y);
        sys.import_vel(vx,vy);

        // initialize bonds
        if(!read_bonds(bfile,sys)) cout << "bonds successfully read!" << endl; 
        else fail=true;

        if(!fail) {

            // choose indexing algorithm
            sys.indexdata.method = 1;
            sys.index();
            sys.network.update=false;

            cout << endl;

            if(makebmps) {
                bitmap bmp(W,H);
                bmp.fillup(BLACK);

                for(ui h=0;h<writesteps;h++)
                {
                    for(ui i=0;i<systemsize;i++) bmp.set(W*sys.particles[i].x[0]/sys.simbox.L[0]+W/2.0,H*sys.particles[i].x[1]/sys.simbox.L[1]+H/2,GREEN);
                    bmp.save_png_seq(const_cast<char *>(outFile));

                    if(!energyonly) {
                        write_points_x(outFile+std::to_string(h)+".pts", sys);
                        write_bonds(outFile+std::to_string(h)+".bds", sys);
                    }

                    if(h>0) append_energies(outFile+std::string(".energies"), sys);
                        else write_energies(outFile+std::string(".energies"), sys);

                    sys.timesteps(ministeps);
                    //        cout << h << '\r' << flush;
                    cout << "\x1b[A" << h+1 <<  ((h>0) ? (" files written.") : (" file written.") )<< endl;

                }

            } else {

                for(ui h=0;h<writesteps;h++)
                {
                    if(!energyonly) {
                        write_points_x(outFile+std::to_string(h)+".pts", sys);
                        write_bonds(outFile+std::to_string(h)+".bds", sys);
                    }

                    if(h>0) append_energies(outFile+std::string(".energies"), sys);
                        else write_energies(outFile+std::string(".energies"), sys);

                    sys.timesteps(ministeps);
                    //        cout << h << '\r' << flush;
                    cout << "\x1b[A" << h+1 <<  ((h>0) ? (" files written.") : (" file written.") )<< endl;
                }

            }
            return EXIT_SUCCESS;

        } else {
            cout << "File reading failed." << endl;
            return EXIT_FAILURE;
        }

    } else { 
        printf("Not enough input files; run with -h to see command line parameters.\n");
        return EXIT_FAILURE;
    }
}

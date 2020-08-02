///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The official libmd regression tester                                                                          //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef LIBMD__LONG_DOUBLE__
#define LIBMD__LONG_DOUBLE__
#endif

#include "../../libmd.h"
#include <ctime>

using namespace std;

const long double eps=sqrt(numeric_limits<ldf>::epsilon());

#define test_success return printf("%s: %s: ",__FILE__,__FUNCTION__)
#define test_fail return !printf("%s: %s: ",__FILE__,__FUNCTION__)

/* Regression test function template (see folder "rtgroups")
 * bool test_group_component()
 * {
 *     //TODO: Write test utility here
 *     if (checks_out)
 *       test_success;
 *     else
 *       test_fail;
 * }
 *
 * Include the function here.
 *
 * Dont forget to add your test to the switch loops in testunit
 * To make a new group (eg group #1) add this switch
 * case 1: switch(j)
 * {
 *     case 0: p=test_group_component();
 *     break;
 *     default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
 * }
 * To make a new component (eg component #2) in an existing group add this to the group switch
 *     case 2: p=test_group_component();
 * If you have doubts mail an AUTHOR
 */

#include "rtgroups.h"

#include "rtgroups/integrator/verlet.cc"
#include "rtgroups/integrator/seuler.cc"

#include "rtgroups/boxshear/shearxy.cc"

#include "rtgroups/orbit/orbit.cc"
#include "rtgroups/orbit/orbit-bf.cc"

#include "rtgroups/curved-orbit/curved-orbit.cc"

#include "rtgroups/indexer/indexing.cc"

#include "rtgroups/network/modify_interactions.cc"
#include "rtgroups/network/modify_sp_interactions.cc"
#include "rtgroups/network/modify_bonds.cc"
#include "rtgroups/network/modify_sp_bonds.cc"
#include "rtgroups/network/remove_particles.cc"
#include "rtgroups/network/clone_remove_sp.cc"

#include "rtgroups/autodiff/autodiff.cc"
#include "rtgroups/autodiff/autodiff2.cc"
#include "rtgroups/autodiff/autodiff2b.cc"

ui groups=7;
ui group_size[]={2,1,2,1,2,6,3};

struct testunit
{
    bool run_all()
    {
        bool retval=true;
        for(ui i=0;i<groups;i++) retval=(retval and run(i));
        return retval;
    }
    bool run(ui i)
    {
        if(i<groups)
        {
            bool retval=true;
            for(ui j=0;j<group_size[i];j++) retval=(retval and run(i,j));
            return retval;
        }
        else
        {
            printf("test_group_not_found(%d): " IO_BOLDRED "failed" IO_RESET ".\n",i);
            return false;
        }
    }
    bool run(ui i,ui j)
    {
        bool retval=false;
        switch(i) //Group switch
        {
            case 0: switch(j) //Integrator Component switch
            {
                case 0: retval=test_integrator_verlet(); break;
                case 1: retval=test_integrator_seuler(); break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
            }
            break;
            case 1: switch(j) //BoxShear Component switch
            {
                case 0: retval=test_boxshear_shearxy(); break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
            }
            break;
            case 2: switch(j) //Orbit Component switch
            {
                case 0: retval=test_orbit_orbit(); break;
                case 1: retval=test_orbit_orbit_bf(); break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
            }
            break;
            case 3: switch(j) //Curved Orbit Component switch
            {
                case 0: retval=test_curved_orbit_orbit(); break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
            }
            break;
            case 4: switch(j) //Indexer Component switch
            {
                case 0: retval=test_indexer_noshear(); break;
                case 1: retval=test_indexer_shear(); break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
            }
            break;
            case 5: switch(j)
            {
                case 0: retval=test_modify_interactions(); break;
                case 1: retval=test_modify_sp_interactions(); break;
                case 2: retval=test_modify_bonds(); break;
                case 3: retval=test_modify_sp_bonds(); break;
                case 4: retval=test_remove_particles(); break;
                case 5: retval=test_clone_remove_sp(); break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
            }
            break;
            case 6: switch(j)
            {
                case 0: retval=test_autodiff(); break;
                case 1: retval=test_autodiff2_gaussian_bump(); break;
                case 2: retval=test_autodiff2_standard_expr(); break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
            }
            break;
            default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return retval;
        }
        if(retval) printf(IO_BOLDGREEN "pass" IO_RESET ".\n");
        else printf(IO_BOLDRED "failed" IO_RESET ".\n");
        return retval;
    }
} tests;

int main(int argc,char *argv[])
{
    __libmd__info();
    printf("\n");
    bool retval=false;
    if(argc<2) retval=tests.run_all();
    else if(argc==2) retval=tests.run(strtoul(argv[1],nullptr,10));
    else if(argc==3) retval=tests.run(strtoul(argv[1],nullptr,10),strtoul(argv[2],nullptr,10));
    else printf("too_many_input_arguments(%d): " IO_BOLDRED "failed" IO_RESET ".\n",argc-1);
    if(retval) return EXIT_SUCCESS;
    else return EXIT_FAILURE;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The official libmd regression tester                                                                         //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The official libmd regression tester                                                                          //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"

using namespace std;

const long double eps=sqrt(numeric_limits<ldf>::epsilon());

#define test_success {printf("%s: %s: ",__FILE__,__FUNCTION__); return true;}
#define test_fail {printf("%s: %s: ",__FILE__,__FUNCTION__); return false;}

#define test_success2 return printf("%s: %s: ",__FILE__,__FUNCTION__) >= 0;
#define test_fail2 return printf("%s: %s: ",__FILE__,__FUNCTION__) < 0;

/* Regression test function template (see folder "rtgroups")
 * bool test_group_component()
 * {
 *     //TODO: Write test utility here
 *     if (checks_out)
 *       test_success
 *     else
 *       test_fail
 *     //No semicolons!
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

#include "rtgroups/integrator/verlet.cc"
#include "rtgroups/integrator/seuler.cc"

#include "rtgroups/boxshear/shearxy.cc"

#include "rtgroups/orbit/orbit.cc"
#include "rtgroups/orbit/orbit-bf.cc"

#include "rtgroups/indexer/indexing.cc"

#include "rtgroups/network/modify_interactions.cc"
#include "rtgroups/network/modify_bonds.cc"
#include "rtgroups/network/remove_particles.cc"

ui groups=5;
ui group_size[]={2,1,2,2,3};

struct testunit
{
    void run_all()
    {
        for(ui i=0;i<groups;i++) run(i);
    }
    void run(ui i)
    {
        for(ui j=0;j<group_size[i];j++) run(i,j);
    }
    void run(ui i,ui j)
    {
        bool p=false;
        switch(i) //Group switch
        {
            case 0: switch(j) //Integrator Component switch
            {
                case 0: p=test_integrator_verlet();
                break;
                case 1: p=test_integrator_seuler();
                break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
            }
            break;
            case 1: switch(j) //BoxShear Component switch
            {
                case 0: p=test_boxshear_shearxy();
                break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
            }
            break;
            case 2: switch(j) //Orbit Component switch
            {
                case 0: p=test_orbit_orbit();
                break;
                case 1: p=test_orbit_orbit_bf();
                break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
            }
            break;
            case 3: switch(j) //Indexer Component switch
            {
                case 0: p=test_indexer_noshear();
                break;
                case 1: p=test_indexer_shear();
                break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
            }
            break;
            default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
            case 4: switch(j)
            {
                case 0: p=test_modify_interactions();
                break;
                case 1: p=test_modify_bonds();
                break;
                case 2: p=test_remove_particles();
                break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
            }
        }
        if(p) printf(IO_BOLDGREEN "pass" IO_RESET ".\n");
        else printf(IO_BOLDRED "failed" IO_RESET ".\n");
    }
} tests;

int main(int argc,char *argv[])
{
    __libmd__info();
    printf("\n");
    if(argc<2) tests.run_all();
    if(argc==2) tests.run(atoi(argv[1]));
    if(argc==3) tests.run(atoi(argv[1]),atoi(argv[2]));
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The official libmd regression tester                                                                         //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

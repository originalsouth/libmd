///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The official libmd regression tester                                                                          //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"

using namespace std;

const long double eps=sqrt(numeric_limits<ldf>::epsilon());

/* Regression template
 * bool test_group_component()
 * {
 *     printf("%s: ",__FUNCTION__)
 *     //TODO: Write test utility here
 * }
 *
 * Dont forget to add your test to the switch loops in testunit
 * To make a new group (eg group #1) add this switch
 * case 1: switch(j)
 * {
 *     case 0: p=test_group_component();
 *     break;
 *     default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
 * }
 * To make a new component (eg compenent #2) in an excisting group add this to the group switch
 *     case 2: p=test_group_component();
 * If you have doubts mail an AUTHOR
 */

bool test_integrator_verlet()
{
    printf("%s: ",__FUNCTION__);
    ldf x[]={0.0};
    ldf y[]={0.0};
    ldf dx[]={0.1};
    ldf dy[]={0.1};
    md<2> test(1);
    test.network.update=false;
    test.integrator.h=1.0;
    test.integrator.method=INTEGRATOR::VVERLET;
    test.import_pos(x,y);
    test.import_vel(dx,dy);
    test.timesteps(10);
    test.export_pos(x,y);
    test.export_vel(dx,dy);
    if(fabs(dx[0]-0.1)<=eps and fabs(dy[0]-0.1)<=eps and fabs(x[0]-1.0)<=eps and fabs(y[0]-1.0)<=eps) return true;
    else return false;
}

bool test_integrator_seuler()
{
    printf("%s: ",__FUNCTION__);
    ldf x[]={0.0};
    ldf y[]={0.0};
    ldf dx[]={0.1};
    ldf dy[]={0.1};
    md<2> test(1);
    test.network.update=false;
    test.integrator.h=1.0;
    test.integrator.method=INTEGRATOR::SEULER;
    test.import_pos(x,y);
    test.import_vel(dx,dy);
    test.timesteps(10);
    test.export_pos(x,y);
    test.export_vel(dx,dy);
    if(fabs(dx[0]-0.1)<=eps and fabs(dy[0]-0.1)<=eps and fabs(x[0]-1.0)<=eps and fabs(y[0]-1.0)<=eps) return true;
    else return false;
}

struct testunit
{
    void run_all()
    {
        ui Ni=1;
        for(ui i=0;i<Ni;i++) run(i);
    }
    void run(ui i)
    {
        ui Nj=2;
        for(ui j=0;j<Nj;j++) run(i,j);
    }
    void run(ui i,ui j)
    {
        bool p=false;
        switch(i) //Group switch
        {
            case 0: switch(j) //Component switch
            {
                case 0: p=test_integrator_verlet();
                break;
                case 1: p=test_integrator_seuler();
                break;
                default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
            }
            break;
            default: printf("test_not_found(%d,%d): " IO_BOLDRED "failed" IO_RESET ".\n",i,j); return;
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

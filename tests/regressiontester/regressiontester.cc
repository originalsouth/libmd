///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The official libmd regression tester                                                                          //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"

using namespace std;

#define IO_RESET   "\033[0m"
#define IO_BLACK   "\033[30m"
#define IO_RED     "\033[31m"
#define IO_GREEN   "\033[32m"
#define IO_YELLOW  "\033[33m"
#define IO_BLUE    "\033[34m"
#define IO_MAGENTA "\033[35m"
#define IO_CYAN    "\033[36m"
#define IO_WHITE   "\033[37m"
#define IO_BOLDBLACK   "\033[1m\033[30m"
#define IO_BOLDRED     "\033[1m\033[31m"
#define IO_BOLDGREEN   "\033[1m\033[32m"
#define IO_BOLDYELLOW  "\033[1m\033[33m"
#define IO_BOLDBLUE    "\033[1m\033[34m"
#define IO_BOLDMAGENTA "\033[1m\033[35m"
#define IO_BOLDCYAN    "\033[1m\033[36m"
#define IO_BOLDWHITE   "\033[1m\033[37m"

const long double eps=sqrt(numeric_limits<ldf>::epsilon());

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
        switch(i)
        {
            case 0: switch(j)
            {
                case 0: p=test_integrator_verlet();
                break;
                case 1: p=test_integrator_seuler();
                break;
                default: printf("test_not_found: " IO_BOLDRED "failed" IO_RESET ".\n"); return;
            }
            break;
            default: printf("test_not_found: " IO_BOLDRED "failed" IO_RESET ".\n"); return;
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

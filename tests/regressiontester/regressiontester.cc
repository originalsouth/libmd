///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The official libmd regression tester                                                                          //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.cc"
#include "../BaX/BaX.h"

using namespace std;

const long double eps=sqrt(numeric_limits<ldf>::epsilon());

bool test_1()
{
    ldf x[]={0.0};
    ldf y[]={0.0};
    ldf dx[]={0.1};
    ldf dy[]={0.1};
    md<2> test(1);
    test.network.update=false;
    test.integrator.h=1.0;
    test.import_pos(x,y);
    test.import_vel(dx,dy);
    test.timesteps(10);
    test.export_pos(x,y);
    //test.export_pos(y);
    test.export_vel(dx,dy);
    printf("%Lf\n",x[0]);
    printf("%Lf\n",y[0]);
    printf("%Lf\n",dx[0]);
    printf("%Lf\n",dy[0]);
    if(fabs(dx[0]-0.1)<=eps and fabs(dy[0]-0.1)<=eps and fabs(x[0]-1.0)<=eps and fabs(y[0]-1.0)<=eps) return true;
    else return false;
}

int main()
{
    __libmd__info();
    printf("%d\n",test_1());
    return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The official libmd regression tester                                                                         //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

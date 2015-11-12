///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../../libmd.h"

using namespace std;

template<class X> X YUKAWA_1(X r,std::vector<ldf> &parameters)
{
    const ldf b=parameters[0];
    const ldf k=parameters[1];
    return b/(r*exp(k*r));
}

template<class X> X YUKAWA_2(X r,std::vector<ldf> &parameters)
{
    const ldf b=parameters[0];
    const ldf k=parameters[1];
    return b*exp(-k*r)/r;
}

template<class X> X YUKAWA_3(X r,std::vector<ldf> &parameters)
{
    const ldf b=parameters[0];
    const ldf k=parameters[1];
    return b*exp(-k*r)*pow(r,-1);
}

template<class X> X YUKAWA_4(X r,std::vector<ldf> &parameters)
{
    const ldf b=parameters[0];
    const ldf k=parameters[1];
    return b*exp(-k*r-log(r));
}

template<class X> X YUKAWA_5(X r,std::vector<ldf> &parameters)
{
    const ldf b=parameters[0];
    const ldf k=parameters[1];
    return b*exp(-(k*r+log(r)));
}

template<class X> X YUKAWA_6(X r,std::vector<ldf> &parameters)
{
    const ldf b=parameters[0];
    const ldf k=parameters[1];
    return b*exp(-k*r)*1.0/r;
}

template<class X> X YUKAWA_7(X r,std::vector<ldf> &parameters)
{
    const ldf b=parameters[0];
    const ldf k=parameters[1];
    return b/r*exp(-k*r);
}

int main(int argc,char *argv[])
{
    ui counter=1;
    const ui I=ui((argc>1)?atof(argv[1]):1e5);
    vector<ldf> param={1.0,10.0};
    vector<potentialptr<dual>> YUKAWAS={YUKAWA_1,YUKAWA_2,YUKAWA_3,YUKAWA_4,YUKAWA_5,YUKAWA_6,YUKAWA_7};
    for(potentialptr<dual> X: YUKAWAS)
    {
        TicToc();
        for(ui i=0;i<I;i++) for(ldf x=1e-8;x<=1.0;x+=0.001)
        {
            const dual ox=X(dual(x,1.0),param);
            printf("%f;%f;%f\n",x,ox.x,ox.dx);
        }
        fprintf(stderr,"YUKAWA_%u:\t%f(s)\n",counter++,TicToc());
    }
    return EXIT_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple test file                                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*this test tests the autodiff structure*/

template<class X> X LINEAR(X x)
{
    const ldf a=20.0;
    const ldf b=10.0;
    return a*x+b;
}

template<class X> X dLINEAR(X x)
{
    (void) x;
    const ldf a=20.0;
    return a;
}

template<class X> X QUOTIENT(X x)
{
    return (x-1.0)/(1.0+x*x);
}

template<class X> X dQUOTIENT(X x)
{
    return (1.0+2.0*x-pow(x,2))/pow(1.0+pow(x,2),2);
}

template<class X> X POWER(X x)
{
    return pow(x,3);
}

template<class X> X dPOWER(X x)
{
    return 3.0*pow(x,2);
}

template<class X> X EXP(X x)
{
    return exp(2.0*pow(x,2)+1.0);
}

template<class X> X dEXP(X x)
{
    return 4.0*x*EXP(x);
}

template<class X> X Th(X x)
{
    return heaviside(x);
}

template<class X> X dTh(X x)
{
    return (x==0.0)?numeric_limits<ldf>::infinity():0.0;
}

template<class X> using tadptr=X (*)(X);

bool test_autodiff()
{
    rseed=rseedb=time(NULL);
    vector<tadptr<dual>> func;
    vector<tadptr<ldf>> dfunc;
    func.push_back(LINEAR<dual>);
    dfunc.push_back(dLINEAR<ldf>);
    func.push_back(QUOTIENT<dual>);
    dfunc.push_back(dQUOTIENT<ldf>);
    func.push_back(POWER<dual>);
    dfunc.push_back(dPOWER<ldf>);
    func.push_back(EXP<dual>);
    dfunc.push_back(dEXP<ldf>);
    func.push_back(Th<dual>);
    dfunc.push_back(dTh<ldf>);
    for(ui i=0;i<1000;i++)
    {
        ui k=randnrb()%5;
        ldf x=randnr();
        if(k==5 and i<100) x=0.0;
        dual y(x,1.0);
        ldf z=func[k](y).dx;
        #if DEBUG_LEVEL>0
        printf("autodiff[debug]: %d " F_LDF " " F_LDF " " F_LDF " \n",k,x,z,dfunc[k](x));
        #endif
        if(fabs(z-dfunc[k](x))>numeric_limits<ldf>::epsilon()) test_fail;
    }
    test_success;
}

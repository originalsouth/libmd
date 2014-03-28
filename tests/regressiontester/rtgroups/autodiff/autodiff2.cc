/*this test tests the autodiff2 structure*/

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
    for(ui i=0;i<1000;i++)
    {
        ui k=randnrb()%4;
        ldf x=randnr();
        dual y(x,1.0);
        ldf z=func[k](y).dx;
        #if DEBUG_LEVEL>0
        printf("autodiff[debug]: %d %Lf %Lf %Lf \n",k,x,z,dfunc[k](x));
        #endif
        if(fabs(z-dfunc[k](x))>numeric_limits<ldf>::epsilon()) test_fail;
    }
    test_success;
}

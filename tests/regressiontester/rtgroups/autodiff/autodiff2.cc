/*this test tests the autodiff2 structure*/

//Gaussian bump
template<ui dim> ldf GAUSSIANBUMP(ldf *x,vector<ldf> *param)
{
    const ldf A=param->at(0);
    const ldf K=param->at(1);
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return A*exp(-retval*K);
}

template<ui dim> ldf dGAUSSIANBUMP(ui i,ldf *x,vector<ldf> *param)
{
    const ldf A=param->at(0);
    const ldf K=param->at(1);
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return -2.0*A*K*x[i]*exp(-retval*K);
}

template<ui dim> ldf ddGAUSSIANBUMP(ui i,ui j,ldf *x,vector<ldf> *param)
{
    const ldf A=param->at(0);
    const ldf K=param->at(1);
    const ldf kdel=kdelta(i,j);
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return 2.0*A*K*(2.0*K*x[i]*x[j]-kdel)*exp(-retval*K);
}

template<class X> using tadptr=X (*)(X);

bool test_autodiff2_gaussian_bump()
{
    rseed=rseedb=time(NULL);
    for(ui i=0;i<1000;i++)
    {
        ldf x[2]={(randnr()-0.5)*(randnrb()%10),(randnr()-0.5)*(randnrb()%10)};
        duals<dim> y[2],z[2];
        for(ui d=0;d<2;d++) y[d]=duals<dim>(x[d],d);

        if(fabs(z-dfunc[k](x))>numeric_limits<ldf>::epsilon()) test_fail;
    }
    test_success;
}

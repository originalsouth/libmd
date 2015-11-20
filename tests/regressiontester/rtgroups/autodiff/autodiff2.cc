/*this test tests the autodiff2 structure*/

//Gaussian bump
template<ui dim,class X> X GAUSSIANBUMP(X *x,vector<ldf> &param)
{
    const ldf A=param[0];
    const ldf K=param[1];
    X retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return A*exp(-retval*K);
}

template<ui dim> ldf dGAUSSIANBUMP(ui i,ldf *x,vector<ldf> &param)
{
    const ldf A=param[0];
    const ldf K=param[1];
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return -2.0*A*K*x[i]*exp(-retval*K);
}

template<ui dim> ldf ddGAUSSIANBUMP(ui i,ui j,ldf *x,vector<ldf> &param)
{
    const ldf A=param[0];
    const ldf K=param[1];
    const ldf kdel=(i==j)?1.0:0.0;
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return 2.0*A*K*(2.0*K*x[i]*x[j]-kdel)*exp(-retval*K);
}

bool test_autodiff2_gaussian_bump()
{
    const ldf eps=numeric_limits<ldf>::epsilon();
    rseed=time(NULL);
    vector<ldf> param(2,1.0);
    for(ui i=0;i<1000000;i++)
    {
        ldf x[2]={(urand()-0.5)*(irand()%10),(urand()-0.5)*(irand()%10)};
        duals<2> y[2],z;
        for(ui d=0;d<2;d++) y[d]=duals<2>(x[d],d);
        z=GAUSSIANBUMP<2>(y,param);
        #if DEBUG_LEVEL>0
        printf("autodiff2[debug]: " F_LDF " " F_LDF " " F_LDF " " F_LDF " \n",x[0],x[1],z.x,GAUSSIANBUMP<2>(x,param));
        printf("autodiff2[debug]: " F_LDF " " F_LDF " 0 " F_LDF " " F_LDF " \n",x[0],x[1],z.dx[0],dGAUSSIANBUMP<2>(0,x,param));
        printf("autodiff2[debug]: " F_LDF " " F_LDF " 1 " F_LDF " " F_LDF " \n",x[0],x[1],z.dx[1],dGAUSSIANBUMP<2>(1,x,param));
        printf("autodiff2[debug]: " F_LDF " " F_LDF " 0 0 " F_LDF " " F_LDF " \n",x[0],x[1],z.dxdy[0][0],ddGAUSSIANBUMP<2>(0,0,x,param));
        printf("autodiff2[debug]: " F_LDF " " F_LDF " 0 1 " F_LDF " " F_LDF " \n",x[0],x[1],z.dxdy[0][1],ddGAUSSIANBUMP<2>(0,1,x,param));
        printf("autodiff2[debug]: " F_LDF " " F_LDF " 1 0 " F_LDF " " F_LDF " \n",x[0],x[1],z.dxdy[1][0],ddGAUSSIANBUMP<2>(1,0,x,param));
        printf("autodiff2[debug]: " F_LDF " " F_LDF " 1 1 " F_LDF " " F_LDF " \n",x[0],x[1],z.dxdy[1][1],ddGAUSSIANBUMP<2>(1,1,x,param));
        #endif
        if(fabs(z.x-GAUSSIANBUMP<2>(x,param))>eps) test_fail;
        else if(fabs(z.dx[0]-dGAUSSIANBUMP<2>(0,x,param))>eps or fabs(z.dx[1]-dGAUSSIANBUMP<2>(1,x,param))>eps) test_fail;
        else if(fabs(z.dxdy[0][0]-ddGAUSSIANBUMP<2>(0,0,x,param))>eps or fabs(z.dxdy[0][1]-ddGAUSSIANBUMP<2>(0,1,x,param))>eps) test_fail;
        else if(fabs(z.dxdy[1][0]-ddGAUSSIANBUMP<2>(1,0,x,param))>eps or fabs(z.dxdy[1][1]-ddGAUSSIANBUMP<2>(1,1,x,param))>eps) test_fail;
    }
    test_success;
}

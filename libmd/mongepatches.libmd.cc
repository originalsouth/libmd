ldf kdelta(ui i,ui j)
{
    return (i==j)?1.0:0.0;
}

//Flat space
template<ui dim> ldf FLATSPACE(ldf *x,vector<ldf> *param)
{
    (void) x;
    (void) param;
    return 0.0;
}

template<ui dim> ldf dFLATSPACE(ui i,ldf *x,vector<ldf> *param)
{
    (void) i;
    (void) x;
    (void) param;
    return 0.0;
}

template<ui dim> ldf ddFLATSPACE(ui i,ui j,ldf *x,vector<ldf> *param)
{
    (void) i;
    (void) j;
    (void) x;
    (void) param;
    return 0.0;
}

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
    return 2.0*A*K*(K*x[i]*x[j]-kdel)*exp(-retval*K);
}

//Egg carton
//ldf EGGCARTON(ldf *x,vector<ldf> *param)
//{

//}

//ldf dEGGCARTON(ui i,ldf *x,vector<ldf> *param)
//{

//}

//ldf ddEGGCARTON(ui i,ui j,ldf *x,vector<ldf> *param)
//{

//}

//Mollifier
//ldf MOLLIFIER(ldf *x,vector<ldf> *param)
//{

//}

//ldf dMOLLIFIER(ui i,ldf *x,vector<ldf> *param)
//{

//}

//ldf ddMOLLIFIER(ui i,ui j,ldf *x,vector<ldf> *param)
//{

//}

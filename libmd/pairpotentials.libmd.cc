pairpotentials::pairpotentials()
{
    add(COULOMB,dCOULOMBdr);
    add(YUKAWA,dYUKAWAdr);
    add(HOOKIAN,dHOOKIANdr);
    add(LJ,dLJdr);
    add(MORSE,dMORSEdr);
}

ui pairpotentials::add(potentialptr p,potentialptr dpdr)
{
    potentials.push_back(p);
    dpotentialsdr.push_back(dpdr);
    return potentials.size()-1;
}

ldf pairpotentials::operator()(ui type,ldf r,ldf rsq,vector<ldf>* parameters)
{
    return (potentials[type])(r,rsq,parameters); 
}

ldf pairpotentials::dr(ui type,ldf r,ldf rsq,vector<ldf>* parameters)
{
    return (dpotentialsdr[type])(r,rsq,parameters);
}

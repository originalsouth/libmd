#ifndef libmd_h
#include "../libmd.h"
#endif

pairpotentials::pairpotentials()
{
    add(COULOMB);
    add(YUKAWA);
    add(HOOKIAN);
    add(LJ);
    add(MORSE);
    add(FORCEDIPOLE);
    add(HOOKEANFORCEDIPOLE);
}

ui pairpotentials::add(potentialptr p)
{
    potentials.push_back(p);
    return potentials.size()-1;
}

ldf pairpotentials::operator()(ui type,ldf r,vector<ldf>* parameters)
{
    dual rdx=r;
    return (potentials[type])(rdx,parameters).x;
}

ldf pairpotentials::dr(ui type,ldf r,vector<ldf>* parameters)
{
    dual rdx=r;
    return (potentials[type])(rdx,parameters).dx;
}

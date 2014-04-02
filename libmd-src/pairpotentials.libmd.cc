#ifndef libmd_h
#include "../libmd.h"
#endif

pairpotentials::pairpotentials()
{
    potentials.reserve(16);
    add(COULOMB<dual>);
    add(YUKAWA<dual>);
    add(HOOKEAN<dual>);
    add(LJ<dual>);
    add(MORSE<dual>);
    add(FORCEDIPOLE<dual>);
    add(HOOKEANFORCEDIPOLE<dual>);
    add(ANHARMONICSPRING<dual>);
}

ui pairpotentials::add(potentialptr<dual> p)
{
    potentials.push_back(p);
    return potentials.size()-1;
}

ldf pairpotentials::operator()(ui type,ldf r,vector<ldf>* parameters)
{
    dual rdx={r,1.0};
    return (potentials[type])(rdx,parameters).x;
}

ldf pairpotentials::dr(ui type,ldf r,vector<ldf>* parameters)
{
    dual rdx={r,1.0};
    return (potentials[type])(rdx,parameters).dx;
}

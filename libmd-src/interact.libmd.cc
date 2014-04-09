#ifndef libmd_h
#include "../libmd.h"
#endif

interactiontype::interactiontype(ui ppot,vector<ldf> *param,ldf Vco)
{
    potential=ppot;
    parameters=*param;
    vco=Vco;
}

interactionneighbor::interactionneighbor(ui noneighbor,ui nointeraction)
{
    neighbor=noneighbor;
    interaction=nointeraction;
}

interact::interact()
{
    update=true;
    rco=ssz=1.0;
}

pair<ui,ui> interact::hash(ui type1,ui type2)
{
    return (type2<type1)?pair<ui,ui>(type2,type1):pair<ui,ui>(type1,type2);
}

bool interact::probe(ui type1,ui type2)
{
    return static_cast<bool>(lookup.count(hash(type1,type2)));
}

interactiontype::interactiontype()
{

}

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

pair<ui,ui> interact::hash(ui type1,ui type2)
{
    if(type2<type1) return pair<ui,ui>(type2,type1);
    else return pair<ui,ui>(type1,type2);
}

#ifndef libmd_h
#include "../libmd.h"
#endif

interactiontype::interactiontype(ui ppot,vector<ldf> &param,ldf Rco,ldf Vco)
{
    //!
    //! Constructor for interactiontype.
    //! Sets <tt>potential</tt> to <tt>ppot</tt>, <tt>parameters</tt> to <tt>&param</tt> and <tt>vco</tt> to <tt>Vco</tt>.
    //!
    potential=ppot;
    parameters=param;
    rco=Rco;
    vco=Vco;
}

interactionneighbor::interactionneighbor(ui noneighbor,ui nointeraction)
{
    //!
    //! Constructor for interactionneighbor.
    //! Sets <tt>neighbor</tt> to <tt>noneighbor</tt> and <tt>interaction</tt> to <tt>nointeraction</tt>.
    //!
    neighbor=noneighbor;
    interaction=nointeraction;
}

interact::interact()
{
    //!
    //! Constructor for interact.
    //! Sets <tt>update</tt> to <tt>true</tt> and both <tt>rco</tt> and <tt>ssz</tt> to 1.
    //!
    update=true;
    rco=ssz=1.0;
}

pair<ui,ui> interact::hash(ui type1,ui type2)
{
    //!
    //! Returns a pair containing <tt>type1</tt> and <tt>type2</tt> with the largest of the two as the first.
    //!
    return (type2>type1)?pair<ui,ui>(type2,type1):pair<ui,ui>(type1,type2);
}

bool interact::probe(ui type1,ui type2)
{
    //!
    //! Returns whether there is an interaction associated with <tt>type1</tt> and <tt>type2</tt>.
    //!
    return static_cast<bool>(lookup.count(hash(type1,type2)));
}

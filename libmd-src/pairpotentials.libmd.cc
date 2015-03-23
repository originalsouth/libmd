#ifndef libmd_h
#include "../libmd.h"
#endif

pairpotentials::pairpotentials()
{
    //!
    //! pairpotentials constuctor. <br>
    //! Reserves 16 slot in the potentials vector. <br>
    //! Adds the builtin pairpotentials. <br>
    //!
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
    //!
    //! This function adds a custom potential function <tt>p</tt> to the potentials vector. <br>
    //!
    potentials.push_back(p);
    return potentials.size()-1;
}

ldf pairpotentials::operator()(ui type,ldf r,vector<ldf> *parameters)
{
    //!
    //! This function evaluates a potential function in <tt>pairpotentials[type]</tt> at <tt>r</tt> with <tt>parameters</tt> <br>
    //! Make sure that if you use a custom potential function <tt>type</tt> and <tt>potentials</tt> align.
    //!
    dual rdx(r,1.0);
    return (potentials[type])(rdx,parameters).x;
}

ldf pairpotentials::dr(ui type,ldf r,vector<ldf> *parameters)
{
    //!
    //! This function evaluates the derivative of a potential function in <tt>pairpotentials[type]</tt> at <tt>r</tt> with <tt>parameters</tt> <br>
    //! Make sure that if you use a custom potential function <tt>type</tt> and <tt>potentials</tt> align.
    //!
    dual rdx(r,1.0);
    return (potentials[type])(rdx,parameters).dx;
}

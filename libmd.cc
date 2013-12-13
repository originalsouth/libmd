///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of LIBRARY source                                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __LIBMD__
#define __LIBMD__

#include "libmd.h"

#include "libmd/threads.libmd.cc"               //This file implements the thread structure
#include "libmd/autodiff.libmd.cc"              //This file implements automatic differentation
#include "libmd/potentials.libmd.cc"            //This file has all the builtin pairpotential functions and derivatives
#include "libmd/particle.libmd.cc"              //This file implements the particle structure
#include "libmd/box.libmd.cc"                   //This file implements the box structure
#include "libmd/interact.libmd.cc"              //This file implements the interact structure
#include "libmd/index.libmd.cc"                 //This file takes care of indexing algorithms
#include "libmd/pairpotentials.libmd.cc"        //This file implements the pairpotential structure
#include "libmd/integrators.libmd.cc"           //This file implements the integration structure
#include "libmd/md.libmd.cc"                    //This file implements the md structure which takes care of molecular dynamics in flat space
#include "libmd/mongepatches.libmd.cc"          //This file has all the builtin monge patch functions and derivatives
#include "libmd/mp.libmd.cc"                    //This file implements the mp structure
#include "libmd/mpmd.libmd.cc"                  //This file implements the mpmd structure which takes care of molecular dynamics on monge patches

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of LIBRARY source                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

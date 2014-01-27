///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of LIBRARY source                                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __LIBMD__
#define __LIBMD__
#if __cplusplus < 201103L
#warning "C++11 not detetected: libmd needs C++11 to work (more) properly."
#define CC11 "NO!"
#else
#define CC11 "yes"
#endif

#include "libmd.h"

#include "libmd/threads.libmd.cc"               //This file implements the thread structure
#include "libmd/autodiff.libmd.cc"              //This file implements automatic differentation
#include "libmd/potentials.libmd.cc"            //This file has all the builtin pairpotential functions
#include "libmd/particle.libmd.cc"              //This file implements the particle structure
#include "libmd/forces.libmd.cc"                //This file has all the builtin externalforces functions
#include "libmd/box.libmd.cc"                   //This file implements the box structure
#include "libmd/interact.libmd.cc"              //This file implements the interact structure
#include "libmd/index.libmd.cc"                 //This file takes care of indexing algorithms
#include "libmd/pairpotentials.libmd.cc"        //This file implements the pairpotential structure
#include "libmd/externalforces.libmd.cc"        //This file implements the externalforces structure
#include "libmd/integrators.libmd.cc"           //This file implements the integration structure
#include "libmd/variadic_vars.libmd.cc"         //This file implements the variadic_vars structure
#include "libmd/md.libmd.cc"                    //This file implements the md structure which takes care of molecular dynamics in flat space
#include "libmd/mongepatches.libmd.cc"          //This file has all the builtin monge patch functions and derivatives
#include "libmd/mp.libmd.cc"                    //This file implements the mp structure
#include "libmd/mpmd.libmd.cc"                  //This file implements the mpmd structure which takes care of molecular dynamics on monge patches

void __libmd__info()
{
    printf("libmd branch: %s\n",BRANCH);
    printf("libmd branch version: 0.%s\n",VER);
    printf("Compiler: %s\n",CC);
    printf("C++11: %s\n",CC11);
    printf("Compilation message: %s\n",CMSG);
}

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of LIBRARY source                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of LIBRARY source                                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __libmd_cc__
#define __libmd_cc__

#ifndef libmd_h
#warning '#include "libmd/libmd.cc"' is deprecated use '#include "libmd/libmd.h"'
#include "libmd.h"
#endif

void __libmd__info()
{
    //!
    //! This function is designed to give the user an overview of the compilation
    //!
    #ifdef BRANCH
    printf("libmd branch: " IO_BOLDCYAN "%s" IO_RESET "\n",BRANCH);
    #endif
    #ifdef VER
    printf("libmd branch version: " IO_BOLDCYAN "%s" IO_RESET "\n",VER);
    #endif
    #ifdef CC
    printf("Compiler: " IO_WHITE "%s" IO_RESET "\n",CC);
    #endif
    printf("Compiler version: " IO_WHITE "%s" IO_RESET "\n",__VERSION__);
    printf("Floating point type: ");
    #ifdef LIBMD__LONG_DOUBLE__
    printf(IO_BLUE "long double" IO_RESET "\n");
    #elif defined LIBMD__FLOAT__
    printf(IO_YELLOW "float" IO_RESET "\n");
    #else
    printf(IO_WHITE "double" IO_RESET "\n");
    #endif
    printf("Thread option: %s\n",THREAD_MODEL);
    #ifdef CMSG
    printf("Compilation message: " IO_YELLOW "%s" IO_RESET "\n",CMSG);
    #endif
}

#include "libmd-src/error.libmd.cc"                 //< This file implements the structure that handles errors/warnings/debug levels
#include "libmd-src/autodiff.libmd.cc"              //< This file implements automatic differentation
#include "libmd-src/potentials.libmd.cc"            //< This file has all the builtin pairpotential functions
#include "libmd-src/particle.libmd.cc"              //< This file implements the particle structure
#include "libmd-src/forces.libmd.cc"                //< This file has all the builtin externalforces functions
#include "libmd-src/box.libmd.cc"                   //< This file implements the box structure
#include "libmd-src/boundaryconditions.libmd.cc"    //< This file implements the periodicity functions
#include "libmd-src/bcond.libmd.cc"                 //< This file implements the periodicity structure
#include "libmd-src/interact.libmd.cc"              //< This file implements the interact structure
#include "libmd-src/index.libmd.cc"                 //< This file takes care of indexing algorithms
#include "libmd-src/pairpotentials.libmd.cc"        //< This file implements the pairpotential structure
#include "libmd-src/externalforces.libmd.cc"        //< This file implements the externalforces structure
#include "libmd-src/hooks.libmd.cc"                 //< This file implemente the hook structure
#include "libmd-src/integrators.libmd.cc"           //< This file implements the integration structure
#include "libmd-src/variadic_vars.libmd.cc"         //< This file implements the variadic_vars structure
#include "libmd-src/md.libmd.cc"                    //< This file implements the md structure which takes care of molecular dynamics in flat space
#include "libmd-src/mongepatches.libmd.cc"          //< This file has all the builtin monge patch functions and derivatives
#include "libmd-src/autodiff2.libmd.cc"             //< This file implements automatic differentation for Monge patches
#include "libmd-src/mp.libmd.cc"                    //< This file implements the mp structure
#include "libmd-src/mpmd.libmd.cc"                  //< This file implements the mpmd structure which takes care of molecular dynamics on monge patches

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of LIBRARY source                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

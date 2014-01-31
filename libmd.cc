///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of LIBRARY source                                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __LIBMD__
#define __LIBMD__

#define IO_RESET   "\033[0m"
#define IO_BLACK   "\033[30m"
#define IO_RED     "\033[31m"
#define IO_GREEN   "\033[32m"
#define IO_YELLOW  "\033[33m"
#define IO_BLUE    "\033[34m"
#define IO_MAGENTA "\033[35m"
#define IO_CYAN    "\033[36m"
#define IO_WHITE   "\033[37m"
#define IO_BOLDBLACK   "\033[1m\033[30m"
#define IO_BOLDRED     "\033[1m\033[31m"
#define IO_BOLDGREEN   "\033[1m\033[32m"
#define IO_BOLDYELLOW  "\033[1m\033[33m"
#define IO_BOLDBLUE    "\033[1m\033[34m"
#define IO_BOLDMAGENTA "\033[1m\033[35m"
#define IO_BOLDCYAN    "\033[1m\033[36m"
#define IO_BOLDWHITE   "\033[1m\033[37m"

#define STRING_ME(x) #x

#if __cplusplus < 201103L
#warning "C++11 not detetected: libmd needs C++11 to work (more) properly."
#define CC11 (IO_BOLDRED "NO!" IO_RESET)
#else
#define CC11 (IO_BOLDGREEN "yes \u263A" IO_RESET)
#endif

#ifdef THREADS
#define THREAD_MODEL (IO_BOLDYELLOW "C++11 STL" IO_RESET)
#elif OPENMP
#define THREAD_MODEL (IO_BOLDYELLOW "OpenMP" IO_RESET)
#else
#define THREAD_MODEL (IO_BOLDYELLOW "Threading disabled" IO_RESET)
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
    printf("libmd branch: " IO_BOLDCYAN "%s" IO_RESET "\n",BRANCH);
    printf("libmd branch version: " IO_BOLDCYAN "0.%s" IO_RESET "\n",VER);
    printf("Compiler: " IO_WHITE "%s" IO_RESET "\n",CC);
    printf("Compiler version: " IO_WHITE "%s" IO_RESET "\n",__VERSION__);
    printf("C++11: %s\n",CC11);
    printf("Thread option: %s\n",THREAD_MODEL);
    printf("Compilation message: " IO_YELLOW "%s" IO_RESET "\n",CMSG);
}

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of LIBRARY source                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

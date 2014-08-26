#ifndef libmd_h
#include "../libmd.h"
#endif

#define BUFFERSIZE 2048

#ifdef TIMER
ldf TicToc()
{
    static chrono::high_resolution_clock::time_point start,stop;
    start=stop,stop=chrono::high_resolution_clock::now();
    return chrono::duration<ldf>(stop-start).count();
}
#endif

/// This structure handles errors/warnings/debug levels
struct t_error
{
    char buffer[BUFFERSIZE];                                ///< Buffer of what needs to be printed (default size 2048 byte; static).
    ui term_level;                                          ///< Terminate level for libmd. The default value is 1.
    FILE *error_file;                                       ///< libmd error output file (default stderr)
    FILE *warning_file;                                     ///< libmd warning output file (default stderr)
    FILE *debug_1_file;                                     ///< libmd debug[1] output file (default stdout)
    FILE *debug_2_file;                                     ///< libmd debug[2] output file (default stdout)
    FILE *debug_3_file;                                     ///< libmd debug[3] output file (default stdout)
    FILE *debug_timer_file;                                 ///< libmd debug[timer] output file (default stdout)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    t_error();                                              ///< Constructor
    ~t_error();                                             ///< Destructor (to close the files)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void set_error_file(const char *fname);                 ///< Sets the error output file
    void set_warning_file(const char *fname);               ///< Sets the warning output file
    void set_debug_1_file(const char *fname);               ///< Sets the debug[1] output file
    void set_debug_2_file(const char *fname);               ///< Sets the debug[2] output file
    void set_debug_3_file(const char *fname);               ///< Sets the debug[3] output file
    void set_debug_timer_file(const char *fname);           ///< Sets the debug[timer] output file
    void print_error();                                     ///< Prints a error to the error output file (for internal use)
    void print_warning();                                   ///< Prints a warning to the warning output file (for internal use)
    void print_debug_1();                                   ///< Prints debug[1] message to the debug[1] output file (for internal use)
    void print_debug_2();                                   ///< Prints debug[2] message to the debug[2] output file (for internal use)
    void print_debug_3();                                   ///< Prints debug[3] message to the debug[3] output file (for internal use)
    void print_debug_timer();                               ///< Prints debug[timer] message to the debug[timer] output file (for internal use)
    void terminate(ui term);                                ///< Terminate if termlevel allows it (for internal use)
} error;

#define MSG_ERROR IO_BOLDRED "libmd-error: " IO_RESET
#define MSG_WARNING IO_BOLDBLUE "libmd-warning: " IO_RESET
#define MSG_DEBUG_1 IO_BOLDYELLOW "libmd-debug[1]: " IO_RESET
#define MSG_DEBUG_2 IO_BOLDMAGENTA "libmd-debug[2]: " IO_RESET
#define MSG_DEBUG_3 IO_BOLDCYAN "libmd-debug[3]: " IO_RESET
#define MSG_DEBUG_T IO_BOLDGREEN "libmd-timer["

t_error::t_error()
{
    #ifdef FE
    #ifdef FE_ALL_EXCEPT
    feenableexcept(FE_ALL_EXCEPT);
    #elif defined FE_EXCEPT
    feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
    #elif defined FE_LOW_EXCEPT
    feenableexcept(FE_DIVBYZERO);
    #else
    feclearexcept(FE_ALL_EXCEPT);
    #endif
    #endif
    term_level=1;
    error_file=stderr;
    warning_file=stderr;
    debug_1_file=stdout;
    debug_2_file=stdout;
    debug_3_file=stdout;
    debug_timer_file=stdout;
}

t_error::~t_error()
{
    fclose(error_file);
    fclose(warning_file);
    fclose(debug_1_file);
    fclose(debug_2_file);
    fclose(debug_3_file);
    fclose(debug_timer_file);
}

void t_error::set_error_file(const char *fname)
{
    if(!strcmp(fname,"stdout")) error_file=stdout;
    else if(!strcmp(fname,"stderr")) error_file=stderr;
    else error_file=fopen(fname,"w");
}

void t_error::set_warning_file(const char *fname)
{
    if(!strcmp(fname,"stdout")) warning_file=stdout;
    else if(!strcmp(fname,"stderr")) warning_file=stderr;
    else warning_file=fopen(fname,"w");
}

void t_error::set_debug_1_file(const char *fname)
{
    if(!strcmp(fname,"stdout")) debug_1_file=stdout;
    else if(!strcmp(fname,"stderr")) debug_1_file=stderr;
    else debug_1_file=fopen(fname,"w");
}

void t_error::set_debug_2_file(const char *fname)
{
    if(!strcmp(fname,"stdout")) debug_2_file=stdout;
    else if(!strcmp(fname,"stderr")) debug_2_file=stderr;
    else debug_2_file=fopen(fname,"w");
}

void t_error::set_debug_3_file(const char *fname)
{
    if(!strcmp(fname,"stdout")) debug_3_file=stdout;
    else if(!strcmp(fname,"stderr")) debug_3_file=stderr;
    else debug_3_file=fopen(fname,"w");
}

void t_error::set_debug_timer_file(const char *fname)
{
    if(!strcmp(fname,"stdout")) debug_timer_file=stdout;
    else if(!strcmp(fname,"stderr")) debug_timer_file=stderr;
    else debug_timer_file=fopen(fname,"w");
}

void t_error::print_error()
{
    fputs(buffer,error_file);
    terminate(0);
}

void t_error::print_warning()
{
    fputs(buffer,warning_file);
    terminate(1);
}

void t_error::print_debug_1()
{
    fputs(buffer,debug_1_file);
    terminate(4);
}

void t_error::print_debug_2()
{
    fputs(buffer,debug_2_file);
    terminate(3);
}

void t_error::print_debug_3()
{
    fputs(buffer,debug_3_file);
    terminate(2);
}

void t_error::print_debug_timer()
{
    fputs(buffer,debug_timer_file);
    terminate(5);
}

void t_error::terminate(ui term)
{
    if(term<term_level) exit(EXIT_FAILURE);
}

#ifdef PASS_ERROR
#define ERROR(str,...) \
{\
    int n=snprintf(error.buffer,BUFFERSIZE,"%s" IO_BOLDWHITE "%s:%d " IO_RESET "in" IO_WHITE " %s: " IO_RESET,MSG_ERROR,__FILE__,__LINE__,__FUNCTION__);\
    snprintf(error.buffer+n,BUFFERSIZE-n,str,##__VA_ARGS__);\
    strcat(error.buffer,"\n");\
    error.print_error();\
}
#else
#define ERROR(str,...) ;
#endif

#ifdef PASS_WARNING
#define WARNING(str,...)\
{\
    int n=snprintf(error.buffer,BUFFERSIZE,"%s" IO_BOLDWHITE "%s:%d " IO_RESET "in" IO_WHITE " %s: " IO_RESET,MSG_WARNING,__FILE__,__LINE__,__FUNCTION__);\
    snprintf(error.buffer+n,BUFFERSIZE-n,str,##__VA_ARGS__);\
    strcat(error.buffer,"\n");\
    error.print_warning();\
}
#else
#define WARNING(...) ;
#endif

#if DEBUG_LEVEL>0
#define DEBUG_1(str,...)\
{\
    int n=snprintf(error.buffer,BUFFERSIZE,"%s" IO_BOLDWHITE "%s:%d " IO_RESET "in" IO_WHITE " %s: " IO_RESET,MSG_DEBUG_1,__FILE__,__LINE__,__FUNCTION__);\
    snprintf(error.buffer+n,BUFFERSIZE-n,str,##__VA_ARGS__);\
    strcat(error.buffer,"\n");\
    error.print_debug_1();\
}
#else
#define DEBUG_1(str,...) ;
#endif

#if DEBUG_LEVEL>1
#define DEBUG_2(str,...)\
{\
    int n=snprintf(error.buffer,BUFFERSIZE,"%s" IO_BOLDWHITE "%s:%d " IO_RESET "in" IO_WHITE " %s: " IO_RESET,MSG_DEBUG_2,__FILE__,__LINE__,__FUNCTION__);\
    snprintf(error.buffer+n,BUFFERSIZE-n,str,##__VA_ARGS__);\
    strcat(error.buffer,"\n");\
    error.print_debug_2();\
}
#else
#define DEBUG_2(str,...) ;
#endif

#if DEBUG_LEVEL>2
#define DEBUG_3(str,...)\
{\
    int n=snprintf(error.buffer,BUFFERSIZE,"%s" IO_BOLDWHITE "%s:%d " IO_RESET "in" IO_WHITE " %s: " IO_RESET,MSG_DEBUG_3,__FILE__,__LINE__,__FUNCTION__);\
    snprintf(error.buffer+n,BUFFERSIZE-n,str,##__VA_ARGS__);\
    strcat(error.buffer,"\n");\
    error.print_debug_3();\
}
#else
#define DEBUG_3(str,...) ;
#endif

#if TIMER
#define DEBUG_TIMER(str,...)\
{\
    int n=snprintf(error.buffer,BUFFERSIZE,"%s%.10Lf]: " IO_RESET IO_BOLDWHITE "%s:%d " IO_RESET "in" IO_WHITE " %s: " IO_RESET,MSG_DEBUG_T,TicToc(),__FILE__,__LINE__,__FUNCTION__);\
    snprintf(error.buffer+n,BUFFERSIZE-n,str,##__VA_ARGS__);\
    strcat(error.buffer,"\n");\
    error.print_debug_timer();\
}
#else
#define DEBUG_TIMER(str,...) ;
#endif

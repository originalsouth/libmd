#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

ldf TicToc()
{
    static std::chrono::high_resolution_clock::time_point start,stop;
    start=stop,stop=std::chrono::high_resolution_clock::now();
    return std::chrono::duration<ldf>(stop-start).count();
}

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
    TicToc();
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

void t_error::print_error(char *buffer)
{
    fprintf(error_file,"%s\n",buffer);
    terminate(0);
}

void t_error::print_warning(char *buffer)
{
    fprintf(warning_file,"%s\n",buffer);
    terminate(1);
}

void t_error::print_debug_1(char *buffer)
{
    fprintf(debug_1_file,"%s\n",buffer);
    terminate(4);
}

void t_error::print_debug_2(char *buffer)
{
    fprintf(debug_2_file,"%s\n",buffer);
    terminate(3);
}

void t_error::print_debug_3(char *buffer)
{
    fprintf(debug_3_file,"%s\n",buffer);
    terminate(2);
}

void t_error::print_debug_timer(char *buffer)
{
    fprintf(debug_timer_file,"%s\n",buffer);
    terminate(5);
}

void t_error::terminate(ui term)
{
    if(term<term_level) exit(EXIT_FAILURE);
}

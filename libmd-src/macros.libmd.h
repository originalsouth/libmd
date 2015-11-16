#ifndef macros_h
#define macros_h

#define SYS ((md<dim>*)sys)

#ifndef BORING
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
#else
#define IO_RESET
#define IO_BLACK
#define IO_RED
#define IO_GREEN
#define IO_YELLOW
#define IO_BLUE
#define IO_MAGENTA
#define IO_CYAN
#define IO_WHITE
#define IO_BOLDBLACK
#define IO_BOLDRED
#define IO_BOLDGREEN
#define IO_BOLDYELLOW
#define IO_BOLDBLUE
#define IO_BOLDMAGENTA
#define IO_BOLDCYAN
#define IO_BOLDWHITE
#endif

#define BUFFERSIZE 2048

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

#ifdef TIMER
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

#define MSG_ERROR IO_BOLDRED "libmd-error: " IO_RESET
#define MSG_WARNING IO_BOLDBLUE "libmd-warning: " IO_RESET
#define MSG_DEBUG_1 IO_BOLDYELLOW "libmd-debug[1]: " IO_RESET
#define MSG_DEBUG_2 IO_BOLDMAGENTA "libmd-debug[2]: " IO_RESET
#define MSG_DEBUG_3 IO_BOLDCYAN "libmd-debug[3]: " IO_RESET
#define MSG_DEBUG_T IO_BOLDGREEN "libmd-timer["

#define THREAD_MODEL (IO_BOLDYELLOW "disabled" IO_RESET)

#define STRING_ME(x) #x

#endif

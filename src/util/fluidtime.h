#ifndef TIME_H
#define TIME_H

#ifdef WIN32
#include "Windows.h"
#endif
#ifdef __linux__
#include <sys/time.h>
#endif
class Time{
    public:
        // returns milliseconds passed since 1970 (as long int?!)
        static long int getMilliseconds(){
            #ifdef WIN32

                    LONGLONG currentCount;
                    QueryPerformanceCounter((LARGE_INTEGER*)&currentCount);

                    LONGLONG frequency;
                    QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);

                    /*
                        Ich hab mal ner anderen Windows time function gesucht, weil ich nicht genau versteh, wie du das machst
                        Also falls es probleme beim casting gibt:

                        SYSTEMTIME time;
                        GetSystemTime(&time);
                        LONG time_ms = (time.wSecond * 1000) + time.wMilliseconds;
                    */

                    return (double)currentCount / frequency;
            #endif

            #ifdef __linux__
                struct timeval tp;
                gettimeofday(&tp, NULL);

                long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
                return ms;
            #endif

            return 0;
        }

};


#endif

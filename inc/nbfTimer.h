#pragma once

/** @file nbfTimer.h
*	Timer.
*/

#include <time.h>

/** Measure execution time in seconds.	
*/
class nbfTimer {

 clock_t start_time, end_time;

public:

 void start()
   {
     start_time = clock();
   }

 void stop()
   {
     end_time = clock();
   }
 
 double elapsedSeconds()
   {
     return ( (double) (end_time - start_time) / CLOCKS_PER_SEC );
   }

};
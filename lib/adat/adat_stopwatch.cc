/*! @file
 * @brief Timer support
 *
 * A stopwatch like timer.
 */


#include <stdlib.h>
#include <iostream>
#include "adat/adat_stopwatch.h"

namespace Util
{

  StopWatch::StopWatch() 
  {
    stoppedP=false;
    startedP=false;
  }

  StopWatch::~StopWatch() {}

  void StopWatch::reset() 
  {
    startedP = false;
    stoppedP = false;
  }

  void StopWatch::start() 
  {
    int ret_val = gettimeofday(&t_start, NULL);
    if( ret_val != 0 ) 
    {
      std::cerr << "Gettimeofday failed in StopWatch::start()" << std::endl;
      exit(1);
    }
    startedP = true;
    stoppedP = false;
  }

  void StopWatch::stop() 
  {
    if( !startedP ) 
    { 
      std::cerr << "Attempting to stop a non running stopwatch in StopWatch::stop()" << std::endl;
      exit(1);
    }

    int ret_val = gettimeofday(&t_end, NULL);
    if( ret_val != 0 ) 
    {
      std::cerr << "Gettimeofday failed in StopWatch::end()" << std::endl;
      exit(1);
    }
    stoppedP = true;
  }

  double StopWatch::getTimeInMicroseconds() 
  {
    long usecs=0;
    if( startedP && stoppedP ) 
    { 
      if( t_end.tv_sec < t_start.tv_sec ) 
      { 
	std::cerr << "Critical timer rollover" << std::endl;
	exit(1);
      }
      else 
      { 
	usecs = (t_end.tv_sec - t_start.tv_sec)*1000000;

	if( t_end.tv_usec < t_start.tv_usec ) 
	{
	  usecs -= 1000000;
	  usecs += 1000000+t_end.tv_usec - t_start.tv_usec;
	}
	else 
	{
	  usecs += t_end.tv_usec - t_start.tv_usec;
	}
      }
    }
    else 
    {
      std::cerr << "Either stopwatch not started, or not stopped" << std::endl;
      exit(1);
    }

    return (double)usecs;
  }
    
  double StopWatch::getTimeInSeconds()  
  {
    long secs=0;
    long usecs=0;
    if( startedP && stoppedP ) 
    { 
      if( t_end.tv_sec < t_start.tv_sec ) 
      { 
	std::cerr << "Critical timer rollover" << std::endl;
	exit(1);
      }
      else 
      { 
	secs = t_end.tv_sec - t_start.tv_sec;

	if( t_end.tv_usec < t_start.tv_usec ) 
	{
	  secs -= 1;
	  usecs = 1000000;
	}
	usecs += t_end.tv_usec - t_start.tv_usec;
      }
    }
    else 
    {
      std::cerr << "Either stopwatch not started, or not stopped" << std::endl;
      exit(1);
    }

    return (double)secs + ((double)usecs / 1e6);
  }


  void StopWatch::splitTime() 
  {
    if( !startedP ) 
    { 
      std::cerr << "Attempting to get a split time from a non running stopwatch in StopWatch::split()" << std::endl;
      exit(1);
    }

    int ret_val = gettimeofday(&t_split, NULL);
    if( ret_val != 0 ) 
    {
      std::cerr << "Gettimeofday failed in StopWatch::splitTime()" << std::endl;
      exit(1);
    }
  }

  double StopWatch::getSplitTimeInMicroseconds() 
  {
    splitTime();

    long usecs=0;
    if( startedP && (! stoppedP) ) 
    { 
      if( t_split.tv_sec < t_start.tv_sec ) 
      { 
	std::cerr << "Critical timer rollover" << std::endl;
	exit(1);
      }
      else 
      { 
	usecs = (t_split.tv_sec - t_start.tv_sec)*1000000;

	if( t_split.tv_usec < t_start.tv_usec ) 
	{
	  usecs -= 1000000;
	  usecs += 1000000+t_split.tv_usec - t_start.tv_usec;
	}
	else 
	{
	  usecs += t_split.tv_usec - t_start.tv_usec;
	}
      }
    }
    else 
    {
      std::cerr << "Either stopwatch not started, or not stopped" << std::endl;
      exit(1);
    }

    return (double)usecs;
  }
    
  double StopWatch::getSplitTimeInSeconds()  
  {
    splitTime();

    long secs=0;
    long usecs=0;
    if( startedP && (! stoppedP) ) 
    { 
      if( t_split.tv_sec < t_start.tv_sec ) 
      { 
	std::cerr << "Critical timer rollover" << std::endl;
	exit(1);
      }
      else 
      { 
	secs = t_split.tv_sec - t_start.tv_sec;

	if( t_split.tv_usec < t_start.tv_usec ) 
	{
	  secs -= 1;
	  usecs = 1000000;
	}
	usecs += t_split.tv_usec - t_start.tv_usec;
      }
    }
    else 
    {
      std::cerr << "Either stopwatch not started, or not stopped" << std::endl;
      exit(1);
    }

    return (double)secs + ((double)usecs / 1e6);
  }


} // namespace Util

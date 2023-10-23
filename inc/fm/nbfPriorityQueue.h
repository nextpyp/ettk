#ifndef FILE_nbfPriorityQueue
#define FILE_nbfPriorityQueue

#include <queue>

#include <nbfDataPoint.h>

template< class T, class C, class Cmp >
class nbfPriorityQueue : public std::priority_queue< T,C,Cmp>
{
public :
  void clear(){ c.clear();};
explicit nbfPriorityQueue() : std::priority_queue< T,C,Cmp>(){};

};

template< class Pixel, int const Dim >
struct greaterPD : public binary_function< nbfDataPoint< Pixel, Dim >, nbfDataPoint< Pixel, Dim >, bool>
{
bool operator()( const nbfDataPoint< Pixel, Dim >& x, const nbfDataPoint< Pixel, Dim >& y) const
    {
      return ( *(x.distance) > *(y.distance) );
    }
};


#endif // FILE_nbfPriorityQueue

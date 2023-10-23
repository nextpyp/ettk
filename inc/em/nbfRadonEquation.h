#pragma once

#include <em/nbfRadonElement.h>
#include <vector>

/** Data structure to handle individual radon rays.
*/
template< class Pixel >
class nbfRadonEquation
{
public:

	// Default constructor
	nbfRadonEquation(){};

	// Constructor
	nbfRadonEquation( Pixel v ){ this->b = v; }

	~nbfRadonEquation(){};

	// Compute and iternally store projection along ray.
	inline void project();

	// Backproject error between current and given projections (times the given factor).
	inline void backProject( Pixel = 1.0 );

	// Store elements contributing to current ray.
	vector< nbfRadonElement< Pixel > > elements;

	// Store target projection.
	Pixel b;

	// Store norm of hyperplane normal.
	Pixel norm;

	// Store current value of projection.
	Pixel proj;

};


template< class Pixel >
inline void nbfRadonEquation< Pixel > :: project()
{
	this->proj = 0;
	typedef typename vector< nbfRadonElement< Pixel > > :: iterator myIterator;
	myIterator iter = this->elements.begin();
	while ( iter != this->elements.end() ){
		this->proj += iter->project();
		++iter;
	}
}

template< class Pixel >
inline void nbfRadonEquation< Pixel > :: backProject( Pixel factor )
{
	// skip if part of background
	if ( this->b == 0 ){
		return;
	}
	//Pixel update = ( this->b - this->proj ) * factor / this->elements.size();
	//Pixel update = ( this->b * this->norm - this->proj ) * factor / this->elements.size();
	Pixel update = ( this->b * this->norm - this->proj ) * factor / this->norm;
	//Pixel update = ( this->b - this->proj ) / this->norm * factor / (float)this->elements.size();
	//cout << update << endl;
	//if ( abs( update ) < 50 ){
	//	return;
	//}
	typedef typename vector< nbfRadonElement< Pixel > > :: iterator myIterator;
	myIterator iter = this->elements.begin();
	while ( iter != this->elements.end() ){
		iter->backProject( update );
		++iter;
	}
}
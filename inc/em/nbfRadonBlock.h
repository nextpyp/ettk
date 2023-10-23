#pragma once

#include <em/nbfRadonEquation.h>

/** Data structure to handle individual radon projections and backprojections.
*/
template< class Pixel >
class nbfRadonBlock
{
public:

	// constructor
	nbfRadonBlock( int e ){ this->equations.resize(e); }

	nbfRadonBlock( int, Array< Pixel, 1 > & );

	~nbfRadonBlock(){};

	void project();

	void backProject(Pixel=1.0);

	vector< nbfRadonEquation< Pixel > > equations;

	void setProjections( Array< Pixel, 1 > & );
};

template< class Pixel >
nbfRadonBlock< Pixel > :: nbfRadonBlock( int e, Array< Pixel, 1 > & A )
{
	this->equations.resize(e);
	this->setProjections( A );
}

template< class Pixel >
void nbfRadonBlock< Pixel > :: setProjections( Array< Pixel, 1 > & A )
{
	int e = this->equations.size();
	Pixel center = ( A.size() - 1.0 ) / 2.0 - ( e - 1.0 ) / 2.0;
	for ( int i = 0; i < e; i++ ){
		this->equations[i].b = A( ceil( center + i ) );
	}
}

template< class Pixel >
inline void nbfRadonBlock< Pixel > :: project()
{
	typedef typename vector< nbfRadonEquation< Pixel > > :: iterator myIterator;
	myIterator iter = this->equations.begin();
	while ( iter != this->equations.end() ){
		iter->project();
		++iter;
	}
}

template< class Pixel >
inline void nbfRadonBlock< Pixel > :: backProject(  Pixel factor )
{
	typedef typename vector< nbfRadonEquation< Pixel > > :: iterator myIterator;
	myIterator iter = this->equations.begin();
	while ( iter != this->equations.end() ){
		iter->backProject(factor);
		++iter;
	}
}
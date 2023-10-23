#pragma once

using namespace blitz;

#include <vector>
#include <em/nbfImageMetric.h>

/** Interface for Clustering Methods.
*/
template< class Pixel >
class nbfClustering
{
public:

	nbfClustering();

	~nbfClustering(){};

	/// Set input list of images.
	void setInput( vector< nbfWedgedSubImage3D< Pixel > > & in ){ this->input =  in; } 

	/// Set metric to use in computations. Reset state to recompute distance matrix.
	void setMetric( nbfCorrelationImageMetric< Pixel , 3 > * m ){ this->metric = m; }

	/// Set metric to use in computations. And set precomputed distance matrix;
	void setMetric( nbfCorrelationImageMetric< Pixel , 3 > * m, Array< Pixel, 4 > d )
	{ 
		this->metric = m; 
		this->alignments.resize( d.shape() );
		this->alignments = d;
		//this->alignments.reference( d );
	}

	/// Run clustering algorithm. Distance matrix is computed only once (expensive operation).
	void execute( Array< Pixel, 3 > & );

	void setReferences( Array< Pixel, 2 > & v ){ this->classes.resize( v.shape() ); this->classes = v;}

protected:

	/// Define interface for different clustering algorithms.
	virtual void doClustering( Array< Pixel, 3 > & ) = 0;

	vector< nbfWedgedSubImage3D< Pixel > > input;
	//vector< vtkImageData * > input;
	nbfCorrelationImageMetric< Pixel , 3 > * metric;

	Array< Pixel, 3 > classes;
	//vector< vector< int > > classes;

	// store distance, ccc, wedge overlap and alignment transformations
	Array< Pixel, 4 > alignments;
};

template< class Pixel >
nbfClustering< Pixel > :: nbfClustering()
{
	this->metric = NULL;
}

template< class Pixel >
void nbfClustering< Pixel > :: execute( Array< Pixel, 3 > & classification )
{
	// IGNORE - WE NO LONGER USE THE DISTANCE MATRIX
	if ( ( this->input.size() == 0 ) || ( this->alignments.rows() != this->input.size() ) ){
		// cerr << "WARNING - Dimensions of volume list and alignment matrix does not match: " << this->input.size() << " != " << this->alignments.rows() << endl;
	}
	// assume distance matrix is up to date
	this->doClustering( classification );
}
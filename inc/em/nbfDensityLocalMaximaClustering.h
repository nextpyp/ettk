#pragma once

using namespace blitz;

#include <em/nbfClustering.h>

/** Local Maximum Clustering Method.
*/
template< class Pixel >
class nbfDensityLocalMaximaClustering : public nbfClustering< Pixel >
{
public:

	nbfDensityLocalMaximaClustering();

	~nbfDensityLocalMaximaClustering(){};

	void setDistanceRadius( Pixel d ){ this->distanceRadius = d; }

protected:

	void doClustering( vector< vector< int > > & );

	Pixel distanceRadius;
};

template< class Pixel >
nbfDensityLocalMaximaClustering< Pixel > :: nbfDensityLocalMaximaClustering()
: nbfClustering< Pixel >()
{
	this->distanceRadius = .5e-4;
}

template< class Pixel >
void nbfDensityLocalMaximaClustering< Pixel > :: doClustering( vector< vector< int > > & classes )
{
	// assume distance matrix already computed
	classes.clear();

	// compute number of neighbors within distanceRadius.
	secondIndex j;
	Array< int, 1 > kNearestNeighbors( this->distanceMatrix.rows() );
	kNearestNeighbors = count( this->distanceMatrix < this->distanceRadius, j );

	vector< int > nearestNeighbors;

	// store flat maxima regions for grouping them together later on
	vector< TinyVector< int, 2 > > repeated;

	// look for local maxima of "density" function and store nearest neighbors in clusters.
	for ( int i = 0; i < this->distanceMatrix.rows(); i++ ){
		bool isLocalMaxima = true;
		nearestNeighbors.clear();
		// store current point as first cluster element.
		nearestNeighbors.push_back(i);
		for ( int j = 0; j < this->distanceMatrix.cols(); j++ ){
			if ( i != j ){
				// if adjacent neighbor
				if ( this->distanceMatrix(i,j) < this->distanceRadius ){
					nearestNeighbors.push_back(j);
					if ( kNearestNeighbors(j) > kNearestNeighbors(i) ){
						isLocalMaxima = false;
						continue;
					}
					else{
						// handle multiple local maxima
						if ( kNearestNeighbors(j) == kNearestNeighbors(i) ) {
							// check if already included
							bool alreadyIncluded = false;
							for ( int k = 0; k < repeated.size(); k++ ){
								if ( ( repeated[k](firstDim) == j ) && ( repeated[k](secondDim) == i ) ){
									alreadyIncluded = true;
									continue;
								}
							}
							// add if not already included
							if ( alreadyIncluded == false ){
								repeated.push_back( TinyVector< int, 2 >(i,j) );
							}
						}
					}
				}
			}
		}
		// store cluster
		if ( isLocalMaxima == true ){
			classes.push_back( nearestNeighbors );
		}
	}

	// post-processing: take care of multiple local maxima (group clusters together)
	for ( int i = 0; i < repeated.size(); i++ ){
		vector< vector< int > > :: iterator iter = classes.begin();
		while ( iter != classes.end() ){
			if ( (*iter)[0] == repeated[i](secondDim) ){
				classes.erase( iter );
			}
			else{
				++iter;
			}
		}
	}
}
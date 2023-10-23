#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkImageReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageReslice.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkImageNonMaximumSuppression.h>
#include <vtkImageMathematics.h>
#include <vtkCylindricalTransform.h>
#include <vtkImageMedian3D.h>
#include <vtkImageButterworthHighPass.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>

#include <nbfTimer.h>

#include <vtkPNGReader.h>

#include <random/normal.h>

#define PIXEL double

using namespace blitz;

void main( int argc, char ** argv )
{
	Array< PIXEL, 2 > D;
	nbfMatlabReader reader;
	reader.setFileName( argv[1] );
	reader.read(D);

	PIXEL distanceRadius = atof(argv[2]);

	vector< vector< int > > classes;

	// compute number of neighbors within distanceRadius.
	secondIndex j;
	Array< int, 1 > kNearestNeighbors( D.rows() );
	kNearestNeighbors = count( D < distanceRadius, j );

	vector< int > nearestNeighbors;

	// store flat maxima regions for grouping them together later on
	vector< TinyVector< int, 2 > > repeated;

	// look for local maxima of "density" function and store nearest neighbors in clusters.
	for ( int i = 0; i < D.rows(); i++ ){
		bool isLocalMaxima = true;
		nearestNeighbors.clear();
		// store current point as first cluster element.
		nearestNeighbors.push_back(i);
		for ( int j = 0; j < D.cols(); j++ ){
			if ( i != j ){
				// if adjacent neighbor
				if ( D(i,j) < distanceRadius ){
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

	for ( int i = 0; i < classes.size(); i++ ){
		cout << "Head: " << classes[i][0] << ", elements = " << classes[i].size() << endl;
	}
}
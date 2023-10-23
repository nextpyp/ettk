#pragma once

#include <vtkTransform.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageMedian3D.h>

#include <vector>

#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfWedgedAverageImage3D.h>
#include <em/nbfImageFilter.h>
#include <em/nbfFourierFilter.h>

#define NBF_IMAGE_METRIC_SPHERICAL	0
#define NBF_IMAGE_METRIC_CORR		1
#define NBF_IMAGE_METRIC_FOURIER	2
#define NBF_IMAGE_METRIC_CORE		3
#define NBF_IMAGE_METRIC_ANGULAR	4
#define NBF_IMAGE_METRIC_PROJECTION 5

/** Interface for VTK-like input-output pipeline filters.
	Update state is kept internally so execution is only done when needed.
	User is responsible for changing the state when doing changes that affect the filter's output.
*/
template< class Pixel, const int Dim >
class nbfImageMetric
{
public:

	nbfImageMetric( nbfImageFilter< Pixel, Dim > * = NULL, nbfFourierFilter< Pixel, Dim > * = NULL );
	virtual ~nbfImageMetric();

	/// Set first input volume + relative wedge orientation.
	virtual void setInput1( vtkImageData * );
	virtual void setInput1( nbfWedgedImage3D< Pixel > * );
	vtkImageData * getInput1(){ return this->input1; }

	/// Set second input volume + relative wedge orientation.
	virtual void setInput2( vtkImageData * );
	virtual void setInput2( nbfWedgedImage3D< Pixel > * );
	vtkImageData * getInput2(){ return this->input2; }

	/// To be defined on childs.
	virtual void execute();

	/// Distance value between inputs (function of the correlation value).
	virtual Pixel getDistance();

	virtual int getId() = 0;

	/** Get value and location of cross correlation peak.
	*/
	void findCorrelationPeak();

	/// Compute distance matrix for input image sequence.
	void getMatrix( vector< vtkImageData * > &, Array< Pixel, 2 > & );
	void getMatrix( vector< nbfWedgedSubImage3D< Pixel > > &, Array< Pixel, 2 > & );
	void getMatrix( char *, vector< nbfWedgedSubImage3D< Pixel > > &, Array< Pixel, 2 > &, Array< Pixel, 2 > * = NULL );

	/// Get value of correlation (assuming execute() was already called).
	virtual Pixel getCorrelationPeak();

	/// Get value of correlation (assuming execute() was already called).
	virtual Pixel getCorrelationScale();

	nbfWedgedImage3D< Pixel > * wedgedInput1;
	nbfWedgedImage3D< Pixel > * wedgedInput2;

	virtual Pixel getWedgeOverlap( vtkTransform * ) = 0;

	nbfImageFilter< Pixel, Dim > * imageFilter;
	nbfFourierFilter< Pixel, Dim > * fourierFilter;


	bool isTransformValid( vtkTransform * );
	//// check is transform is within valid limits ( translation allowance, rotation allowance )
	//bool isTransformValid( vtkTransform *, Pixel, Pixel );

	// Search for rotations in restricted area of space. Most commonly, search is restricted in the Y angle
	// and left unrestricted in the two Z angles. Y search space is [0,180]. Default = 180 degrees;
	void setRotationSearchRestriction( Pixel a ){ this->restrictRotationSearch = blitz::extrema::max( a, 1.0 ); } // take a restriction of at least 1 degree
	Pixel getRotationSearchRestriction(){ return this->restrictRotationSearch; }
	void setTranslationSearchRestriction( Pixel a ){ this->restrictTranslationSearch = a; }
	Pixel getTranslationSearchRestriction(){ return this->restrictTranslationSearch; }
	
	bool getToComputeOverlapNormalizedDistances(){ return this->overlapNormalizedDistances; }
	void setToComputeOverlapNormalizedDistances( bool i ){ this->overlapNormalizedDistances = i; }
	void computeOverlapNormalizedDistancesOn(){ this->overlapNormalizedDistances = true; }
	void computeOverlapNormalizedDistancesOff(){ this->overlapNormalizedDistances = false; }

	bool getToUseMutualCorrelation(){ return this->useMutualCorrelation; }
	void setToUseMutualCorrelation( bool i ){ this->useMutualCorrelation = i; }
	void useMutualCorrelationOn(){ this->useMutualCorrelation = true; }
	void useMutualCorrelationOff(){ this->useMutualCorrelation = false; }

	void usingMissingWedgeCompensationOn(){ this->useMissingWedgeCompensation = true; }
	void usingMissingWedgeCompensationOff(){ this->useMissingWedgeCompensation = false; }
	bool getMissingWedgeCompensation(){ return this->useMissingWedgeCompensation; }
	void setMissingWedgeCompensation( bool use ){ this->useMissingWedgeCompensation = use; }

protected:

	/// Store value of correlation.
	Pixel correlationPeak;

	Pixel correlationScale;

	// input images
	vtkImageData * input1;
	vtkImageData * input2;

	vtkImageCast * castVtk;

	bool input1Changed, input2Changed;

	/// Store image normalization values to avoid unnecesary re-computations.
	Pixel norm1, norm2;

	Pixel restrictTranslationSearch;
	Pixel restrictRotationSearch;

	// set whether to normalize distances by the overlap between both images in fourier space
	bool overlapNormalizedDistances;

	bool useMutualCorrelation;

	bool computeOverlaps;

	bool useMissingWedgeCompensation;
};

//template< class Pixel, const int Dim > 
//vector< nbfWedgedSubImage3D< Pixel > * > nbfImageMetric< Pixel, Dim > :: list1 = vector< nbfWedgedSubImage3D< Pixel > * >();
//
//template< class Pixel, const int Dim > 
//vector< nbfWedgedSubImage3D< Pixel > * > nbfImageMetric< Pixel, Dim > :: list2 = vector< nbfWedgedSubImage3D< Pixel > * >();

template< class Pixel, const int Dim >
nbfImageMetric< Pixel, Dim > :: nbfImageMetric( nbfImageFilter< Pixel, Dim > * imageFilter, nbfFourierFilter< Pixel, Dim > * fourierFilter )
: input1Changed(false), input2Changed(false)
{
	this->input1 = vtkImageData::New();
	this->input2 = vtkImageData::New();
	this->wedgedInput1 = NULL;
	this->wedgedInput2 = NULL;
	this->castVtk = vtkImageCast::New();
	this->castVtk->SetOutputScalarTypeToDouble();
	this->imageFilter = imageFilter;
	this->fourierFilter = fourierFilter;

	this->restrictRotationSearch = 180;
	this->restrictTranslationSearch = numeric_limits< Pixel > :: max();

	this->overlapNormalizedDistances = true;
	this->useMutualCorrelation = false;

	this->computeOverlaps = false;
	this->useMissingWedgeCompensation = false;
}

template< class Pixel, const int Dim >
nbfImageMetric< Pixel, Dim > :: ~nbfImageMetric(){
	this->castVtk->Delete();
	this->input1->Delete();
	this->input2->Delete();
}

template< class Pixel, const int Dim >
void nbfImageMetric< Pixel, Dim > :: setInput1( vtkImageData * i1 )
{
	this->input1->DeepCopy( i1 );
	this->input1Changed = true;
}

//template< class Pixel, const int Dim >
//void nbfImageMetric< Pixel, Dim > :: setInput1( nbfSubVolume< Pixel > * i1 )
//{
//	i1->getSubVolume( this->input1 );
//	this->setInput1( this->input1 );
//}

template< class Pixel, const int Dim >
void nbfImageMetric< Pixel, Dim > :: setInput1( nbfWedgedImage3D< Pixel > * i1 )
{
	this->wedgedInput1 = i1;
	i1->getImage( this->input1 );
	this->input1Changed = true;
}

template< class Pixel, const int Dim >
void nbfImageMetric< Pixel, Dim > :: setInput2( vtkImageData * i2 )
{
	this->input2->DeepCopy( i2 );
	this->input2Changed = true;
}

//template< class Pixel, const int Dim >
//void nbfImageMetric< Pixel, Dim > :: setInput2( nbfSubVolume< Pixel > * i2 )
//{
//	i2->getSubVolume( this->input2 );
//	this->setInput2( this->input2 );
//}

template< class Pixel, const int Dim >
void nbfImageMetric< Pixel, Dim > :: setInput2( nbfWedgedImage3D< Pixel > * i2 )
{
	this->wedgedInput2 = i2;
	i2->getImage( this->input2 );
	this->input2Changed = true;
	//	Array< double, 3 > image;
	//	nbfVTKInterface::vtkToBlitzReference( this->input2, image );
	//	cout << mean(image) << endl;
	//this->setInput2( this->input2 );
}

template< class Pixel, const int Dim >
void nbfImageMetric< Pixel, Dim > :: getMatrix( vector< vtkImageData * > & A, Array< Pixel, 2 > & D )
{
	// output will be square distance matrix with unitary diagonal
	D.resize( A.size(), A.size() );
	for ( int i = 0; i < A.size(); i++ ){
		for ( int j = i + 1; j < A.size(); j++ ){
			this->setInput1( A[i] );
			this->setInput2( A[j] );
			D(i,j) = this->getDistance();
			// make symmetric
			D(j,i) = D(i,j);
		}
	}
	for ( int i = 0; i < A.size(); i++ ){
		D(i,i) = 0;
	}
}

template< class Pixel, const int Dim >
void nbfImageMetric< Pixel, Dim > :: getMatrix( char * fileName, vector< nbfWedgedSubImage3D< Pixel > > & A, Array< Pixel, 2 > & D, Array< Pixel, 2 > * wOverlap )
{
	stringstream fileName1;
	fileName1 << fileName << ".distances.matlab";

	nbfMatlabReader mreader;
	mreader.setFileName( fileName1.str().c_str() );
	mreader.read( D );

	stringstream fileName2;
	fileName2 << fileName << ".overlaps.matlab";

	mreader.setFileName( fileName2.str().c_str() );

	if ( wOverlap != NULL ){
		mreader.read( *wOverlap );
	}

	if ( wOverlap != NULL ){
		if ( sum( abs( wOverlap->shape() - D.shape() ) ) != 0 ){
			cerr << "Incompatible sizes for distance and overlap matrices." << endl;
			return;
		}
	}

	double N = D.rows();

	nbfMatlabWriter mw;

	if( N != A.size() ){
	
		// output will be square distance matrix with unitary diagonal
		cout << "Computing Distance Matrix" << endl;
		
		if ( ( N > 0 ) && ( N < A.size() ) ){
            Array< Pixel, 2 >  Daux( D.shape() );
			Daux = D;
			D.resize( A.size(), A.size() );
			D = -1;
			D( Range(0, N-1) , Range( 0 , N-1)) = Daux;

			if( wOverlap != NULL){
				// Daux = -1;
				Daux = (*wOverlap);
				wOverlap->resize( A.size(), A.size());
				(*wOverlap) = 1;
				(*wOverlap)( Range(0, N-1) , Range( 0 , N-1) ) = Daux;
			}

		} else{
			// if N is biger than A.size() or non positive -> recompute distance matrix
			N = 0;
			D.resize( A.size(), A.size() );
			D = -1;

			if( wOverlap != NULL){
				wOverlap->resize( A.size(), A.size());
				(*wOverlap) = 1;
			}
		}
	}

	for ( int i = 0; i < A.size(); i++ ){
		//for ( int j = max((double)(i + 1) , N); j < A.size(); j++ ){
		for ( int j = i + 1; j < A.size(); j++ ){
			if ( D(i,j) == - 1 ){
				this->setInput1( &A[i] );
				this->setInput2( &A[j] );
				D(i,j) = this->getDistance();
				// make symmetric
				D(j,i) = D(i,j);
				cout << "D(" << i << "," << j << ") = " << D(i,j) << ",\t";

				mw.setFileName( fileName1.str().c_str() );
				mw.write( D );

				if ( wOverlap != NULL ){
					(*wOverlap)(i,j) =  this->getWedgeOverlap();
					(*wOverlap)(j,i) = (*wOverlap)(i,j);
					cout << "O(" << i << "," << j << ") = " << (*wOverlap)(i,j) << endl;

					mw.setFileName( fileName2.str().c_str() );
					mw.write( *wOverlap );

				}
			}
		}
		//nbfMatlabWriter w;
		//stringstream filename;
		//filename << "" << i << ".matlab";
		//w.setFileName( filename.str().c_str() );
		//w.write( D(i,Range::all() ));
	}
	// fill-in diagonal
	for ( int i = 0; i < A.size(); i++ ){
			if ( D(i,i) == - 1 ){
				D(i,i) = 0;
			}
	}
}

template< class Pixel, int const Dim >
void nbfImageMetric< Pixel, Dim > :: execute()
{
	Array< double, Dim > S1;
	nbfVTKInterface::vtkToBlitzReference( this->input1, S1 );

	if ( this->input1Changed == true ){

		// apply filter in real space
		if ( this->imageFilter != NULL ){
			this->imageFilter->execute( this->input1 );
		}

		// apply filter in real space
		if ( this->fourierFilter != NULL ){
			this->fourierFilter->execute( S1 );
		}

		S1 = S1 - mean(S1);
		this->norm1 = sum( pow2( S1 ) );

		this->input1Changed = false;
	}

	Array< double, Dim > S2;
	nbfVTKInterface::vtkToBlitzReference( this->input2, S2 );

	if ( this->input2Changed == true ){

		// apply filter in real space
		if ( this->imageFilter != NULL ){
			this->imageFilter->execute( this->input2 );
		}

		// apply filter in real space
		if ( this->fourierFilter != NULL ){
			this->fourierFilter->execute( S2 );
		}

		S2 = S2 - mean(S2);
		this->norm2 = sum( pow2( S2 ) );

		this->input2Changed = false;
	}

	this->correlationPeak = ( sum( S1 * S2 ) / sqrt( this->norm1 * this->norm2 ) );
}

template< class Pixel, const int Dim >
Pixel nbfImageMetric< Pixel, Dim > :: getDistance()
{
	this->getCorrelationPeak();
	if ( this->correlationPeak == numeric_limits< Pixel > :: max() ){
		cout << "Distance computation failed." << endl;
		return -1;
	} else {
		Pixel epsilon = 0e0;
		Pixel factor = 1.0;
		return ( factor * fabs( 1.0 - this->correlationPeak ) + epsilon );
	}
}

template< class Pixel, const int Dim >
Pixel nbfImageMetric< Pixel, Dim > :: getCorrelationPeak(){
	if ( this->input1Changed || this->input2Changed ){
		this->execute();
	}
	return this->correlationPeak; 
}

template< class Pixel, const int Dim >
Pixel nbfImageMetric< Pixel, Dim > :: getCorrelationScale(){
	if ( this->input1Changed || this->input2Changed ){
		this->execute();
	}
	return this->correlationScale; 
}

template< class Pixel, const int Dim >
bool nbfImageMetric< Pixel, Dim > :: isTransformValid( vtkTransform * t )
{
	// check if translation and rotation allowances are meet

	// check translation first
	double shifts[3];
	t->GetPosition(shifts);

	//if ( vtkMath::Norm( shifts ) > this->restrictTranslationSearch ){
	if ( ( fabs(shifts[0]) > this->restrictTranslationSearch ) ||
		 ( fabs(shifts[1]) > this->restrictTranslationSearch ) ||
		 ( fabs(shifts[2]) > this->restrictTranslationSearch ) ){
		//cout << "translation not conforming = " << vtkMath::Norm( shifts ) << endl;
		return false;
	}

	// transform the intrinsic normal vector
	double pin[3];
	pin[0] = pin[1] = 0; pin[2] = 1;
	double pout[3];
	t->TransformPoint( pin, pout );
		
	// substract shift before computing angle
	pout[0] -= shifts[0];
	pout[1] -= shifts[1];
	pout[2] -= shifts[2];

	vtkMath::Normalize( pin );
	vtkMath::Normalize( pout );

	double angle = vtkMath::Dot( pin, pout );
	if ( angle >= 1 ){
		angle = 0;
	} else {
		angle = fabs( vtkMath::DegreesFromRadians( acos( angle ) ) );
	}

	//if ( angle > this->restrictRotationSearch ){
	//	cout << "rotation not conforming" << endl;
	//}

	// check if angle is within allowable limits
	return ( angle <= this->restrictRotationSearch );
}
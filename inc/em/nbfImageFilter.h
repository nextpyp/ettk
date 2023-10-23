#pragma once

#include <vtkImageGradientMagnitude.h>
#include <vtkImageMedian3D.h>
#include <vtkImageConstantPad.h>

template< class Pixel > class nbfWedgedAverageImage3D;
// #include <em/nbfWedgedAverageImage3D.h>

/** Interface for VTK-like input-output pipeline filters.
	Update state is kept internally so execution is only done when needed.
	User is responsible for changing the state when doing changes that affect the filter's output.
*/
template< class Pixel, const int Dim >
class nbfImageFilter
{
public:

	nbfImageFilter();
	~nbfImageFilter();

	/// Apply filter in-place.
	void execute( vtkImageData *, bool = false );
	void execute( Array< double, Dim > &, bool = false );

	/// Apply symmetry
	static void symmetrize( vtkImageData *, int fold );

	/// Use window to mask image points. Default: ON.
	void maskOn(){ this->useMask = true; this->useFileMask = false; this->usePolarMask = false; this->mask.free(); }
	void maskOff(){ this->useMask = false; }

	/// Use window to mask image points. Default: ON.
	void fileMaskOn(){ 
		if ( this->maskFile.str().size() > 0 ){
			this->useMask = false; this->useFileMask = true; this->usePolarMask = false; this->mask.free(); 
		}
	}
	void fileMaskOff(){ this->useMask = true; this->useFileMask = false; this->mask.free(); }

	/// Set external image mask
	void setMaskFile( stringstream & s, Pixel & th ){ 
		//this->useMask = false;
		//this->useFileMask = true;
		//this->usePolarMask = false;
		//this->mask.free();
		this->maskFile.clear();
		this->maskFile << s.str(); 
		this->maskFileTh = th;
	}

	/// Use one dimensional window for computation of CCC in polar coordinates (this will force peridicity in the \rho dimension). Default: OFF.
	void polarMaskOn(){ this->usePolarMask = true; this->useMask = false; this->useFileMask = false; this->mask.free(); }
	void polarMaskOff(){ this->usePolarMask = false; }

	/// Use gaussian filter on inputs. Default: OFF.
	void gaussianFilterOn( Pixel s ){ this->useGaussianFilter = true; this->gaussianFactor = s; }
	void gaussianFilterOff(){ this->useGaussianFilter = false; }

	/// Use edge filter on inputs. Default: OFF.
	void edgeFilterOn(){ this->useEdgeFilter = true; }
	void edgeFilterOff(){ this->useEdgeFilter = false; }

	/// Use window to mask image points. Default: OFF.
	void medianFilterOn( int size ){ this->useMedianFilter = size; }
	void medianFilterOff(){ this->useMedianFilter = 0; }

	/// Use padding to increase accuracy (slower). Default: OFF.
	void paddingOn( int p = 1 ){ this->padding = p; }
	void paddingOff(){ this->paddingOn(1); }
	int getPaddingFactor(){ return this->padding; }

	// void setMaskSize( Pixel inner, Pixel outer ){ this->innerMaskSize = inner; this->outerMaskSize = outer; }
	void getMaskSize( Pixel & x, Pixel & y, Pixel & z, Pixel & sigma, bool & cylinder ){
		x = this->maskx;
		y = this->masky;
		z = this->maskz;
		sigma = this->sigma;
		cylinder = this->cylindricalMask;
	}

	void setMaskSize( Pixel x, Pixel y, Pixel z, Pixel sigma, bool cylinder = false ){
		this->useMask = true; this->useFileMask = false;
		this->maskx = x;
		this->masky = y;
		this->maskz = z;
		this->sigma = sigma;
		this->cylindricalMask = cylinder;
	}

	/// Store window image
	Array< Pixel, Dim > mask;

	void softenMask( Array< Pixel, Dim > &, Pixel = 0 );

	stringstream maskFile;
	Pixel maskFileTh;

	/// Build window image using Blackman's formula.
	void buildMask( TinyVector< int, Dim > & );
	void buildMaskFromFile( TinyVector< int, Dim > & );

	Pixel getSigma(){ return this->sigma; }
	void setSigma( Pixel s ){ this->sigma = s; this->mask.free(); }

	Pixel squareMaskSize;

	void setSymmetryFactor( int sym ){
		this->symmetryFactor = sym;
	}
	
protected:

	/// Build one-dimensional Blackman window image for correlation in polar coordinates.
	void buildPolarMask( vtkImageData * );

	bool useMask, useFileMask;

	bool usePolarMask;

	bool useGaussianFilter;
	double gaussianFactor;

	bool useEdgeFilter;

	int useMedianFilter;

	vtkImageConstantPad * pad;

	int padding;

	Pixel maskx, masky, maskz, sigma, cylindricalMask;
	
	int symmetryFactor;

};


template< class Pixel, const int Dim >
nbfImageFilter< Pixel, Dim > :: nbfImageFilter()
: useMask(false), useFileMask(false), maskFileTh(0), usePolarMask(false), useGaussianFilter(false), gaussianFactor(2.0), useEdgeFilter(false), useMedianFilter(0), padding(1), maskx(0), masky(0), maskz(0), sigma(.05), squareMaskSize(6), symmetryFactor(0)
{
	this->pad = vtkImageConstantPad::New();
}

template< class Pixel, const int Dim >
nbfImageFilter< Pixel, Dim > :: ~nbfImageFilter(){
	this->pad->Delete();
}

template< class Pixel, const int Dim >
void nbfImageFilter< Pixel, Dim > :: execute( Array< double, Dim > & A, bool bypassMask )
{
	vtkImageData * data = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( A, data );
	this->execute( data, bypassMask );
	nbfVTKInterface::vtkToBlitz( data, A );
	data->Delete();
}

template< class Pixel, const int Dim >
void nbfImageFilter< Pixel, Dim > :: execute( vtkImageData * data, bool bypassMask )
{
	// get blitz reference to VTK data
	Array< double, Dim > image;
	bool notUsingReference = false;
	if ( data->GetScalarType() == VTK_DOUBLE ){
		nbfVTKInterface::vtkToBlitzReference( data, image );
	} else {
		notUsingReference = true;
		nbfVTKInterface::vtkToBlitz( data, image );
	}

	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");
	//w.write(image);

	// apply gaussian filter if applicable
	if ( this->useGaussianFilter == true ){
		vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();
		filter->SetDimensionality(Dim);
		filter->SetInput( data );
		filter->SetStandardDeviations( this->gaussianFactor, this->gaussianFactor, this->gaussianFactor );
		filter->SetRadiusFactors( this->gaussianFactor * .75, this->gaussianFactor * .75, this->gaussianFactor * .75 );
		filter->Update();
		Array< double, Dim > reference;
		if ( filter->GetOutput()->GetScalarType() == VTK_DOUBLE ){
			nbfVTKInterface::vtkToBlitzReference( filter->GetOutput(), reference );
		} else {
			nbfVTKInterface::vtkToBlitz( filter->GetOutput(), reference );
		}
		image = reference;
		filter->Delete();
	}

	// apply edge filter if applicable
	if ( this->useEdgeFilter == true ){

		vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();
		filter->SetInput( data );
		filter->SetRadiusFactors(.5,.5,.5);
		filter->Update();

		vtkImageGradientMagnitude * gradient = vtkImageGradientMagnitude::New();
		gradient->SetDimensionality(Dim);
		gradient->SetInput( filter->GetOutput() );
		gradient->Update();
		Array< double, Dim > reference;
		if ( gradient->GetOutput()->GetScalarType() == VTK_DOUBLE ){
			nbfVTKInterface::vtkToBlitzReference( gradient->GetOutput(), reference );
		} else {
			nbfVTKInterface::vtkToBlitz( gradient->GetOutput(), reference );
		}
		image = reference;
		gradient->Delete();
		
		filter->Delete();
	}

	// apply median filter if applicable
	if ( this->useMedianFilter > 0 ){
		vtkImageMedian3D * median = vtkImageMedian3D::New();
		median->SetKernelSize( this->useMedianFilter, this->useMedianFilter, this->useMedianFilter );
		median->SetInput( data );
		median->Update();
		Array< double, Dim > reference;
		if ( median->GetOutput()->GetScalarType() == VTK_DOUBLE ){
			nbfVTKInterface::vtkToBlitzReference( median->GetOutput(), reference );
		} else {
			nbfVTKInterface::vtkToBlitz( median->GetOutput(), reference );
		}
		image = reference;
		median->Delete();
	}

	Pixel meanValue = 0.0;

	// apply mask if applicable
	if ( ( ( this->useMask == true ) || ( this->useFileMask == true ) || ( this->usePolarMask == true ) ) ){
		
		// initialize mask if nedded
		if ( data->GetNumberOfPoints() != this->mask.size() ){

			if ( this->useFileMask == true ){
				// try loading external mask (if available)
				TinyVector< int, Dim > s = image.shape();
				this->buildMaskFromFile( s );
			}
            
			if ( this->useMask == true ){
				TinyVector< int, Dim > s = image.shape();
				this->buildMask( s );
			}

			if ( this->usePolarMask == true ){
				this->buildPolarMask(data);
			}
		}

		if ( ( bypassMask == false ) && ( this->mask.size() == data->GetNumberOfPoints() ) ){
			// compute mean restricted to inside of mask
			//meanValue = sum( where( this->mask > 0, image, 0 ) ) / sum( where( this->mask > 0, 1, 0 ) );

			//float v = sum( pow2( image - meanValue ) ) / image.size() / 2.0;
			//cout << "m=" << meanValue << ", v=" << v << endl;
			//Pixel lower = ( meanValue - v );
			//Pixel upper = ( meanValue + v );
			//image = where( image < lower, lower, image );
			//image = where( image > upper, upper, image );
			meanValue = mean(image);

			// overwrite input with windowed data
			image = ( image - meanValue ) * this->mask + meanValue;
			//w.write(this->mask);
			//image = image * this->mask; // + meanValue;
			//w.write(image);
			//if ( B.size() > 0 ){
			//	B.resize( image.shape() );
			//	B = pow2( image - meanValue ) * this->mask + meanValue;
			//}
		}
	}

	//w.write(mask);

	// if not using a reference we need to update the input
	if ( notUsingReference == true ){
		nbfVTKInterface::blitzToVtk( image, data );
	}

	// PAD INPUT IMAGES
	if ( this->padding > 1 ){
		int extent[6];
		data->GetExtent(extent);
		if ( this->padding == 1 ){
			pad->SetOutputWholeExtent( extent[0], extent[1]+0, extent[2], extent[3]+0, extent[4], extent[5]+0 );
		} else if ( this->padding == 2 ){
			pad->SetOutputWholeExtent( extent[0], this->padding*extent[1]+1, extent[2], this->padding*extent[3]+1, extent[4], this->padding*extent[5]+1 );
		} else {
			cerr << "ERROR - Padding different than 2 not implemented." << endl;
		}
		pad->SetInput( data );
		pad->SetConstant( data->GetScalarComponentAsDouble(0,0,0,meanValue) );
		pad->Update();
		data->ShallowCopy( pad->GetOutput() );

		//if ( B.size() > 0 ){
		//	vtkImageData * dataSqr = vtkImageData :: New();
		//	nbfVTKInterface :: blitzToVtk( B, dataSqr );
		//	pad->SetInput( dataSqr );
		//	pad->SetConstant( dataSqr->GetScalarComponentAsDouble(0,0,0,meanValue) );
		//	pad->Update();
		//	nbfVTKInterface :: vtkToBlitz( pad->GetOutput(), B );
		//	dataSqr->Delete();
		//}
	}
}

template< class Pixel, const int Dim >
void nbfImageFilter< Pixel, Dim > :: buildMaskFromFile( TinyVector< int, Dim > & size )
{
	int length = this->maskFile.str().length();
	if ( length > 4 ){
		if ( ( this->maskFile.str().at( length - 4 ) == '.' ) &&
			( this->maskFile.str().at( length - 3 ) == 'm' ) &&
			( this->maskFile.str().at( length - 2 ) == 'r' ) &&
			( this->maskFile.str().at( length - 1 ) == 'c' ) ){
				nbfMrcReader r;
				r.setFileName( this->maskFile.str().c_str() );
				vtkImageData * data = vtkImageData :: New();
				r.read( data );
				nbfVTKInterface :: vtkToBlitz( data, this->mask );
				data->Delete();
			} else { // assume matlab format
				nbfMatlabReader r;
				r.setFileName( this->maskFile.str().c_str() );
				r.read( this->mask );
			}
			// assume mask is binary
			if ( sum( abs( this->mask.shape() - size ) ) == 0 ){
				// cerr << "Using image mask from file " << this->maskFile.str().c_str() << endl;
				this->mask = where( this->mask > this->maskFileTh, 1, 0 );
				this->softenMask( this->mask );
			} else {
				if ( this->mask.size() == 0 ){
					cerr << "Failed to read mask image from " << this->maskFile.str().c_str() << endl;
					this->mask.free();
				} else {
					cerr << "Specified mask: " << this->maskFile.str().c_str() << ", does not match current volume size." << endl;
					Array< Pixel, Dim > A( size );
					A = 0;
					int diff = ( size[0] - this->mask.rows() ) / 2.0;
					if ( Dim == 2 ){
						if ( diff > 0 ){
							Range I( size[0] / 2 + 1 - this->mask.rows() / 2 - 1, size[0] / 2 + 1 - this->mask.rows() / 2 - 1 + this->mask.rows() - 1 );
							A(I,I) = this->mask;
						} else {
							Range I( this->mask.rows() / 2 + 1 - size[0] / 2 - 1, this->mask.rows() / 2 + 1 - size[0] / 2 - 1 + size[0] - 1 );
							A = this->mask(I,I);
						}
					}
					if ( Dim == 3 ){
						if ( diff > 0 ) {
							Range I( size[0] / 2 + 1 - this->mask.rows() / 2 - 1, size[0] / 2 + 1 - this->mask.rows() / 2 - 1 + this->mask.rows() - 1 );
							A(I,I,I) = this->mask;
						} else {
							Range I( this->mask.rows() / 2 + 1 - size[0] / 2 - 1,this->mask.rows() / 2 + 1 - size[0] / 2 - 1 + size[0] - 1 );
							A = this->mask(I,I,I);
						}
					}
					this->mask.resize( A.shape() );
					this->mask = A;
					this->mask = where( this->mask > this->maskFileTh, 1, 0 );
					this->softenMask( this->mask );
					//nbfMatlabWriter w;
					//w.setFileName("p.matlab");
					//w.write(this->mask);
				}
			}
	}
}

template< class Pixel, const int Dim >
void nbfImageFilter< Pixel, Dim > :: buildMask( TinyVector< int, Dim > & size )
{
	this->mask.resize( size );
	
	TinyVector< Pixel, Dim > center = this->mask.shape() / 2.0;
	
	firstIndex in;
	secondIndex jn;
	thirdIndex kn;

	if ( this->maskx == 0 ){
		this->mask = 1;
		return;
	}
	// use circular mask
	//this->mask = pow2( in - center[0] ) + pow2( jn - center[1] ) + pow2( kn - center[2] ) < pow2( center[0] );

	if ( this->cylindricalMask == true ){
		// use cylindrical mask
			switch ( Dim ){
		case 2:
			this->mask = pow2( ( in - center[0] ) / this->maskx );
			break;
		case 3:
			this->mask = pow2( ( in - center[0] ) / this->maskx ) + pow2( ( jn - center[1] ) / this->masky );
			break;
		default:
			cerr << "ERROR - image mask not implemented for dimension other than 2 or 3"<< endl;
			}

			this->mask = where( this->mask <= 1, 1, 0 );
			switch ( Dim ){
				case 2:
					this->mask( Range :: all(), Range( 0, blitz :: extrema :: max( 0, center[1] - this->masky ) ) ) = 0;
					this->mask( Range :: all(), Range( blitz :: extrema :: min( center[1] - this->masky, this->mask.ubound(secondDim) ), this->mask.ubound(secondDim) ) ) = 0;
					break;
				case 3:
					this->mask( Range :: all(), Range :: all(), Range( 0, blitz :: extrema :: max( 0, center[2] - fabs(this->maskz) ) ) ) = 0;
					this->mask( Range :: all(), Range :: all(), Range( blitz :: extrema :: min( center[2] + fabs(this->maskz), this->mask.ubound(thirdDim) ), this->mask.ubound(thirdDim) ) ) = 0;
					break;
				default:
					cerr << "ERROR - image mask not implemented for dimension other than 2 or 3"<< endl;
			}

			// soften the mask
			this->softenMask( this->mask );
	} 
	else {
		// use spherical mask
		if ( ( this->maskx == this->masky ) && ( this->maskx == this->maskz ) ){
			// use raised cosine
			Pixel r2 = blitz :: extrema :: min( min( size ) - 1, ( this->maskx + sigma / 2.0 ) ); // outer radius
			Pixel r1 = blitz :: extrema :: max( 0.0, ( this->maskx - sigma / 2.0 ) ); // inner radius
			Pixel kappa = vtkMath::Pi() / ( r2 - r1 );

			switch ( Dim ){
		case 2:
			this->mask = sqrt( pow2( in - center[0] ) + pow2( jn - center[1] ) );
			break;
		case 3:
			this->mask = sqrt( pow2( in - center[0] ) + pow2( jn - center[1] ) + pow2( kn - center[2] ) );
			break;
		default:
			cerr << "ERROR - image mask not implemented for dimension other than 2 or 3"<< endl;
			}

			this->mask = where( this->mask <= r1, 1, this->mask );
			this->mask = where( this->mask >= r2, 0, this->mask );
			this->mask = where( ( this->mask > r1 ) && ( this->mask <= r2 ), ( 1 + cos( kappa * ( this->mask - r1 ) ) ) / 2.0, this->mask );
		} else {
			// use ellipsoidal mask
			switch ( Dim ){
		case 2:
			this->mask = pow2( ( in - center[0] ) / this->maskx ) + pow2( ( jn - center[1] ) / this->masky );
			break;
		case 3:
			this->mask = pow2( ( in - center[0] ) / this->maskx ) + pow2( ( jn - center[1] ) / this->masky ) + pow2( ( kn - center[2] ) / this->maskz );
			break;
		default:
			cerr << "ERROR - image mask not implemented for dimension other than 2 or 3"<< endl;
			}

			this->mask = where( this->mask <= 1, 1, 0 );

			// soften the mask
			this->softenMask( this->mask );
		}
	}
}

template< class Pixel, const int Dim >
void nbfImageFilter< Pixel, Dim > :: softenMask( Array< Pixel, Dim > & M, Pixel s )
{
	// TinyVector< Pixel, Dim > center = M.shape() / 2.0;
	if ( s == 0 ){
		s = this->sigma;
	}
	if ( s > 0 ){
		vtkImageData * window = vtkImageData::New();
		vtkImageGaussianSmooth * filter = vtkImageGaussianSmooth::New();
		nbfVTKInterface :: blitzToVtk( M, window );
		filter->SetInput( window );
		//Pixel factor = this->sigma * min( center );
		filter->SetStandardDeviations( s, s, s );
		filter->SetRadiusFactors( s * .75, s * .75, s * .75 );
		filter->Update();
		nbfVTKInterface :: vtkToBlitz( filter->GetOutput(), M );
		filter->Delete();
		window->Delete();
	}
}

template< class Pixel, const int Dim >
void nbfImageFilter< Pixel, Dim > :: buildPolarMask( vtkImageData * data )
{
	int dimensions[3];
	data->GetDimensions(dimensions);
	
	TinyVector< int, Dim > size;
	for ( int i = 0; i < Dim; i++ ){
		size[i] = dimensions[i];
	}
	this->mask.resize( size );
	
	TinyVector< Pixel, Dim > center = this->mask.shape() / 2.0;
	
	typename Array< Pixel, Dim > :: iterator iter = this->mask.begin();

	firstIndex in;
	secondIndex jn;
	thirdIndex kn;

	if ( Dim == 2 ){ // polar coordinates (apply window only in \rho direction)
		this->mask = where( fabs( in - center[0] ) < center[0], .42 + .5 * cos( vtkMath::Pi() * fabs( in - center[0] ) / center[0] ) + .08 * cos( 2 * vtkMath::Pi() * fabs( in - center[0] ) / center[0] ), 0 );
	}
	
	if ( Dim == 3 ){ // spherical coordinates (apply window in \rho and \phi directions)
		this->mask = where( fabs( in - center[0] ) < center[0], .42 + .5 * cos( vtkMath::Pi() * fabs( in - center[0] ) / center[0] ) + .08 * cos( 2 * vtkMath::Pi() * fabs( in - center[0] ) / center[0] ), 0 );
		this->mask *= where( fabs( kn - center[2] ) < center[2], .42 + .5 * cos( vtkMath::Pi() * fabs( kn - center[2] ) / center[2] ) + .08 * cos( 2 * vtkMath::Pi() * fabs( kn - center[2] ) / center[2] ), 0 );
	}
}


template< class Pixel, const int Dim >
void nbfImageFilter< Pixel, Dim > :: symmetrize( vtkImageData * input, int fold )
{
	if ( fold > 1 ){
		Array< double, 3 > result, A;
		vtkImageData * data = vtkImageData :: New();
		nbfVTKInterface :: vtkToBlitzReference( input, A );
		nbfWedgedAverageImage3D< Pixel > dummy;
		dummy.setAverageImage(A);
		// assume symmetry axis is Z direction
		for ( int i = 0; i < fold; i++ ){
			vtkTransform * t = vtkTransform :: New();
			t->RotateZ( i * 360.0 / fold );
			dummy.getImage( data, t );
			t->Delete();
			Array< double, 3 > current;
			nbfVTKInterface :: vtkToBlitzReference( data, current );
			if ( result.size() == 0 ){
				result.resize( current.shape() );
				result = 0;
			}
			result += current;
		}
		result /= fold;
		nbfVTKInterface :: blitzToVtk( result, input );
		data->Delete();
	}
}

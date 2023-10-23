#ifndef FILE_nbfMaximalFlow
#define FILE_nbfMaximalFlow

////////////////////////!!!!!!!!!
// REPLACE EDGEFILTER //!!!!!!!!!
////////////////////////!!!!!!!!!

#include <vtkExtractVOI.h>
#include <vtkImageResample.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageMathematics.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkMath.h>

#include <io/nbfVTKInterface.h>

#include <algorithm>
#include <vector>

template< class Pixel >
class nbfMaximalFlow
{
public:

	// constructor takes weight array as input
	nbfMaximalFlow(){};

	~nbfMaximalFlow(){};

	void execute(Array< Pixel, 2 > &, Array< Pixel, 2 > &, int );
	void execute(Array< Pixel, 3 > &, Array< Pixel, 3 > &, int );

	void weightMetric( Array< Pixel, 2 > &, Array< Pixel, 2 > & );
	void weightMetric( Array< Pixel, 3 > &, Array< Pixel, 3 > & );

protected:

	Pixel flood;

};

template< class Pixel >
void nbfMaximalFlow< Pixel > :: execute( Array< Pixel, 2 > & P, Array< Pixel, 2 > & G, int iterations )
{
	this->weightMetric(P,G);

	Array< Pixel, 2 > S( 2 * P.shape() + 1 );
	Pixel * dataZero = S.dataZero();
	S = 0;

	int rows = S.rows();
	int cols = S.cols();

	S( Range(1,rows-2,2), Range(1,cols-2,2) ) = P;

	// store only relevant positions
	vector< Pixel * > pPos;
	vector< Pixel * > FxPos;
	vector< Pixel * > FyPos;
	vector< Pixel > gPos;
	typename Array< Pixel, 2 > :: iterator iG = G.begin(); 
	// traverse P values in extended array S
	for ( int x = 1; x < rows - 1; x+=2 ){
		for ( int y = 1; y < cols - 1; y+=2 ){
			Pixel * current = dataZero + x * cols + y;
			if ( abs(*current) < 1 ){
				pPos.push_back( current );
				gPos.push_back( *iG );
				if ( x < rows - 2 ){ 
					FxPos.push_back( current + cols );
				}
				if ( y < cols - 2 ){
					FyPos.push_back( current + 1 );
				}
			}
			else{
				if ( ( x < rows - 2 ) && ( abs(*(current+2*cols)) < 1 ) ){
					FxPos.push_back( current + cols );
				}
				if ( ( y < cols - 2 ) && ( abs(*(current+2)) < 1 ) ){
					FyPos.push_back( current + 1 );
				}
			}
			++iG;
		}
	}

	typename vector< Pixel * > :: iterator iPos;
	typename vector< Pixel > :: iterator iGPos;

	Pixel dt = 1.0 / sqrt(2.0);

	// use large constant to ensure we have enough fuild to push through
	this->flood = 1000;

	P = this->flood * where( P == -1, 0, P );

	for ( int iteration = 0; iteration < iterations; iteration++ ){
	
		// update P
		iPos = pPos.begin();
		while ( iPos != pPos.end() ){
			(*(*iPos)) -= dt * ( *((*iPos)+cols) - *((*iPos)-cols) + *((*iPos)+1) - *((*iPos)-1) );
			++iPos;
		}

		// update Fx
		iPos = FxPos.begin();
		while ( iPos != FxPos.end() ){
			// update Fx
			(*(*iPos)) -= dt * ( *( (*iPos) + cols ) - *( (*iPos) - cols ) );
			++iPos;
		}

		// update Fy
		iPos = FyPos.begin();
		while ( iPos != FyPos.end() ){
			// update Fy
			(*(*iPos)) -= dt * ( *( (*iPos) + 1 ) - *( (*iPos) - 1 ) );
			++iPos;
		}

		// reinforce capacity constraint
		iPos = pPos.begin();
		iGPos = gPos.begin();
		Pixel PxMax, PyMax, norm, gVal;
		while ( iPos != pPos.end() ){
			gVal = *iGPos;
			Pixel * fMx = (*iPos) - cols;
			Pixel * fPx = (*iPos) + cols;
			Pixel * fMy = (*iPos) - 1;
			Pixel * fPy = (*iPos) + 1;
			PxMax = blitz::extrema::max( -(*fMx), blitz::extrema::max( (*fPx), 0.0f ) );
			PyMax = blitz::extrema::max( -(*fMy), blitz::extrema::max( (*fPy), 0.0f ) );
			norm = sqrt( PxMax * PxMax + PyMax * PyMax );
			if ( norm > gVal ){
				PxMax *= gVal / norm;
				PyMax *= gVal / norm;
				(*fMx) = blitz::extrema::max( (*fMx), - PxMax );
				(*fPx) = blitz::extrema::min( (*fPx),   PxMax );
				(*fMy) = blitz::extrema::max( (*fMy), - PyMax );
				(*fPy) = blitz::extrema::min( (*fPy),   PyMax );
			}
			++iPos;
			++iGPos;
		}
	}

	// extract useful part + re-scale wrt used constant
	P = S( Range(1,rows-2,2), Range(1,cols-2,2) ) / this->flood;
}


template< class Pixel >
void nbfMaximalFlow< Pixel > :: execute( Array< Pixel, 3 > & P, Array< Pixel, 3 > & G, int numIterations )
{
	//nbfTimer t;
	//t.start();
	this->weightMetric(P,G);
	//t.stop();
	//cout << "weighting time = " << t.elapsedSeconds() << endl;

	//nbfMatlabWriter writer;
	//writer.setFileName( "test.blitz" );
	//writer.write(P);

	//return;

	nbfVTKInterface converter;
#ifdef _DEBUG
	cout << P.shape() << endl;
#endif
	// Multiresolution

	int scale = 4;

	//Array< Pixel, 3 > P0( Pfull.shape() );

	//for ( int multres = 2; multres > 0; multres-- ){

		//Range I( 0, Pfull.rows() - 3, multres );
		//Range J( 0, Pfull.cols() - 3, multres );
		//Range K( 0, Pfull.depth() - 3, multres );

		//Array< Pixel, 3 > P( Pfull(I,J,K) );
		//Array< Pixel, 3 > G( Gfull(I,J,K) );

		//cout << P.shape() << endl;

		Array< Pixel, 3 > MW( P.shape() );
		MW = 0;

		Array< Pixel, 3 > S( 2 * P.shape() + 1 );
		Pixel * dataZero = S.dataZero();
		S = 0;

		int rows = S.rows();
		int cols = S.cols();
		int depth = S.depth();

		int offsetX = cols * depth;

		S( Range(1,rows-2,2), Range(1,cols-2,2), Range(1,depth-2,2) ) = P;

		// store only relevant positions
		vector< Pixel * > pPos;
		vector< Pixel * > FxPos;
		vector< Pixel * > FyPos;
		vector< Pixel * > FzPos;
		vector< Pixel > gPos;
		vector< Pixel * > sinkPos;

		typename Array< Pixel, 3 > :: iterator iG = G.begin(); 

		vector< Pixel * > MWPos;
		Pixel MG = mean( G );
		cout << "MG = [ " << min(G) << ", " << max(G) << ": " << mean(G) << "]" << endl;
		
		// traverse P values in extended array S
		for ( int x = 1; x < rows - 1; x+=2 ){
			for ( int y = 1; y < cols - 1; y+=2 ){
				for ( int z = 1; z < depth - 1; z+=2 ){
					Pixel * current = dataZero + x * offsetX + y * depth + z;
					Pixel * currentNew = MW.dataZero() + ( x - 1 ) / 2 * MW.cols() * MW.depth() + ( y - 1 ) / 2 * MW.depth() + ( z - 1 ) / 2;
					if ( abs(*current) < 1 ){
						pPos.push_back( current );
						gPos.push_back( *iG );
						MWPos.push_back( currentNew );
						if ( x < rows - 2 ){ 
							FxPos.push_back( current + offsetX );
						}
						if ( y < cols - 2 ){
							FyPos.push_back( current + depth );
						}
						if ( z < depth - 2 ){
							FzPos.push_back( current + 1 );
						}
					}
					else{
						if ( ( x < rows - 2 ) && ( abs(*(current+2*offsetX)) < 1 ) ){
							FxPos.push_back( current + offsetX );
						}
						if ( ( y < cols - 2 ) && ( abs(*(current+2*depth)) < 1 ) ){
							FyPos.push_back( current + depth );
						}
						if ( ( z < depth - 2 ) && ( abs(*(current+2)) < 1 ) ){
							FzPos.push_back( current + 1 );
						}
					}
					// store sink positions
					//if ( ( x == 3 ) || ( x == rows - 4 ) || ( y == 3 ) || ( y == cols - 4 ) || ( z == 3 ) || ( z == depth - 4 ) ){
					if ( ( x == 3 ) || ( x == rows - 4 ) || ( y == 3 ) || ( y == cols - 4 ) || ( z == 3 ) ){
						sinkPos.push_back( current );
					}
					++iG;
				}
			}
		}

		// set source = 10, sink = 0;
		//if ( multres == 1 ){
		//	P = P0;
		//}

		// use large constant to ensure we have enough fuild to push through
		this->flood = 1000;
	
		P = this->flood * where( P == -1, 0, P );
		S( Range(1,rows-2,2), Range(1,cols-2,2), Range(1,depth-2,2) ) = P;

		typename vector< Pixel * > :: iterator iPos;
		typename vector< Pixel > :: iterator iGPos;

		Pixel dt = 1.0 / sqrt(3.0);

		//int stopCount;
		//Pixel gamma = .1;
		//int usedElements = pPos.size();
#if 0
		cout << "iterations = " << numIterations << endl;
#endif
		bool counting = false;

		for ( int iteration = 0; iteration < numIterations; iteration++ ){

			// update P
			iPos = pPos.begin();
			while ( iPos != pPos.end() ){
				(*(*iPos)) -= dt * ( *((*iPos) + offsetX ) - *((*iPos) - offsetX ) + *((*iPos) + depth ) - *((*iPos) - depth ) + *((*iPos) + 1 ) - *((*iPos) - 1 ) );
				++iPos;
			}

			// update Fx
			iPos = FxPos.begin();
			while ( iPos != FxPos.end() ){
				// update Fx
				(*(*iPos)) -= dt * ( *( (*iPos) + offsetX ) - *( (*iPos) - offsetX ) );
				++iPos;
			}

			// update Fy
			iPos = FyPos.begin();
			while ( iPos != FyPos.end() ){
				// update Fy
				(*(*iPos)) -= dt * ( *( (*iPos) + depth ) - *( (*iPos) - depth ) );
				++iPos;
			}

			// update Fz
			iPos = FzPos.begin();
			while ( iPos != FzPos.end() ){
				// update Fy
				(*(*iPos)) -= dt * ( *( (*iPos) + 1 ) - *( (*iPos) - 1 ) );
				++iPos;
			}

			// reinforce capacity constraint
			iPos = pPos.begin();
			iGPos = gPos.begin();
			
			typename vector< Pixel * > :: iterator iMWPos;
			iMWPos = MWPos.begin();

			Pixel PxMax, PyMax, PzMax, gVal;
			while ( iPos != pPos.end() ){
				Pixel * fMx = (*iPos) - offsetX;
				Pixel * fPx = (*iPos) + offsetX;
				Pixel * fMy = (*iPos) - depth;
				Pixel * fPy = (*iPos) + depth;
				Pixel * fMz = (*iPos) - 1;
				Pixel * fPz = (*iPos) + 1;
				PxMax = max( -(*fMx), max( (*fPx), 0.0f ) );
				PyMax = max( -(*fMy), max( (*fPy), 0.0f ) );
				PzMax = max( -(*fMz), max( (*fPz), 0.0f ) );
				gVal = (*iGPos) / sqrt( PxMax * PxMax + PyMax * PyMax + PzMax * PzMax );
				
				// find out what the velocity field is
				Pixel fx = ( (*fPx) + (*fMx) ) / 2.0;
				Pixel fy = ( (*fPy) + (*fMy) ) / 2.0;
				Pixel fz = ( (*fPz) + (*fMz) ) / 2.0;
				Pixel orientation = 90;
				if ( fabs(fx) > 0 ){
					orientation = vtkMath::DegreesFromRadians(atan( fz / fx ));
				}
				*(*iMWPos) = orientation;

				// if ( ( iteration > 250 ) && ( fabs(orientation) > 80 ) ){
					// gVal = MG / sqrt( PxMax * PxMax + PyMax * PyMax + PzMax * PzMax );
				// }
				
				// if ( iteration == 400 ){
					// printf( "F = [ %02f, %02f, %02f ] = %02f\n", fx, fy, fz, *(*iMWPos) );
				// }

				// if ( ( gVal < 1 ) && ( fabs(orientation) < 60 ) ){
				if ( gVal < 1 ){
					PxMax *= gVal;
					PyMax *= gVal;
					PzMax *= gVal;
					(*fMx) = max( (*fMx), - PxMax );
					(*fPx) = min( (*fPx),   PxMax );
					(*fMy) = max( (*fMy), - PyMax );
					(*fPy) = min( (*fPy),   PyMax );
					(*fMz) = max( (*fMz), - PzMax );
					(*fPz) = min( (*fPz),   PzMax );
				}
				++iPos;
				++iGPos;
				++iMWPos;
			}
#ifdef _DEBUG
			if ( iteration % 25 == 0 ){
				Pixel accum = 0;
				iPos = sinkPos.begin();
				while ( iPos != sinkPos.end() ){
					accum += (*(*iPos));
					++iPos;
				}

				cout << "iter = " << iteration << ", net output = " << accum << " of " << sinkPos.size() <<endl;

				if ( ( !counting ) && ( accum > 0 ) ){
					counting = true;
				}

				//if ( counting && accum < 0 ){
				//	cout << iteration << " iterations needed." << endl;
				//	iteration = numIterations;
				//}
			}
			//cout << "Iter " << iteration << ", Contour size = " << saturated / 4.0 << endl;
			//// check stop condition
			//if ( ( iteration > 100 ) && ( stopCount > .98 * usedElements ) ){
			//	cout << iteration << " iterations needed." << endl;
			//	iteration = numIterations;
			//}
#else
			if ( ( iteration > 0 ) && ( iteration % 1000 == 0 ) ){
				cout << "Iteration " << iteration << " of " << numIterations << endl;
			}
#endif

		}
		//if ( multres == 2 ){
		//	Range I1( 0, Pfull.rows() - 2, multres );
		//	Range J1( 0, Pfull.cols() - 2, multres );
		//	Range K1( 0, Pfull.depth() - 2, multres );

		//	P0( I1, J1, K1 ) = where( P < .5, 0, 1 );

		//	P0( I1+1, J1+1, K1+1 ) = where( P < .5, 0, 1 );
		//}
		//else{
			P = S( Range(1,rows-2,2), Range(1,cols-2,2), Range(1,depth-2,2) ) / this->flood;
		//}
	//}
					nbfMrcWriter w;
					w.setFileName("iter_400.mrc");
					cout << "MW = [ " << min(MW) << ", " << max(MW) << "]" << endl;
					w.write(MW);
}

template< class Pixel >
void nbfMaximalFlow< Pixel > :: weightMetric( Array< Pixel, 3 > & P, Array< Pixel, 3 > & G )
{
	//nbfMatlabWriter writer;
	//writer.setFileName( "test.blitz" );

	// W = | \nabla \rho \conv \phi |

	// convolution is computed in Fourier domain as the product of FFTs
	// To avoid aliasing, the \rho image is downsampled a factor of 2

	int uboundX = floor( P.rows()  / 2.0 ) - 1;
	int uboundY = floor( P.cols()  / 2.0 ) - 1;
	int uboundZ = floor( P.depth()  / 2.0 ) - 1;
	Range I( 0, uboundX );
	Range J( 0, uboundY );
	Range K( 0, uboundZ );

	Array< Pixel, 3 > Pbig( P.shape() );
	Pbig = 0;

	// down sample by 2 to avoid aliasing
	Pbig( I, J, K ) = P( Range(0,2*uboundX,2), Range(0,2*uboundY,2), Range(0,2*uboundZ,2));
	
	// only consider source points
	Pbig = where( Pbig > 0, 1, 0 );

	//writer.write(Pbig);

	nbfVTKInterface converter;
	vtkImageData * data = vtkImageData::New();

	// convert to VTK
	converter.blitzToVtk( Pbig, data );

	// FFT( \rho )
	vtkImageFFT * fft1 = vtkImageFFT::New();
	fft1->SetDimensionality(3);
	fft1->SetInput( data );
	fft1->Update();

	// Build \phi = Green function in \Re^2
	Array< Pixel, 3 > Phi( P.shape() );
	Phi = 0;
	firstIndex i;
	secondIndex j;
	thirdIndex k;
	TinyVector< Pixel, 3 > center = Phi.shape() / 2.0;
	Pixel maxim = - sqrt( center(0) * center(0) + center(1) * center(1) + center(2) * center(2) ) / 100.0;
	Phi = maxim / ( sqrt( pow2( (i - center(0)) ) + pow2( (j - center(1)) ) + pow2( (k - center(2)) ) ) + 1e-1 );

	//writer.write(Phi);

	// convert to VTK
	vtkImageData * phiImage = vtkImageData::New();
	converter.blitzToVtk( Phi, phiImage );

	// FFT( \phi )
	vtkImageFFT * fft2 = vtkImageFFT::New();
	fft2->SetDimensionality(3);
	fft2->SetInput( phiImage );
	fft2->Update();

	// FFT( \rho ) * FFT( \phi )
	vtkImageMathematics * prod = vtkImageMathematics::New();
	prod->SetOperationToComplexMultiply();
	prod->SetInput1( fft1->GetOutput() );
	prod->SetInput2( fft2->GetOutput() );
	prod->Update();

	// go back to image domain
	vtkImageRFFT * rfft = vtkImageRFFT::New();
	rfft->SetDimensionality(3);
	rfft->SetInput( prod->GetOutput() );
	
	// correct FFT shift
	vtkImageFourierCenter * centerfft = vtkImageFourierCenter::New();
	centerfft->SetDimensionality(3);
	centerfft->SetInput( rfft->GetOutput() );

	// extract real component
	vtkImageExtractComponents * real = vtkImageExtractComponents::New();
	real->SetComponents(0);
	real->SetInput( centerfft->GetOutput() );
	real->Update();

	// keep only region of interest
	vtkExtractVOI * voi = vtkExtractVOI::New();
	voi->SetInput( real->GetOutput() );
	voi->SetVOI(0,P.rows()/2,0,P.cols()/2,0,P.depth()/2);
	voi->Update();

	int dimensions[3];
	voi->GetOutput()->GetDimensions(dimensions);

	// upsample a factor of 2 to go back to original size
	vtkImageResample * resa = vtkImageResample::New();
	resa->SetDimensionality(3);
	resa->SetInput( voi->GetOutput() );
	Pixel factorX = ( P.rows() + 1.0 ) / (Pixel)dimensions[0];
	Pixel factorY = ( P.cols() + 1.0 ) / (Pixel)dimensions[1];
	Pixel factorZ = ( P.depth() + 1.0 ) / (Pixel)dimensions[2];
	resa->SetAxisMagnificationFactor(0,factorX);
	resa->SetAxisMagnificationFactor(1,factorY);
	resa->SetAxisMagnificationFactor(2,factorZ);
	resa->SetInterpolationModeToCubic();
	resa->Update();

	// compute gradient magnitude
	vtkImageGradientMagnitude * gradMagnitude = vtkImageGradientMagnitude::New();
	gradMagnitude->SetDimensionality(3);
	//gradMagnitude->HandleBoundariesOn();
	gradMagnitude->SetInput( resa->GetOutput() );
	gradMagnitude->Update();

	//vtkStructuredPointsWriter * w = vtkStructuredPointsWriter::New();
	//w->SetFileName("w.vtk");
	//w->SetInput( gradMagnitude->GetOutput() );
	//w->Write();
	//w->Delete();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToFloat();
	cast->SetInput( gradMagnitude->GetOutput() );
	cast->Update();

	// apply weight to image
	converter.vtkToBlitz( cast->GetOutput(), Pbig );
	//cout << min(Pbig) << ", " << max(Pbig) << endl;
	G = G * Pbig;
	//writer.write(Pbig);
	//writer.write(G);

	cast->Delete();
	gradMagnitude->Delete();
	resa->Delete();
	voi->Delete();
	real->Delete();
	centerfft->Delete();
	rfft->Delete();
	prod->Delete();
	fft1->Delete();
	fft2->Delete();

}


template< class Pixel >
void nbfMaximalFlow< Pixel > :: weightMetric( Array< Pixel, 2 > & P, Array< Pixel, 2 > & G )
{
	//nbfMatlabWriter writer;
	//writer.setFileName( "test.blitz" );

	// W = | \nabla \rho \conv \phi |

	// convolution is computed in Fourier domain as the product of FFTs
	// To avoid aliasing, the \rho image is downsampled a factor of 2

	int uboundX = floor( P.rows()  / 2.0 ) - 1;
	int uboundY = floor( P.cols()  / 2.0 ) - 1;
	Range I( 0, uboundX );
	Range J( 0, uboundY );

	Array< Pixel, 2 > Pbig( P.shape() );
	Pbig = 0;

	// down sample by 2 to avoid aliasing
	Pbig( I, J ) = P( Range(0,2*uboundX,2), Range(0,2*uboundY,2) );
	
	// only consider source points
	Pbig = where( Pbig > 0, 1, 0 );

	//writer.write(Pbig);

	nbfVTKInterface converter;
	vtkImageData * data = vtkImageData::New();

	// convert to VTK
	converter.blitzToVtk( Pbig, data );

	// FFT( \rho )
	vtkImageFFT * fft1 = vtkImageFFT::New();
	fft1->SetDimensionality(2);
	fft1->SetInput( data );
	fft1->Update();

	// Build \phi = Green function in \Re^2
	Array< Pixel, 2 > Phi( P.shape() );
	Phi = 0;
	firstIndex i;
	secondIndex j;
	TinyVector< Pixel, 2 > center = Phi.shape() / 2.0;
	Pixel maxim = - sqrt( center(0) * center(0) + center(1) * center(1) ) / 100.0;
	Phi = maxim / ( sqrt( pow2( (i - center(0)) ) + pow2( (j - center(1)) ) ) + 1e-1 );

	//writer.write(Phi);

	// convert to VTK
	vtkImageData * phiImage = vtkImageData::New();
	converter.blitzToVtk( Phi, phiImage );

	// FFT( \phi )
	vtkImageFFT * fft2 = vtkImageFFT::New();
	fft2->SetDimensionality(2);
	fft2->SetInput( phiImage );
	fft2->Update();

	// FFT( \rho ) * FFT( \phi )
	vtkImageMathematics * prod = vtkImageMathematics::New();
	prod->SetOperationToComplexMultiply();
	prod->SetInput1( fft1->GetOutput() );
	prod->SetInput2( fft2->GetOutput() );
	prod->Update();

	// go back to image domain
	vtkImageRFFT * rfft = vtkImageRFFT::New();
	rfft->SetDimensionality(2);
	rfft->SetInput( prod->GetOutput() );
	
	// correct FFT shift
	vtkImageFourierCenter * centerfft = vtkImageFourierCenter::New();
	centerfft->SetDimensionality(2);
	centerfft->SetInput( rfft->GetOutput() );

	// extract real component
	vtkImageExtractComponents * real = vtkImageExtractComponents::New();
	real->SetComponents(0);
	real->SetInput( centerfft->GetOutput() );
	real->Update();

	// keep only region of interest
	vtkExtractVOI * voi = vtkExtractVOI::New();
	voi->SetInput( real->GetOutput() );
	voi->SetVOI(0,P.rows()/2,0,P.cols()/2,0,0);
	voi->Update();

	int dimensions[3];
	voi->GetOutput()->GetDimensions(dimensions);

	// upsample a factor of 2 to go back to original size
	vtkImageResample * resa = vtkImageResample::New();
	resa->SetDimensionality(2);
	resa->SetInput( voi->GetOutput() );
	Pixel factorX = ( P.rows() + 1.0 ) / (Pixel)dimensions[0];
	Pixel factorY = ( P.cols() + 1.0 ) / (Pixel)dimensions[1];
	resa->SetAxisMagnificationFactor(0,factorX);
	resa->SetAxisMagnificationFactor(1,factorY);
	resa->SetInterpolationModeToCubic();
	resa->Update();

	// compute gradient magnitude
	vtkImageGradientMagnitude * gradMagnitude = vtkImageGradientMagnitude::New();
	gradMagnitude->SetDimensionality(2);
	//gradMagnitude->HandleBoundariesOn();
	gradMagnitude->SetInput( resa->GetOutput() );
	gradMagnitude->Update();

	//vtkStructuredPointsWriter * w = vtkStructuredPointsWriter::New();
	//w->SetFileName("w.vtk");
	//w->SetInput( gradMagnitude->GetOutput() );
	//w->Write();
	//w->Delete();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToFloat();
	cast->SetInput( gradMagnitude->GetOutput() );
	cast->Update();

	// apply weight to image
	converter.vtkToBlitz( cast->GetOutput(), Pbig );
	//cout << min(Pbig) << ", " << max(Pbig) << endl;
	G = G * Pbig;
	//writer.write(Pbig);
	//writer.write(G);

	cast->Delete();
	gradMagnitude->Delete();
	resa->Delete();
	voi->Delete();
	real->Delete();
	centerfft->Delete();
	rfft->Delete();
	prod->Delete();
	fft1->Delete();
	fft2->Delete();
}

#endif /* FILE_nbfMaximalFlow */

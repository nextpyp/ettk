#ifndef FILE_nbf3DAlignment
#define FILE_nbf3DAlignment

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageReslice.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkBMPWriter.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageShrink3D.h>
#include <vtkImageEllipsoidSource.h>

#include <io/nbfVTKInterface.h>

using namespace blitz;

template< class Pixel >
class nbf3DAlignment
{
public:

	// constructor takes weight array as input
	nbf3DAlignment();

	~nbf3DAlignment();

	// set first image to align
	void setImage1( vtkImageData * i1, vtkTransform * t1 = NULL )
	{ 
		this->myImage1 = i1;
		this->t1 = t1;
	};

	// set second image to align
	void setImage2( vtkImageData * i2, vtkTransform * t2 = NULL )
	{
		this->myImage2 = i2;
		this->t2 = t2;
	};

	void setWedge( float );

	void execute( vtkImageData * );

	Pixel computeCCC( vtkImageData *, vtkImageData *, TinyVector< int, 3 > & );

	void normalizeFourier( vtkImageData * );

	void applyWindow( vtkImageData *, bool = false );
	void applyWindowComplex( vtkImageData * );

	void computeWedgeFilter( vtkTransform *, vtkImageData * );

protected:

	vtkImageData * myImage1;
	vtkImageData * myImage2;

	Array< Pixel, 3 > myWedge;
	vtkImageData * myWedgeVtk;
	vtkImageData * myCurrentWedge;

	vtkImageData * mask;

	vtkImageChangeInformation * change;
	vtkImageReslice * reslice;
	vtkImageFourierCenter * centerfft;

	vtkTransform * t1;
	vtkTransform * t2;

	vtkImageMathematics * math;
	vtkImageMathematics * finalConj;
	vtkImageMathematics * finalProd;
	vtkImageRFFT * ifft;
	vtkImageExtractComponents * real;
	vtkImageMathematics * norm;

	vtkImageMathematics * conj;	
	vtkImageMathematics * mult;	
	vtkImageMathematics * sqrt;	
	vtkImageShrink3D * shrink;
	vtkImageMathematics * prem;	

	// store canonical wedge position as plane normals
	float wnCan1[3], wnCan2[3];
};

template< class Pixel >
nbf3DAlignment< Pixel > :: nbf3DAlignment()
{
	this->myImage1 = NULL;
	this->myImage2 = NULL;
	this->myWedgeVtk = NULL;
	this->myCurrentWedge = vtkImageData::New();
	this->change = vtkImageChangeInformation::New();
	this->reslice = vtkImageReslice::New();
	this->reslice->WrapOff();
	this->reslice->MirrorOff();
	this->reslice->SetBackgroundLevel(0);

	this->centerfft = vtkImageFourierCenter::New();

	this->math = vtkImageMathematics::New();
	this->finalConj = vtkImageMathematics::New();
	this->finalProd = vtkImageMathematics::New();
	this->ifft = vtkImageRFFT::New();
	this->real = vtkImageExtractComponents::New();
	this->norm = vtkImageMathematics::New();

	this->conj = vtkImageMathematics::New();
	this->mult = vtkImageMathematics::New();
	this->shrink = vtkImageShrink3D::New();
	this->prem = vtkImageMathematics::New();	
	this->sqrt = vtkImageMathematics::New();	

	this->mask = vtkImageData::New();
}


template< class Pixel >
nbf3DAlignment< Pixel > :: ~nbf3DAlignment()
{
	if ( this->myWedgeVtk != NULL ){
		this->myWedgeVtk->Delete();
	}
	if ( this->myCurrentWedge != NULL ){
		this->myCurrentWedge->Delete();
	}
	this->change->Delete();
	this->reslice->Delete();
	this->centerfft->Delete();

	this->math->Delete();
	this->finalConj->Delete();
	this->finalProd->Delete();
	this->ifft->Delete();
	this->real->Delete();
	this->norm->Delete();

	this->conj->Delete();
	this->mult->Delete();
	this->shrink->Delete();
	this->prem->Delete();
	this->sqrt->Delete();

	this->mask->Delete();
}

//template< class Pixel >
//void nbf3DAlignment< Pixel > :: applyWindow( Array< Pixel, 3 > & in )
//{
//	Array< Pixel, 3 > window( in.shape() );
//	firstIndex i;
//	secondIndex j;
//	thirdIndex k;
//	TinyVector< int, 3 > center = floor( in.shape() / 2.0 );
//	window = blitz::sqrt( pow2( i - center[0] ) + 
//		                  pow2( j - center[1] ) + 
//				          pow2( k - center[2] ) + 0.0 );
//	window = where( window < floor( in.rows() / 2.0 ) - 1, window, floor( in.rows() / 2.0 ) - 1 );
//	Pixel alpha = 1.0 * max(window);
//	cout << min(window) << ", " << max(window) << endl;
//	window = .5 * ( 1 + cos( vtkMath::Pi() * window / alpha ) );
//	in = in * window;
//
//	nbfMatlabWriter w;
//	w.setFileName("ipath");
//	w.write(in);
//}

template< class Pixel >
void nbf3DAlignment< Pixel > :: applyWindow( vtkImageData * in, bool complex  )
{
	vtkImageFourierCenter * center = vtkImageFourierCenter::New();
	center->SetInput( this->mask );
	
	vtkImageMathematics * prod = vtkImageMathematics::New();
	prod->SetOperationToMultiply();
	prod->SetInput1( in );
	
	if ( complex == true ){
		center->Update();
		prod->SetInput2( center->GetOutput() );
	}
	else{
		prod->SetInput2( this->mask );
	}

	prod->Update();

	in->DeepCopy( prod->GetOutput() );

	prod->Delete();
	center->Delete();
}

template< class Pixel >
void nbf3DAlignment< Pixel > :: applyWindowComplex( vtkImageData * in )
{
	vtkImageData * data = vtkImageData::New();

	data->CopyStructure( in );
	data->SetNumberOfScalarComponents(1);
	data->SetScalarTypeToDouble();
	data->AllocateScalars();
	data->GetPointData()->GetScalars()->CopyComponent(0, in->GetPointData()->GetScalars(), 0 );

	this->applyWindow( data, true );
	in->GetPointData()->GetScalars()->CopyComponent(0, data->GetPointData()->GetScalars(), 0 );

	data->GetPointData()->GetScalars()->CopyComponent(0, in->GetPointData()->GetScalars(), 1 );
	this->applyWindow( data, true );
	in->GetPointData()->GetScalars()->CopyComponent(1, data->GetPointData()->GetScalars(), 0 );
	data->Delete();
}


template< class Pixel >
void nbf3DAlignment< Pixel > :: setWedge( float wedge )
{
	//nbfVTKInterface converter;
	//this->myWedge = vtkImageData::New();
	//this->myCurrentWedge = vtkImageData::New();
	//converter.blitzToVtk( in, this->myWedge );

	//// pre-center for rotation
	//this->change->SetInput( this->myWedge );
	//this->change->CenterImageOn();
	//this->change->Update();
	//this->myWedge->DeepCopy( change->GetOutput() );

	// check for valid wedge
	if ( ( wedge > 90 ) || ( wedge < 0 ) ){
		wedge = 90.0;
	}

	// store canonical wedge position as plane normals
	this->wnCan1[0] =   cos( ( 90.0 - wedge ) * vtkMath::DegreesToRadians() );
	this->wnCan1[1] =   0;
	this->wnCan1[2] = - sin( ( 90.0 - wedge ) * vtkMath::DegreesToRadians() );
	this->wnCan2[0] =   this->wnCan1[0];
	this->wnCan2[1] =   this->wnCan1[1];
	this->wnCan2[2] = - this->wnCan1[2];
}


// DUMMY
template< class Pixel >
void nbf3DAlignment< Pixel > :: execute( vtkImageData * output ){

	// create mask
	vtkImageEllipsoidSource * sphere = vtkImageEllipsoidSource::New();
	sphere->SetWholeExtent( this->myImage1->GetWholeExtent() );
	int dims[3];
	this->myImage1->GetDimensions(dims);
	sphere->SetCenter( floor( dims[0] / 2.0 ), floor( dims[1] / 2.0 ), floor( dims[2] / 2.0 ) );
	sphere->SetRadius( ceil( dims[0] / 2.1 ), ceil( dims[1] / 2.1 ), ceil( dims[2] / 2.1 ) );
	sphere->SetInValue(1.0);
	sphere->SetOutValue(0.0);
	sphere->SetOutputScalarTypeToDouble();

	vtkImageGaussianSmooth * gauss = vtkImageGaussianSmooth::New();
	gauss->SetDimensionality(3);
	gauss->SetRadiusFactors(5,5,5);
	gauss->SetInput( sphere->GetOutput() );
	gauss->Update();

	this->mask->SetScalarTypeToDouble();
	this->mask->ShallowCopy( gauss->GetOutput() );

	sphere->Delete();
	gauss->Delete();

	// initialize blitz wedge image
	this->myWedge.resize( dims[0], dims[1], dims[2] );

	vtkImageData * storeSecondImage = vtkImageData::New();
	storeSecondImage->DeepCopy( this->myImage2 );

	// BYPASS
	this->applyWindow( this->myImage1 );
	this->applyWindow( this->myImage2 );

	this->myWedgeVtk = vtkImageData::New();

	vtkImageData * imf1 = vtkImageData::New();
	vtkImageData * imf2 = vtkImageData::New();

	// pre-center for rotation
	this->change->SetInput( this->myImage2 );
	this->change->CenterImageOn();
	this->change->Update();
	// overwrite with centered copy
	this->myImage2->DeepCopy( change->GetOutput() );

	this->reslice->SetInput( this->myImage2 );

	TinyVector< int, 3 > translation;
	Pixel currentCCCmaxima = - numeric_limits< Pixel > :: max();

	float bestAngleX = 0;
	float bestAngleY = 0;
	float bestAngleZ = 0;

	vtkImageMathematics * math = vtkImageMathematics::New();
	vtkImageFFT * fft1 = vtkImageFFT::New();
	vtkImageFFT * fft2 = vtkImageFFT::New();

	this->myCurrentWedge->CopyStructure( this->myImage1 );
	this->myCurrentWedge->SetScalarTypeToDouble();
	this->myCurrentWedge->SetNumberOfScalarComponents(2);
	this->myCurrentWedge->AllocateScalars();

	for ( float angleX = 0; angleX <= 10; angleX++ ){
		for ( float angleY = 0; angleY <= 0; angleY++ ){
			for ( float angleZ = 0; angleZ <= 0; angleZ++ ){

				// concatenate transforms
				vtkTransform * concat = vtkTransform :: New();
				concat->RotateX(-angleX);
				concat->RotateY(-angleY);
				concat->RotateZ(-angleZ);

				// Compute intersection of both wedges
				// the first one is the canonical transformed by t1
				// the second is the canonical transformed by ( t2 concatenated with concat)
				this->computeWedgeFilter( concat, this->myWedgeVtk );

				// build wedge filter in complex format
				// vtkImageData * filter = vtkImageData::New();
				this->myCurrentWedge->GetPointData()->GetScalars()->CopyComponent(0, this->myWedgeVtk->GetPointData()->GetScalars(), 0 );
				this->myCurrentWedge->GetPointData()->GetScalars()->CopyComponent(1, this->myWedgeVtk->GetPointData()->GetScalars(), 0 );

				//vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
				//writer->SetFileName("wedge.vtk");
				//writer->SetInput( this->myWedgeVtk );
				//writer->Write();

				// shift FFT before applying filter
				this->centerfft->SetInput( this->myCurrentWedge );
				this->centerfft->Update();

				// overwrite with shifted version
				this->myCurrentWedge->ShallowCopy( this->centerfft->GetOutput() );

				// FFT of input image
				fft1->SetDimensionality(3);
				fft1->SetInput( this->myImage1 );
				fft1->Update();
				imf1->ShallowCopy( fft1->GetOutput() );

				// rotate template image
				this->reslice->SetResliceTransform( concat );
				this->reslice->Update();

				// FFT of input image
				fft2->SetDimensionality(3);
				fft2->SetInput( this->reslice->GetOutput() );
				fft2->Update();
				imf2->ShallowCopy( fft2->GetOutput() );

				TinyVector< int, 3 > maximaPosition;
				Pixel maxima = this->computeCCC( imf1, imf2, maximaPosition );
				if ( maxima > currentCCCmaxima ){
					currentCCCmaxima = maxima;
					translation = maximaPosition;
					bestAngleX = angleX;
					bestAngleY = angleY;
					bestAngleZ = angleZ;
				}

				cout << "theta = [" << angleX << "," << angleY << "," << angleZ <<"] , score = " << maxima << ", t = " << translation << endl;

				concat->Delete();
			}
		}
	}

	this->change->SetInput( storeSecondImage );
	this->reslice->SetInput( this->change->GetOutput() );

	// generate output aligned to first image
	if ( output != NULL ){
		vtkTransform * t = vtkTransform::New();
		t->RotateX(bestAngleX);
		t->RotateY(bestAngleY);
		t->RotateZ(bestAngleZ);
		t->Translate( translation[0], translation[1], translation[2] );
		this->reslice->SetResliceTransform( t );
		this->reslice->Update();
		output->DeepCopy( this->reslice->GetOutput() );
		t->Delete();
	}

	storeSecondImage->Delete();

	math->Delete();
	fft1->Delete();
	fft2->Delete();

	cout << "Alignment = [" << bestAngleX << "," << bestAngleY << "," << bestAngleZ <<"] , t = " << translation << endl;

	imf1->Delete();
	imf2->Delete();
}

template< class Pixel >
Pixel nbf3DAlignment< Pixel > :: computeCCC( vtkImageData * imf1, vtkImageData * imf2, TinyVector< int, 3 > & p )
{
	// apply wedge filter to I1
	this->math->SetOperationToComplexMultiply();
	this->math->SetInput1( this->myCurrentWedge );
	this->math->SetInput2( imf1 );
	this->math->Update();

	// bypass wedge filter
	imf1->DeepCopy( this->math->GetOutput() );

	// apply wedge filter to I2
	math->SetInput2( imf2 );
	math->Update();

	// bypass wedge filter
	imf2->ShallowCopy( this->math->GetOutput() );

	this->normalizeFourier( imf1 );
	this->normalizeFourier( imf2 );

	// conj(imf2)
	this->finalConj->SetOperationToConjugate();
	this->finalConj->SetInput1( imf2 );
	
	// imf1 .* conj(imf2)
	this->finalProd->SetOperationToComplexMultiply();
	this->finalProd->SetInput1( imf1 );
	this->finalProd->SetInput2( this->finalConj->GetOutput() );
	this->finalProd->Update();

	// BYPASS
	// this->applyWindowComplex( this->finalProd->GetOutput() );

	// ifft( imf1 .* conj(imf2) )
	this->ifft->SetDimensionality(3);
	this->ifft->SetInput( this->finalProd->GetOutput() );
	
	// ifftshift( ifft( imf1 .* conj(imf2) ) )
	this->centerfft->SetInput( this->ifft->GetOutput() );
	
	// real( ifftshift( ifft( imf1 .* conj(imf2) ) ) )
	this->real->SetInput( this->centerfft->GetOutput() );
	this->real->SetComponents(0);
	this->real->Update();

	// Save CCC volume
	//this->norm->SetOperationToMultiplyByK();
	//this->norm->SetConstantK(1);
	//this->norm->SetInput1( this->real->GetOutput() );

	//vtkStructuredPointsWriter * w = vtkStructuredPointsWriter::New();
	//w->SetFileName("maxima.vtk");
	//w->SetInput( norm->GetOutput() );
	//w->Write();

	// get maxima of CCC
	double scalarRange[2];
	this->real->GetOutput()->GetScalarRange( scalarRange );

	bool allDone = false;

	int size[3];
	this->real->GetOutput()->GetDimensions(size);

	// get maxima location
	for ( int i = 0; i < size[0]; i++ ){
		for ( int j = 0; j < size[1]; j++ ){
			for ( int k = 0; k < size[2]; k++ ){
				if ( this->real->GetOutput()->GetScalarComponentAsDouble(i,j,k,0) == scalarRange[1] ){
					p(firstDim)  = 1 + size[0] / 2 - i;
					p(secondDim) = 1 + size[1] / 2 - j;
					p(thirdDim)  = 1 + size[2] / 2 - k;
					allDone = true;
				}
				if ( allDone ){
					break;
				}
			}
			if ( allDone ){
				break;
			}
		}
		if ( allDone ){
			break;
		}
	}

	// more elaborate maxima selection

	// find all local maxima locations
	// sort acording to increasing CCC value
	// if first/second ~1 => not a good match
	// if translation too far away from center also discard

	Pixel firstBest = numeric_limits< Pixel > :: max();
	Pixel secondBest = numeric_limits< Pixel > :: max();
	
	for ( int i = 1; i < size[0] - 1; i++ ){
		for ( int j = 1; j < size[1] - 1; j++ ){
			for ( int k = 1; k < size[2] - 1; k++ ){
				if (



	//for ( int i = 0; i < real->GetOutput()->GetNumberOfPoints(); i++ ){
	//	if ( real->GetOutput()->GetPointData()->GetComponent(i,0) == scalarRange[1] ){
	//		double coords[3];
	//		coords = real->GetOutput()->GetPoint(i);
	//		p(firstDim)  = coords[0];
	//		p(secondDim) = coords[1];
	//		p(thirdDim)  = coords[2];
	//		break;
	//	}
	//}

	//real->Delete();
	//ifft->Delete();
	//finalProd->Delete();
	//finalConj->Delete();
	//math->Delete();

	return scalarRange[1];
}

// Assume complex input and compute normalized version for correlation computation
// This is an in-place filter, overwrites the input
template< class Pixel >
void nbf3DAlignment< Pixel > :: normalizeFourier( vtkImageData * im )
{
	// kill DC value
	im->SetScalarComponentFromDouble(0,0,0,0,0);
	im->SetScalarComponentFromDouble(0,0,0,1,0);

	// math pipeline
	// conj(im)
	//vtkImageMathematics * conj = vtkImageMathematics::New();	
	this->conj->SetOperationToConjugate();
	this->conj->SetInput1( im );
	this->conj->Update();

	// im .* conj(im)
	//vtkImageMathematics * mult = vtkImageMathematics::New();	
	this->mult->SetOperationToComplexMultiply();
	this->mult->SetInput1( im );
	this->mult->SetInput2( this->conj->GetOutput() );
	this->mult->Update();

	//// sqrt( im .* conj(im) )
	////vtkImageMathematics * sqrt = vtkImageMathematics::New();	
	//this->sqrt->SetOperationToSquareRoot();
	//this->sqrt->SetInput1( this->mult->GetOutput() );
	//this->sqrt->Update();

	// sum( sqrt( im .* conj(im) ) )
	//vtkImageShrink3D * shrink = vtkImageShrink3D::New();
	int size[3];
	im->GetDimensions(size);
	this->shrink->SetShrinkFactors( size );
	this->shrink->MeanOn();
	this->shrink->SetInput( this->mult->GetOutput() );
	this->shrink->Update();

	// im / sum( sqrt( im .* conj(im) ) )
	//vtkImageMathematics * prem = vtkImageMathematics::New();	
	this->prem->SetOperationToMultiplyByK();
	double constantK = 1.0f / sqrtf( this->shrink->GetOutput()->GetScalarComponentAsDouble(0,0,0,0) );
	this->prem->SetConstantK( constantK );
	this->prem->SetInput1( im );
	this->prem->Update();

	// store result in original image
	im->ShallowCopy( this->prem->GetOutput() );

	//prem->Delete();
	//shrink->Delete();
	//sqrt->Delete();
	//mult->Delete();
	//conj->Delete();
}


// Assume complex input and compute normalized version for correlation computation
// This is an in-place filter, overwrites the input
template< class Pixel >
void nbf3DAlignment< Pixel > :: computeWedgeFilter( vtkTransform * t, vtkImageData * out )
{
	firstIndex i; secondIndex j; thirdIndex k;

	// apply t1 to canonical vectors
	float wnCanT1[3], wnCanT2[3];
	if ( this->t1 != NULL ){
		this->t1->TransformVector(this->wnCan1,wnCanT1);
		this->t1->TransformVector(this->wnCan2,wnCanT2);
	}
	else{
		wnCanT1[0] = this->wnCan1[0];
		wnCanT1[1] = this->wnCan1[1];
		wnCanT1[2] = this->wnCan1[2];
		wnCanT2[0] = this->wnCan2[0];
		wnCanT2[1] = this->wnCan2[1];
		wnCanT2[2] = this->wnCan2[2];
	}

	// apply ( t2 * t ) to canonical vectors
	vtkTransform * transform = vtkTransform::New();
	if ( this->t2 != NULL ){
		transform->DeepCopy( this->t2 );
		transform->Concatenate( t );
	}
	else{
		transform->DeepCopy( t );
	}

	float nt1[3], nt2[3];
	transform->TransformVector(this->wnCan1,nt1);
	transform->TransformVector(this->wnCan2,nt2);

	int dims[3];
	this->myImage1->GetDimensions(dims);
	float center[3];
	center[0] = ( dims[0] - 1.0 ) / 2.0;
	center[1] = ( dims[1] - 1.0 ) / 2.0;
	center[2] = ( dims[2] - 1.0 ) / 2.0;

	// APPLY: t1 ( wedge ) & [ t2 * t ]( wedge )
	this->myWedge = ( ( wnCanT1[0] * (i - center[0] ) + 
		                wnCanT1[1] * (j - center[1] ) + 
					    wnCanT1[2] * (k - center[2] ) ) *
			          ( wnCanT2[0] * (i - center[0] ) + 
		                wnCanT2[1] * (j - center[1] ) + 
					    wnCanT2[2] * (k - center[2] ) ) > 0 ) *
					( ( nt1[0] * (i - center[0] ) + 
		                nt1[1] * (j - center[1] ) + 
					    nt1[2] * (k - center[2] ) ) *
			          ( nt2[0] * (i - center[0] ) + 
		                nt2[1] * (j - center[1] ) + 
					    nt2[2] * (k - center[2] ) ) > 0 );

	if ( out == NULL ){
		out = vtkImageData::New();
	}

	transform->Delete();

	nbfVTKInterface converter;
	converter.blitzToVtk( this->myWedge, out );
}

#endif /* FILE_nbf3Dalignment */
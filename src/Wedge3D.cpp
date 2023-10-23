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
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageMathematics.h>
#include <vtkImageFFT.h>
#include <vtkImageRFFT.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageCast.h>
#include <vtkImageShrink3D.h>
#include <vtkImageEllipsoidSource.h>
#include <vtkImageNoiseSource.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <nbfVeselnessFilter.h>
#include <bs/nbfBordStrategyMirror.h>

#include <nbf3Dalignment.h>

void main( int argv, char ** argc )
{
	vtkStructuredPointsReader * spreader = vtkStructuredPointsReader::New();
	spreader->SetFileName( argc[1] );
	spreader->Update();

	vtkImageData * myImage = vtkImageData::New();
	myImage->DeepCopy( spreader->GetOutput() );

	nbfVTKInterface converter;

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName("wedge.vtk");

	//Array< float, 2 > image;
	//nbfMatlabReader mreader;
	//mreader.setFileName("C:/home/project/EM/matlab/image1.blitz");
	//mreader.read(image);

	//int depth = image.rows();
	//Array< float, 3 > myImage1Blitz( image.rows(), image.cols(), depth );
	//for ( int ii = 0; ii < depth; ii++ ){
	//	myImage1Blitz( Range::all(), Range::all(), ii ) = image;
	//}

	//vtkImageData * myImage1 = vtkImageData::New();
	//converter.blitzToVtk( myImage1Blitz, myImage1 );

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToDouble();
	cast->SetInput( myImage );
	cast->Update();

	// read 3D image data
	myImage->DeepCopy( cast->GetOutput() );

	vtkImageNoiseSource * noise = vtkImageNoiseSource::New();
	noise->SetMinimum(0.0);
	noise->SetMaximum(0.0);
	int extent[6];
	myImage->GetDimensions(extent);
	noise->SetWholeExtent( 0, extent[0], 0, extent[1], 0, extent[2] );
	vtkImageMathematics * addNoise = vtkImageMathematics::New();
	addNoise->SetInput1( noise->GetOutput() );
	addNoise->SetInput2( myImage );
	addNoise->SetOperationToAdd();
	addNoise->Update();
	myImage->DeepCopy( addNoise->GetOutput() );

	//writer->SetInput( myImage1 );
	//writer->Write();

	//// second image
	//mreader.setFileName("C:/home/project/EM/matlab/image2.blitz");
	//mreader.read(image);
	//Array< float, 3 > myImage2Blitz( image.rows(), image.cols(), depth );
	//for ( int ii = 0; ii < depth; ii++ ){
	//	myImage2Blitz( Range::all(), Range::all(), ii ) = image;
	//}

	//vtkImageData * myImage2 = vtkImageData::New();
	//converter.blitzToVtk( myImage2Blitz, myImage2 );

	//cast->SetInput( myImage2 );
	//cast->Update();
	//myImage2->DeepCopy( cast->GetOutput() );

	// specify cropped region size
	int size[3];

	//myImage1->GetDimensions(size);
	myImage->GetDimensions(size);

	// build wedge filter in reciprocal space
	Array< float, 3 > myWedge( 1.0*size[0], 1.0*size[1], 1.0*size[2] );
	
	int wedge = 40;

	firstIndex i; secondIndex j; thirdIndex k;

	// new wedge
	float n1[3], n2[3];
	n1[0] = cos( ( 90.0 - wedge ) * vtkMath::DegreesToRadians() );
	n1[1] = 0;
	n1[2] = -sin( ( 90.0 - wedge ) * vtkMath::DegreesToRadians() );
	n2[0] =   n1[0];
	n2[1] =   n1[1];
	n2[2] = - n1[2];
	
	vtkTransform * transform0 = vtkTransform::New();
	transform0->RotateX(0);
	transform0->RotateY(90);
	transform0->RotateZ(0);
	float n1t[3];
	float n2t[3];
	transform0->TransformVector(n1,n1t);
	transform0->TransformVector(n2,n2t);

	myWedge = ( n1[0] * (i - myWedge.ubound(0) / 2.0) + n1[1] * (j - myWedge.ubound(1) / 2.0) + n1[2] * (k - myWedge.ubound(2) / 2.0) ) 
			  * 
			  ( n2[0] * (i - myWedge.ubound(0) / 2.0) + n2[1] * (j - myWedge.ubound(1) / 2.0) + n2[2] * (k - myWedge.ubound(2) / 2.0) ) > 0;
	//myWedge = ( ( n1[0] * (i - myWedge.ubound(0) / 2.0) + n1[1] * (j - myWedge.ubound(1) / 2.0) + n1[2] * (k - myWedge.ubound(2) / 2.0) ) 
	//		  * 
	//		  ( n2[0] * (i - myWedge.ubound(0) / 2.0) + n2[1] * (j - myWedge.ubound(1) / 2.0) + n2[2] * (k - myWedge.ubound(2) / 2.0) ) < 0 ) *
	//		  ( ( n1t[0] * (i - myWedge.ubound(0) / 2.0) + n1t[1] * (j - myWedge.ubound(1) / 2.0) + n1t[2] * (k - myWedge.ubound(2) / 2.0) ) 
	//		  * 
	//		  ( n2t[0] * (i - myWedge.ubound(0) / 2.0) + n2t[1] * (j - myWedge.ubound(1) / 2.0) + n2t[2] * (k - myWedge.ubound(2) / 2.0) )
	//		  < 0 );

	//myWedge = atan( ( k - myWedge.ubound(2) / 2.0 ) / ( i - myWedge.ubound(0) / 2.0 ) );
	//myWedge = where( fabs( myWedge ) < wedge * vtkMath::DegreesToRadians(), 1, 0 );
	//Array< float, 3 > tmp( myWedge.shape() );
	//tmp = myWedge.transpose(firstDim,thirdDim,secondDim);
	//myWedge = tmp.transpose(secondDim,firstDim,thirdDim);

	//myWedge = 1;

	// convert to VTK format
	vtkImageData * myWedgeVtk = vtkImageData::New();
	converter.blitzToVtk( myWedge, myWedgeVtk );

	//writer->SetInput( myWedgeVtk );
	//writer->Write();
	//return;

	// build wedge filter in complex format
	vtkImageData * filter = vtkImageData::New();
	filter->CopyStructure( myWedgeVtk );
	filter->SetScalarTypeToDouble();
	filter->SetNumberOfScalarComponents(2);
	filter->AllocateScalars();
	filter->GetPointData()->GetScalars()->CopyComponent(0, myWedgeVtk->GetPointData()->GetScalars(), 0 );
	filter->GetPointData()->GetScalars()->CopyComponent(1, myWedgeVtk->GetPointData()->GetScalars(), 0 );

	vtkImageFourierCenter * centerfft = vtkImageFourierCenter::New();
	centerfft->SetInput( filter );
	centerfft->Update();

	vtkImageFFT * fft1 = vtkImageFFT::New();
	fft1->SetDimensionality(3);
	fft1->SetInput( myImage );
	fft1->Update();

	// apply wedge filter in reciprocal space
	vtkImageMathematics * math = vtkImageMathematics::New();
	math->SetOperationToComplexMultiply();
	math->SetInput1( fft1->GetOutput() );
	math->SetInput2( centerfft->GetOutput() );
	math->Update();

	vtkImageRFFT * ifft = vtkImageRFFT::New();
	ifft->SetDimensionality(3);
	ifft->SetInput( math->GetOutput() );
	ifft->Update();
	
	// ifftshift( ifft( imf1 .* conj(imf2) ) )
	// centerfft->SetInput( ifft->GetOutput() );
	
	// real( ifftshift( ifft( imf1 .* conj(imf2) ) ) )
	vtkImageExtractComponents * real = vtkImageExtractComponents::New();
	real->SetInput( ifft->GetOutput() );
	real->SetComponents(0);
	real->Update();

	vtkImageData * myImage1 = vtkImageData::New();
	myImage1->DeepCopy( real->GetOutput() );

	// bypass wedge filter
	// myImage1->DeepCopy( myImage );

	double range[2];
	myImage1->GetPointData()->GetScalars()->GetRange( range );

	// scale
	vtkImageMathematics * prod = vtkImageMathematics::New();
	prod->SetOperationToMultiplyByK();
	prod->SetConstantK(10.0 / range[1]);
	prod->SetInput1( myImage1 );
	prod->Update();
	myImage1->DeepCopy( prod->GetOutput() );

	myImage1->GetPointData()->GetScalars()->GetRange( range );
	cout << range[0] << ", " << range[1] << endl;

	//writer->SetInput( myImage1 );
	//writer->Write();

	vtkImageChangeInformation * change = vtkImageChangeInformation::New();
	change->CenterImageOn();
	change->SetInput( myImage );
	change->Update();

	vtkImageReslice * reslice1 = vtkImageReslice::New();
	reslice1->SetInput( change->GetOutput() );
	reslice1->WrapOff();
	vtkTransform * transform1 = vtkTransform::New();
	transform1->RotateX(5);
	transform1->RotateY(0);
	transform1->RotateZ(0);
	transform1->Translate(0,0,0);
	reslice1->SetResliceTransform( transform1 );
	reslice1->Update();

	vtkImageData * myImage2 = vtkImageData::New();
	myImage2->DeepCopy( reslice1->GetOutput() );

	// rotate wedge to simulate volume cropping
	change->SetInput( myWedgeVtk );
	change->Update();

	vtkImageReslice * reslice2 = vtkImageReslice::New();
	reslice2->SetInput( change->GetOutput() );
	reslice2->WrapOn();
	reslice2->SetBackgroundLevel(1);
	//reslice2->SetInterpolationModeToNearestNeighbor();
	vtkTransform * transform2 = vtkTransform::New();
	transform2->RotateX(0);
	transform2->RotateY(0);
	transform2->RotateZ(0);
	//transform2->RotateX(-30);
	//transform2->RotateY(49);
	//transform2->RotateZ(33);
	reslice2->SetResliceTransform( transform2 );
	reslice2->Update();

	//writer->SetInput( reslice2->GetOutput() );
	//writer->Write();
	//return;

	// BYPASS - USE SAME WEDGE
	//filter->CopyStructure( reslice2->GetOutput() );
	//filter->SetScalarTypeToDouble();
	//filter->SetNumberOfScalarComponents(2);
	//filter->AllocateScalars();
	//filter->GetPointData()->GetScalars()->CopyComponent(0, reslice2->GetOutput()->GetPointData()->GetScalars(), 0 );
	//filter->GetPointData()->GetScalars()->CopyComponent(1, reslice2->GetOutput()->GetPointData()->GetScalars(), 0 );

	centerfft->Update();

	fft1->SetInput( myImage2 );
	math->Update();
	real->Update();

	myImage2->DeepCopy( real->GetOutput() );
	myImage2->GetPointData()->GetScalars()->GetRange( range );

	prod->SetInput1( myImage2 );
	prod->SetConstantK( 10.0 / range[1] );
	prod->Update();
	myImage2->DeepCopy( prod->GetOutput() );

	myImage2->GetPointData()->GetScalars()->GetRange( range );
	cout << range[0] << ", " << range[1] << endl;

	noise->SetMaximum(0.0);
	addNoise->SetInput1( noise->GetOutput() );
	addNoise->SetInput2( myImage2 );
	addNoise->SetOperationToAdd();
	addNoise->Update();
	myImage2->DeepCopy( addNoise->GetOutput() );

	writer->SetInput( myImage1 );
	writer->Write();

	writer->SetInput( myImage2 );
	writer->Write();

#if 1
	nbf3DAlignment< float > align;
	align.setImage1( myImage1 );
	align.setImage2( myImage2 );
	//align.setWedgeImage( myWedge );
	align.setWedge(wedge);
	align.execute( myImage1 );
#else

	return;

	// pre-center for rotation
	vtkImageChangeInformation * change = vtkImageChangeInformation::New();
	change->SetInput( myWedgeVtk );
	change->CenterImageOn();

	// rotate wedge
	vtkTransform * transform = vtkTransform::New();
	transform->RotateZ(-28.0);
	// transform->RotateY(23.0);
	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->SetResliceTransform( transform );
	reslice->WrapOn();
	reslice->SetInput( change->GetOutput() );

	// compute composition of wedges (given+rotated)
	vtkImageMathematics * math = vtkImageMathematics::New();
	math->SetOperationToMultiply();
	math->SetInput1( reslice->GetOutput() );
	math->SetInput2( change->GetOutput() );
	math->Update();

	myWedgeVtk->DeepCopy( math->GetOutput() );

	// build wedge filter in complex format
	vtkImageData * complexFilter = vtkImageData::New();
	complexFilter->CopyStructure( myWedgeVtk );
	complexFilter->SetScalarTypeToDouble();
	complexFilter->SetNumberOfScalarComponents(2);
	complexFilter->AllocateScalars();
	complexFilter->GetPointData()->GetScalars()->CopyComponent(0, myWedgeVtk->GetPointData()->GetScalars(), 0 );
	complexFilter->GetPointData()->GetScalars()->CopyComponent(1, myWedgeVtk->GetPointData()->GetScalars(), 0 );

	// shift FFT before applying filter
	vtkImageFourierCenter * centerfft = vtkImageFourierCenter::New();
	centerfft->SetInput( complexFilter );
	centerfft->Update();
	complexFilter->DeepCopy( centerfft->GetOutput() );

	// FFT of input image
	vtkImageFFT * fft = vtkImageFFT::New();
	fft->SetDimensionality(3);
	fft->SetInput( myImage1 );
	fft->Update();

	// multiply input and filter in frquency space
	math->SetOperationToComplexMultiply();
	math->SetInput1( complexFilter );
	math->SetInput2( fft->GetOutput() );
	math->Update();

	// kill DC value
	vtkImageData * myImage1f = vtkImageData::New();
	myImage1f->DeepCopy( math->GetOutput() );
	myImage1f->SetScalarComponentFromDouble(0,0,0,0,0);
	myImage1f->SetScalarComponentFromDouble(0,0,0,1,0);

	// math pipeline
	vtkImageMathematics * conj = vtkImageMathematics::New();	
	conj->SetOperationToConjugate();
	conj->SetInput1( myImage1f );
	conj->Update();

	vtkImageMathematics * mult = vtkImageMathematics::New();	
	mult->SetOperationToComplexMultiply();
	mult->SetInput1( myImage1f );
	mult->SetInput2( conj->GetOutput() );
	mult->Update();

	vtkImageMathematics * sqrt = vtkImageMathematics::New();	
	sqrt->SetOperationToSquareRoot();
	sqrt->SetInput1( mult->GetOutput() );
	sqrt->Update();

	vtkImageShrink3D * shrink = vtkImageShrink3D::New();
	shrink->SetShrinkFactors( size );
	shrink->AveragingOn();
	shrink->SetInput( sqrt->GetOutput() );
	shrink->Update();

	vtkImageMathematics * prem = vtkImageMathematics::New();	
	prem->SetOperationToMultiplyByK();
	double constantK = 1.0 / shrink->GetOutput()->GetScalarComponentAsDouble(0,0,0,0);
	prem->SetConstantK( constantK );
	prem->SetInput1( myImage1f );
	prem->Update();

	// store result in original image
	myImage1f->DeepCopy( prem->GetOutput() );

	// now get second correlation component

	change->SetInput( myImage2 );
	change->Update();
	reslice->SetInput( change->GetOutput() );
	reslice->Update();
	myImage2->DeepCopy( reslice->GetOutput() );

	writer->SetInput( myImage2 );
	writer->Write();

	// FFT of input image
	vtkImageFFT * fft2 = vtkImageFFT::New();
	fft2->SetDimensionality(3);
	fft2->SetInput( myImage2 );
	fft2->Update();

	// multiply input and filter in frquency space
	vtkImageMathematics * math2 = vtkImageMathematics::New();
	math2->SetOperationToComplexMultiply();
	math2->SetInput1( complexFilter );
	math2->SetInput2( fft2->GetOutput() );
	math2->Update();

	// kill DC value
	vtkImageData * myImage2f = vtkImageData::New();
	myImage2f->DeepCopy( math2->GetOutput() );
	myImage2f->SetScalarComponentFromDouble(0,0,0,0,0);
	myImage2f->SetScalarComponentFromDouble(0,0,0,1,0);

	// math pipeline
	vtkImageMathematics * conj2 = vtkImageMathematics::New();	
	conj2->SetOperationToConjugate();
	conj2->SetInput1( myImage2f );
	conj2->Update();

	vtkImageMathematics * mult2 = vtkImageMathematics::New();	
	mult2->SetOperationToComplexMultiply();
	mult2->SetInput1( myImage2f );
	mult2->SetInput2( conj2->GetOutput() );
	mult2->Update();

	vtkImageMathematics * sqrt2 = vtkImageMathematics::New();	
	sqrt2->SetOperationToSquareRoot();
	sqrt2->SetInput1( mult2->GetOutput() );
	sqrt2->Update();

	vtkImageShrink3D * shrink2 = vtkImageShrink3D::New();
	shrink2->SetShrinkFactors( size );
	//shrink2->AveragingOn();
	shrink2->SetInput( sqrt2->GetOutput() );
	shrink2->Update();

	vtkImageMathematics * prem2 = vtkImageMathematics::New();	
	prem2->SetOperationToMultiplyByK();
	double constantK2 = 1.0 / shrink2->GetOutput()->GetScalarComponentAsDouble(0,0,0,0);
	prem2->SetConstantK( constantK2 );
	prem2->SetInput1( myImage2f );
	prem2->Update();

	// store result in original image
	myImage2f->DeepCopy( prem2->GetOutput() );

	vtkImageMathematics * finalConj = vtkImageMathematics::New();
	finalConj->SetOperationToConjugate();
	finalConj->SetInput1( myImage2f );
	finalConj->Update();
	
	vtkImageMathematics * finalProd = vtkImageMathematics::New();
	finalProd->SetOperationToComplexMultiply();
	finalProd->SetInput1( myImage1f );
	finalProd->SetInput2( finalConj->GetOutput() );
	finalProd->Update();

	// bring back to image domain
	vtkImageRFFT * ifft = vtkImageRFFT::New();
	ifft->SetDimensionality(3);
	ifft->SetInput( finalProd->GetOutput() );
	ifft->Update();

	centerfft->SetInput( ifft->GetOutput() );
	centerfft->Update();

	// extract real part
	vtkImageExtractComponents * real = vtkImageExtractComponents::New();
	real->SetInput( centerfft->GetOutput() );
	real->SetComponents(0);
	real->Update();

	real->GetOutput()->GetPointData()->GetScalars()->GetRange( range );
	//cout << range[0] << endl;
	cout << "score = " << range[1] << endl;

	// write result
	writer->SetInput( real->GetOutput() );
	writer->Write();

	writer->Delete();
	myImage1->Delete();
	//reader->Delete();
#endif
}
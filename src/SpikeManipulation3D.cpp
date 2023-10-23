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
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>

#include <io/nbfVTKInterface.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfMatlabReader.h>
#include <bs/nbfBordStrategyMirror.h>

#include <em/nbfImageMetric.h>
#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfProjectionRotationMetric3D.h>
#include <em/nbfWedgedSubImage3D.h>

#define PIXEL double

void main( int argc, char ** argv )
{
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( argv[1] );
	reader->Update();

	Array< PIXEL,3 > Dummy;
	nbfVTKInterface::vtkToBlitzReference( reader->GetOutput(), Dummy );

	vtkTransform * transform = vtkTransform::New();
	transform->RotateZ(-30);
	transform->RotateX(90);
	transform->Translate(4,2,0);

	vtkImageChangeInformation * change = vtkImageChangeInformation::New();
	change->SetInput( reader->GetOutput() );
	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->SetInterpolationModeToCubic();
	reslice->SetResliceTransform( transform );
	reslice->SetInput( change->GetOutput() );
	reslice->SetBackgroundLevel( max(Dummy) );
	reslice->Update();

	vtkImageData * rotated = vtkImageData::New();
	rotated->DeepCopy( reslice->GetOutput() );

	nbfVTKInterface::vtkToBlitzReference( rotated, Dummy );

	vtkCylindricalTransform * cylindrical = vtkCylindricalTransform::New();

	int dims[3];
	reslice->GetOutput()->GetDimensions(dims);

	//change->SetInput( rotated );
	//change->Update();
	//reslice->SetResliceTransform( cylindrical );
	//reslice->SetOutputExtent( 0, ( dims[0] - 1 ) / 2.0, 0, 179, 0, ( dims[2] - 1 ) / 2.0 );
	//reslice->SetOutputOrigin(0,0,0);
	//reslice->SetOutputSpacing( 1, 2.0 * vtkMath::Pi() / 179.0, 1 );	
	//reslice->Update();

	Array< PIXEL, 3 > C;
	nbfVTKInterface::vtkToBlitz( reslice->GetOutput(), C );
	
	firstIndex i;
	secondIndex j;
	thirdIndex k;

	Array< PIXEL, 2 > Pp1( C.rows(), C.depth() );
	Pp1 = sum( C(i,k,j),k);

	nbfMatlabWriter w;
	w.setFileName("p.matlab");
	w.write(Pp1);

	Array< PIXEL,2 > Pp( Pp1.rows()/2, Pp1.cols()  );
	Pp = Pp1( Range( 0, Pp1.rows()/2-1 ), Range::all() );
	Array< PIXEL, 2 > P( Pp.transpose(secondDim,firstDim) );
	P.reverseSelf(firstDim);
	w.write(P);

	// interpolate maxima value and assign to attribute
	nbfLinearInterpolator< PIXEL, 2 > interp( P );
	
	TinyVector< PIXEL, 3 > center( dims[0]/2, dims[1]/2 );

	Array< PIXEL,3 > model( dims[0], dims[1], dims[2]/2 );
	Array< PIXEL, 3 > :: iterator iter = model.begin();
	while ( iter != model.end() ){
		TinyVector< int, 3 > p = iter.position();
		PIXEL r = sqrt( 0.0 + pow2(p[0]-center[0]) + pow2(p[1]-center[1]) );
		(*iter) = interp.interpolateSingle( TinyVector< PIXEL, 2 >(r,p[2]) );
		++iter;
	}

	w.setFileName("p.matlab");
	w.write(Dummy);

	vtkImageData * final = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( model, final );

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName( argv[2] );
	writer->SetInput( final );
	writer->Write();
	writer->Delete();
}

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageCast.h>
#include <vtkImageReader.h>
#include <vtkPiecewiseFunction.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkImageChangeInformation.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkTransformFilter.h>
#include <vtkColorTransferFunction.h>
#include <vtkImageMathematics.h>
#include <vtkImageContinuousDilate3D.h>
#include <vtkImageMathematics.h>
#include <vtkPNGWriter.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkVolume.h>
#include <vtkVolumeMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkProbeFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>


void main( int argv, char ** argc )
{
	vtkImageReader * mrcReader = vtkImageReader::New();
	mrcReader->SetFileName( argc[1] );
	mrcReader->SetDataByteOrderToLittleEndian();
	mrcReader->SetFileDimensionality(3);

	// read header manually
	ifstream inFile( argc[1], ios::in | ios::binary );
	if ( inFile.is_open() != 1 ){
		inFile.close();
		return FALSE;		
	}
	long type;
	inFile.read((char *)&Nx, sizeof (Nx));
	inFile.read((char *)&Ny, sizeof (Ny));
	inFile.read((char *)&Nz, sizeof (Nz));
	inFile.read((char *)&type, sizeof (type));

	long fileSize;
	inFile.seekg(0, ios::beg);
	fileSize = inFile.tellg();
	inFile.seekg(0, ios::end);
	fileSize = inFile.tellg() - fileSize;
	inFile.close();

	if ( type == 1 ){				
		mrcReader->SetDataScalarTypeToShort();
		mrcReader->SetHeaderSize( fileSize - Nx * Ny * Nz * 2 );
	}
	mrcReader->SetDataExtent(0,Nx-1,0,Ny-1,0,Nz-1);
	mrcReader->Update();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToFloat();
	cast->SetInput( mrcReader->GetOutput() );

	vtkImageChangeInformation * change = vtkImageChangeInformation::New();
	change->CenterImageOn();
	change->SetINput( cast->GetOutput() );

	vtkTransform * rotate = vtkTransform :: New();
	rotate->RotateX(-30);

	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->SetInput( change->GetOutput() );
	reslice->SetResliceTransform( rotate );

	vtkImageShrink3D * shrink = vtkImageShrink3D::New();
	shrink->SetInput( reslice->GetOutput() );
	shrink->SetShrinkFactors(Nx,Ny,1);
	shrink->AveragingOn();
}
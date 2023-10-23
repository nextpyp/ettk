
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
	vtkStructuredPointsReader * reader = vtkStructuredPointsReader::New();
	reader->SetFileName( argc[1] );
	reader->Update();

	vtkImageCast * cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToUnsignedShort();
	cast->SetInput( reader->GetOutput() );

	vtkRenderer * renderer = vtkRenderer::New();
	renderer->GetActiveCamera()->ParallelProjectionOn();

	vtkRenderWindow * window = vtkRenderWindow::New();

	window->AddRenderer( renderer );

	vtkPiecewiseFunction * opacityTransferFunction = vtkPiecewiseFunction::New();
	opacityTransferFunction->AddPoint(20, 0.0);
	opacityTransferFunction->AddPoint(255, 0.2);

	vtkColorTransferFunction * colorTransferFunction = vtkColorTransferFunction::New();
	colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	colorTransferFunction->AddRGBPoint(64.0, 1.0, 0.0, 0.0);
	colorTransferFunction->AddRGBPoint(128.0, 0.0, 0.0, 1.0);
	colorTransferFunction->AddRGBPoint(192.0, 0.0, 1.0, 0.0);
	colorTransferFunction->AddRGBPoint(255.0, 0.0, 0.2, 0.0);

	vtkVolumeProperty * prop = vtkVolumeProperty::New();
	prop->SetColor(colorTransferFunction);
	prop->SetScalarOpacity(0.5,opacityTransferFunction);
	prop->SetInterpolationTypeToLinear();
	prop->ShadeOn();

	vtkVolumeRayCastMapper * mapper = vtkVolumeRayCastMapper::New();
	mapper->SetInput( cast->GetOutput() );
	vtkVolumeRayCastCompositeFunction * compositeFunction = vtkVolumeRayCastCompositeFunction::New();
	mapper->SetVolumeRayCastFunction(compositeFunction);

	vtkVolume * volume = vtkVolume::New();
	volume->SetMapper(mapper);
	volume->SetProperty(prop);

	renderer->AddVolume( volume );
	renderer->SetBackground(1, 1, 1);
	int size[3];
	reader->GetOutput()->GetDimensions(size);
	window->SetSize( size[0], size[1]);
	window->Render();

	vtkRenderWindowInteractor * iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow( window );
	iren->Initialize();

	vtkWindowToImageFilter * toImage = vtkWindowToImageFilter::New();
	toImage->SetInput( window );

	vtkPNGWriter * writer = vtkPNGWriter::New();
	writer->SetInput( toImage->GetOutput() );
	writer->SetFileName( argc[2] );
	writer->Write();
	writer->Delete();

	iren->Start();
}
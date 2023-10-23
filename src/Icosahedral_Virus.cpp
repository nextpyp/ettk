#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

#include <nbfTimer.h>

#include <io/nbf3DReader.h>
#include <bs/nbfBordStrategyMirror.h>
#include <io/nbfMatlabReader.h>
#include <io/nbfMatlabWriter.h>
#include <io/nbfVTKInterface.h>

#include <vtkPolyDataReader.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkGlyph3D.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkSuperquadricSource.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkImplicitModeller.h>
#include <vtkStructuredPointsWriter.h>

void main( int argc, char ** argv )
{
	// get spike geometry
	vtkPolyDataReader * preader = vtkPolyDataReader :: New();
	preader->SetFileName( argv[1] );
	preader->Update();

	float phi = ( 1.0 + sqrt( 5.0 ) ) / 2.0;

	// isosahedral vertices
	vtkPoints * points = vtkPoints :: New();
	points->InsertNextPoint( 0,  1,  phi );
	points->InsertNextPoint( 0,  1, -phi );
	points->InsertNextPoint( 0, -1,  phi );
	points->InsertNextPoint( 0, -1, -phi );
	
	points->InsertNextPoint(  1,  phi, 0 );
	points->InsertNextPoint(  1, -phi, 0 );
	points->InsertNextPoint( -1,  phi, 0 );
	points->InsertNextPoint( -1, -phi, 0 );
	
	points->InsertNextPoint(  phi, 0,  1 );
	points->InsertNextPoint(  phi, 0, -1 );
	points->InsertNextPoint( -phi, 0,  1 );
	points->InsertNextPoint( -phi, 0, -1 );
	
	vtkFloatArray * normals = vtkFloatArray :: New();
	normals->SetNumberOfComponents(3);
	for ( int i = 0; i < points->GetNumberOfPoints(); i++ ){
		double p[3];
		points->GetPoint(i,p);
		normals->InsertNextTuple3( p[0], p[1], p[2] );
	}

	vtkPolyData * pdata = vtkPolyData :: New();
	pdata->SetPoints( points );
	pdata->GetPointData()->SetNormals( normals );

	normals->Delete();

	vtkSuperquadricSource * ellipse = vtkSuperquadricSource::New();
	ellipse->ToroidalOff();
	ellipse->SetThickness(10);
	ellipse->SetScale(3,1,1);

	vtkGlyph3D * glyph = vtkGlyph3D :: New();
	glyph->SetSource( ellipse->GetOutput() );
	glyph->SetInput( pdata );
	glyph->SetScaleFactor(.35);
	glyph->SetVectorModeToUseNormal();

	vtkSphereSource * sphere = vtkSphereSource :: New();
	sphere->SetRadius(1.75);
	sphere->SetThetaResolution(10);
	sphere->SetPhiResolution(10);

	vtkAppendPolyData * virus = vtkAppendPolyData :: New();
	virus->AddInput( sphere->GetOutput() );
	virus->AddInput( glyph->GetOutput() );

	vtkImplicitModeller * implicit = vtkImplicitModeller :: New();
	implicit->SetInput( virus->GetOutput() );
	implicit->SetMaximumDistance(500);
	implicit->SetSampleDimensions(64,64,64);
	implicit->AdjustBoundsOn();
	// implicit->SetCapValue(5);
	implicit->SetAdjustDistance(1.2);
	implicit->Update();

	Array< float, 3 > A;
	nbfVTKInterface :: vtkToBlitz( implicit->GetOutput(), A );
	cout << "A = [ " << min(A) << "," << max(A) << "]" << endl;
	A = where( abs(A) > 10, 10, A );
	nbfVTKInterface :: blitzToVtk( A, implicit->GetOutput() );

	vtkStructuredPointsWriter * swriter = vtkStructuredPointsWriter :: New();
	swriter->SetInput( implicit->GetOutput() );
	swriter->SetFileName( argv[3] );
	swriter->Write();

	vtkXMLDataSetWriter * pwriter = vtkXMLDataSetWriter :: New();
	pwriter->SetInput( virus->GetOutput() );
	pwriter->SetFileName( argv[2] );
	pwriter->Write();
	pwriter->Delete();

	pdata->Delete();
	points->Delete();
	glyph->Delete();
	preader->Delete();
}
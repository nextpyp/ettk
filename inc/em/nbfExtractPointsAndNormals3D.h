/** @file nbfCutSubVolumes.h
*	Template search inside a volume constrained to a given surface.
*	Given an image, a template and a surface, this class will produce the
*   locations _on_ the surface where the cross-correlation locally maximizes.
*/

#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkContourFilter.h>
#include <vtkDecimatePro.h>


#include <io/nbfVTKInterface.h>

using namespace blitz;

/** Extract triangulated surface from implicit surface representation.
	Polygonal data is smoothed after contouring.
    Normals are also computed for each point on the surface.
	
	Output produced by this filter can be used to perform template search to locate
	positions of interest where CCC wrt a template is maximized.

	@see nbfTemplateSearchEM
	
*/

template< class Pixel >
class nbfExtractPointsAndNormals3D
{
public:

	nbfExtractPointsAndNormals3D();

	~nbfExtractPointsAndNormals3D();

	/// Produce triangulated surface representation with normals.
	void execute( vtkPolyData * );

	/// Set input implicit surface and level set to extract.
	void setSurface( vtkImageData *, Pixel );
	
	/// The implicit surface is usualy computed on a binned image.
	/// To produce correcly scaled 3D positions, the binning factor must be applied back.
	/// If the original image was binned by 2, a magnification factor of 2 has to be applied.
	void setMagnification( Pixel m ){ this->magnification = m; }

protected:

	/// Implicit surface image
	vtkImageData * surfaceImageVtk;

	Pixel magnification, level;
};

template< class Pixel >
nbfExtractPointsAndNormals3D< Pixel > :: nbfExtractPointsAndNormals3D()
{
	this->surfaceImageVtk = NULL;
	this->magnification = 1;
	this->level = 0;
}


template< class Pixel >
nbfExtractPointsAndNormals3D< Pixel > :: ~nbfExtractPointsAndNormals3D()
{
}


template< class Pixel >
void nbfExtractPointsAndNormals3D< Pixel > :: setSurface( vtkImageData * in, Pixel level )
{
	this->surfaceImageVtk = in;
	this->level = level;
}

template< class Pixel >
void nbfExtractPointsAndNormals3D< Pixel > :: execute( vtkPolyData * points )
{	
	// get contour of implicit surface
	vtkContourFilter * contour = vtkContourFilter::New();
	contour->SetInput( this->surfaceImageVtk );
	contour->SetValue(0, this->level );
	contour->Update();

	//vtkDecimatePro * simplify = vtkDecimatePro :: New();
	//simplify->SetInput( contour->GetOutput() );
	//simplify->SetTargetReduction(.9);
	//simplify->Update();

	// smooth polygonal data
	vtkSmoothPolyDataFilter * smooth = vtkSmoothPolyDataFilter::New();
	smooth->SetInput( contour->GetOutput() );
	smooth->SetNumberOfIterations(200);
	smooth->FeatureEdgeSmoothingOff();
	smooth->Update();

	//vtkXMLPPolyDataWriter * writer = vtkXMLPPolyDataWriter::New();
	//writer->SetInput( smooth->GetOutput() );
	//writer->SetFileName( "before.xml" );
	//writer->Write();
	//writer->SetInput( simplify->GetOutput() );
	//writer->SetFileName( "after.xml" );
	//writer->Write();

	vtkPolyDataNormals * cnormals = vtkPolyDataNormals::New();
	cnormals->SetInput( smooth->GetOutput() );
	cnormals->SplittingOff();
	cnormals->Update();

	vtkTransform * transform = vtkTransform::New();
	transform->Scale( this->magnification, this->magnification, this->magnification );
	vtkTransformPolyDataFilter * rescale = vtkTransformPolyDataFilter::New();
	rescale->SetInput( cnormals->GetOutput() );
	rescale->SetTransform( transform );
	rescale->Update();

	points->DeepCopy( rescale->GetOutput() );

	double normal[3];
	for ( int i = 0; i < points->GetNumberOfPoints(); i++ ){

		points->GetPointData()->GetNormals()->GetTuple(i,normal);

		// project normal into horizontal direction
		// normal[2] = 0;

		points->GetPointData()->GetNormals()->SetTuple(i,normal);
	}

	rescale->Delete();
	transform->Delete();
	cnormals->Delete();
	//simplify->Delete();
	smooth->Delete();
	contour->Delete();

	//int dims[3];
	//this->surfaceImageVtk->GetDimensions(dims);

	//for ( int i = 0; i < contour->GetOutput()->GetNumberOfPoints(); i++ ){
	//	Pixel position[3];
	//	Pixel normal[3];
	//	contour->GetOutput()->GetPoint(i,position);

	//	cnormals->GetOutput()->GetPointData()->GetNormals()->GetTuple(i,normal);

	//	normal[1] = - normal[1];
	//	normal[2] = 0;

	//	// compute alignment angles with Z and Y axes
	//	double theta = - atan2( normal[1], -normal[0] ) * vtkMath::RadiansToDegrees();
	//	double phi = atan2( sqrt( pow2(normal[0]) + pow2(normal[1]) ), -normal[2] ) * vtkMath::RadiansToDegrees();


	//	// add point to list of candidates (rescale back to full resolution)
	//	rootPoints.push_back( TinyVector<Pixel,3>( position[0], dims[1] - position[1], position[2] ) / this->magnification );
	//	angles.push_back( TinyVector< Pixel, 3 >( 0.0, phi, theta ) );
	//	normals.push_back( TinyVector< Pixel, 3 >( normal[0], normal[1], normal[2] ) );

	//}
}
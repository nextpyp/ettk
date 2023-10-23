/** @file nbfCutSubVolumes.h
*	Template search inside a volume constrained to a given surface.
*	Given an image, a template and a surface, this class will produce the
*   locations _on_ the surface where the cross-correlation locally maximizes.
*/

#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageReslice.h>
#include <vtkImageChangeInformation.h>
#include <vtkTransform.h>
#include <vtkImageDilateErode3D.h>

#include <io/nbfVTKInterface.h>

using namespace blitz;

/** Extract sub-volumes at given positions and orientations.
    
	Detailed description. Particle extraction.
	
	@see nbfFastMarching
	
	@todo Provide methods for handling correctly.	
*/

template< class Pixel >
class nbfCutSubVolumes
{
public:

	/// constructor takes weight array as input
	nbfCutSubVolumes();

	~nbfCutSubVolumes();

	// cut out volumes
	// inputs:
	//		1. input image
	//		2. implicit membrane representation
	// output:
	//		3. vector containig individual volumes
	void execute( vector< Array< short, 3 > > &, vector< vtkTransform * > & );

	void setTemplate( Array< short, 3 > & );
	void setTemplate( vtkImageData * );

	void setTemplateFullResolution( Array< short, 3 > & );

	void setSurface( vtkImageData * );
	
	void setImage( vtkImageData * );

	/// Set full resolution image to cut the sub-volumes from.
	void setFullResolutionImage( vtkImageData * );

	/// Sets how many voxels below the surface to include in the cropped images
	void setOffsetUnderSurface( int );

	void setOffsetUnderSurfaceFullResolution( int );

	/// Correlation threshold, only keep peaks above threshold
	void setCorrelationTh( Pixel in ){ this->correlationTh = in; }

	void setMinimumSeparation( int in ){ this->neighborhoodSize = in; }

	void setSurfaceLbound( Pixel lbound ){ this->lbound = lbound; }
	void setSurfaceUbound( Pixel ubound ){ this->ubound = ubound; }

	void setHeightRange( int a, int b, int c ){ this->dimBound = a; this->hlbound = b; this->hubound = c; }

	/// Choose whether to ignore image values below the surface for computations.
	void setToIgnoreValuesBelowSurface(){ this->ignoreValuesBelowSurface = true; }
	void setNotToIgnoreValuesBelowSurface(){ this->ignoreValuesBelowSurface = false; }

protected:

	/// Keep local maxima of correlation function
	void nonLocalMaximumSuppression( vector< TinyVector< int, 3 > > &, vector< Pixel > & );

	/// Extract volume at specified coordinates
	void extractSubVolume( TinyVector< int, 3 > &, Array< short, 3 > &, vtkTransform * = NULL, bool = false );

	// input image data
	vtkImageData * inputImageVtk;

	// input full resolution image data
	vtkImageData * inputFullResolutionImageVtk;

	/// Template reference
	Array< short, 3 > templateImageBlitz;
	vtkImageData * templateImageVtk;
	vtkImageData * templateImageFullResolutionVtk;
	
	/// Implicit surface image
	vtkImageData * surfaceImageVtk;

	/// How many voxels to look below the surface
	int offsetUnderSurface;
	int offsetUnderSurfaceFullResolution;

	/// Neighborhood size for the computation of local maxima of the correlation function
	int neighborhoodSize;

	/// Flag to track destruction of template image (destroy only if created in this class)
	bool destroyTemplateImage;

	/// Correlation threshold, only keep local maxima above threshold value
	Pixel correlationTh;

	/// Set surface range to compute correlation
	Pixel lbound, ubound;

	int dimBound, hlbound, hubound;

	/// Object for performing the cutting of volumes
	vtkImageChangeInformation * change;
	vtkImageChangeInformation * change2;
	//vtkTransform * transform;
	vtkImageReslice * reslice;

	vtkImageChangeInformation * changeSurface;
	vtkImageReslice * resliceSurface;

	bool ignoreValuesBelowSurface;
};

template< class Pixel >
nbfCutSubVolumes< Pixel > :: nbfCutSubVolumes()
{
	this->templateImageVtk = NULL;
	this->templateImageFullResolutionVtk = NULL;
	this->surfaceImageVtk = NULL;
	this->inputImageVtk = NULL;
	this->inputFullResolutionImageVtk = NULL;
	this->offsetUnderSurface = 0;
	this->offsetUnderSurfaceFullResolution = 0;
	this->neighborhoodSize = 20;
	this->correlationTh = 0.5;
	this->destroyTemplateImage = false;

	this->ignoreValuesBelowSurface = true;

	this->lbound = -1;
	this->ubound = 0;

	// initialize auxiliary VTK objects and pipeline
	this->change = vtkImageChangeInformation::New();
	this->change2 = vtkImageChangeInformation::New();
	this->reslice = vtkImageReslice::New();

	this->changeSurface = vtkImageChangeInformation::New();
	this->resliceSurface = vtkImageReslice::New();

}


template< class Pixel >
nbfCutSubVolumes< Pixel > :: ~nbfCutSubVolumes()
{
	if ( this->destroyTemplateImage != NULL ){
		this->templateImageVtk->Delete();
	}

	if ( this->templateImageFullResolutionVtk != NULL ){
		this->templateImageFullResolutionVtk->Delete();
	}

	this->reslice->Delete();
	this->change2->Delete();
	this->change->Delete();

	this->resliceSurface->Delete();
	this->changeSurface->Delete();

}


template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setTemplate( Array< short, 3 > & in )
{
	this->templateImageBlitz.reference( in );
	nbfVTKInterface converter;
	this->templateImageVtk = vtkImageData::New();
	converter.blitzToVtk( in, this->templateImageVtk );
	this->destroyTemplateImage = true;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setTemplate( vtkImageData * in )
{
	this->templateImageVtk = in;
	nbfVTKInterface converter;
	converter.vtkToBlitz( in, this->templateImageBlitz );
	this->destroyTemplateImage = false;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setTemplateFullResolution( Array< short, 3 > & in )
{
	this->templateImageFullResolutionVtk = vtkImageData::New();
	nbfVTKInterface::blitzToVtk( in, this->templateImageFullResolutionVtk );
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setSurface( vtkImageData * in )
{
	this->surfaceImageVtk = in;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setImage( vtkImageData * in )
{
	this->inputImageVtk = in;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setFullResolutionImage( vtkImageData * in )
{
	this->inputFullResolutionImageVtk = in;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setOffsetUnderSurface( int in )
{
	this->offsetUnderSurface = in;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: setOffsetUnderSurfaceFullResolution( int in )
{
	this->offsetUnderSurfaceFullResolution = in;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: execute( vector< Array< short, 3 > > & volumeList, vector< vtkTransform * > & transformList )
{
	// check if all necessary parameters are available
	if ( ( this->templateImageBlitz.size() == 0 ) ||
		 ( this->inputImageVtk == NULL ) ||
		 ( this->surfaceImageVtk == NULL ) )
	{
		// not all inputs specified, cannot continue
		return;
	}

	// initialize VTK pipeline
	this->change->SetInput( this->inputImageVtk );

	this->change2->SetInput( this->templateImageVtk );
	int size[3];
	this->templateImageVtk->GetDimensions( size );

	// set relative position of template wrt the image volume	
	this->change2->SetOriginTranslation( - size[0] / 2.0, - size[1] / 2.0, - size[2] + this->offsetUnderSurface );
	
	this->reslice->SetInput( this->change->GetOutput() );
	this->reslice->SetInformationInput( this->change2->GetOutput() );

	this->changeSurface->SetInput( this->surfaceImageVtk );
	this->resliceSurface->SetInput( this->changeSurface->GetOutput() );
	this->resliceSurface->SetInformationInput( this->change2->GetOutput() );

	// store position of candidate volumes
	vector< TinyVector< int, 3 > > rootPoints;
	
	// store position of candidate volumes
	vector< Pixel > correlation;

	// build blitz input image from vtk
	Array< float, 3 > surfaceImageBlitz;
	nbfVTKInterface converter;
	converter.vtkToBlitz( this->surfaceImageVtk, surfaceImageBlitz );

	// temporary store of cut-out volumes
	Array< short, 3 > S;

	cout << "traversing VOI and computing ccc..." << endl;
	
	vtkImageData * ccc3d = vtkImageData::New();
	ccc3d->CopyStructure( this->inputImageVtk );
	ccc3d->SetScalarTypeToDouble();
	ccc3d->SetNumberOfScalarComponents(1);
	ccc3d->AllocateScalars();

	//int counter = 0;

	typename Array< float, 3 > :: iterator iter = surfaceImageBlitz.begin();
	while ( iter != surfaceImageBlitz.end() ){
		TinyVector< int, 3 > p = iter.position();
		
		// if close to the surface  && ( iter.position()[thirdDim] == 148 )
		if ( ( (*iter) < this->ubound ) && ( (*iter) > this->lbound ) && ( p[this->dimBound] > this->hlbound ) && ( p[this->dimBound] < this->hubound ) ){
		//if ( sum( abs( iter.position() - TinyVector<int,3>(30,95,80) ) ) == 0 ){
		//if ( sum( abs( iter.position() - TinyVector<int,3>(31,97,83) ) ) == 0 ){
			
			TinyVector< int, 3 > p( iter.position() );

			// add point to list of candidates
			rootPoints.push_back( p );

			// cout << p << endl;

			vtkTransform * transform = vtkTransform::New();

			// extract sub volume
			this->extractSubVolume( p, S, transform, true );

			//if ( sum( abs( p - TinyVector<int,3>(166,375,178) ) ) < 5 ){
				//nbfMatlabWriter writer;
				//writer.setFileName("joder.array");
				//writer.write(S);
			//}

			// compute and store correlation
			Pixel md = mean( S );
			Pixel mc = mean( this->templateImageBlitz ); 
			Pixel ccc = sum( ( S - md ) * ( this->templateImageBlitz - mc ) ) / sqrt( sum( pow2( S - md ) ) * sum( pow2( this->templateImageBlitz - mc ) ) );
			//ccc = -md;

			double dccc = ccc;
			ccc3d->SetScalarComponentFromDouble( iter.position()[0], iter.position()[1], iter.position()[2], 0, ccc );

			//counter++;
			//cout << "candidate " << counter << ", ccc = " << ccc << endl;

			//nbfMatlabWriter mwr;
			//mwr.setFileName("ipath");
			//mwr.write(S);
			correlation.push_back( ccc );
		}
		++iter;
	}

	vtkStructuredPointsWriter * owriter = vtkStructuredPointsWriter::New();
	owriter->SetFileName("ccc3d.vtk");
	owriter->SetInput( ccc3d );
	owriter->SetFileTypeToBinary();
	owriter->Write();
	//ccc3d->Delete();

	cout << "computing local maxima of ccc of total " << correlation.size() << " points..." << endl;

	// only keep local maximums of correlation
	this->nonLocalMaximumSuppression( rootPoints, correlation );

	typename vector< TinyVector< int, 3 > > :: iterator candidate = rootPoints.begin();
	while ( candidate != rootPoints.end() ){
		vtkTransform * transform = vtkTransform::New();
		transformList.push_back( transform );

		//cout << *candidate << endl;
		// temporary store of cut-out volumes
		this->extractSubVolume( *candidate, S, transform, true );

		volumeList.push_back( Array< short, 3 >() );
		volumeList.back().resize( S.shape() );
		volumeList.back() = S;

		++candidate;
	}

	//for ( int i = 0; i < volumeList.size(); i++ ){
	//	cout << mean(volumeList[i]) << endl;
	//}

	Array< float, 3 > CCC( surfaceImageBlitz.shape() );
	CCC = -1;
	int nsize = 2;

	vtkImageData * output = vtkImageData::New();
	output->SetDimensions( CCC.rows(), CCC.cols(), CCC.depth() );
	output->AllocateScalars();
	output->GetPointData()->GetScalars()->FillComponent(0,-1);
	for ( int i = 0; i < rootPoints.size(); i++ ){
		for ( int x = -nsize; x <= nsize; x++ ){
		for ( int y = -nsize; y <= nsize; y++ ){
		for ( int z = -nsize; z <= nsize; z++ ){
			int xi = rootPoints[i][firstDim];
			int yi = rootPoints[i][secondDim];
			int zi = rootPoints[i][thirdDim];
			output->SetScalarComponentFromDouble( xi+x, yi+y, zi+z, 0, 1 );
		}
		}
		}
		cout << i+1 << " - " << rootPoints[i] << endl;
		//Range I( max( rootPoints[i][0] - nsize, 0 ), min( rootPoints[i][0] + nsize, CCC.ubound(0) ) );
		//Range J( max( rootPoints[i][1] - nsize, 0 ), min( rootPoints[i][1] + nsize, CCC.ubound(1) ) );
		//Range K( max( rootPoints[i][2] - nsize, 0 ), min( rootPoints[i][2] + nsize, CCC.ubound(2) ) );
		//CCC( I, J, K ) = 1; 
	}
	//CCC[ rootPoints ] = 1;

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName("ipath.vtk");
	writer->SetInput( output );
	writer->SetFileTypeToBinary();
	writer->Write();
	writer->Delete();
	cout << "File written." << endl;

	//converter.blitzToVtk( CCC, output );

	//nbfMatlabWriter mwr;
	//mwr.setFileName("ipath");
	//mwr.write(CCC);
	cout << "ipath file written with CCC positions" << endl;
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: nonLocalMaximumSuppression( vector< TinyVector< int, 3 > > & points, vector< Pixel > & correlation ){

	// compute local maxima of correlation in 3D
	typename vector< TinyVector< int, 3 > > :: iterator candidate = points.begin();

	typename vector< Pixel > :: iterator corrs = correlation.begin();
	
	// store current maxima in 1D array
	Array< float, 1 > maxVector( points.size() );

	// initialize to -1
	maxVector = - numeric_limits< float > :: max();

	typename Array< float, 1 > :: iterator maxIterator = maxVector.begin();
	
	// compute local maxima at each point
	while ( candidate != points.end() ){
		for ( int i = 0; i < points.size(); i++ ){
			// if inside ball of given distance
			if ( sqrt( sum( ( points[i] - (*candidate) + 0.0 ) * ( points[i] - (*candidate) + 0.0 ) ) ) < this->neighborhoodSize ){
			//if ( ( fabs( sum( points[i] - (*candidate) + 0.0 ) ) < this->neighborhoodSize ) ){
				if ( correlation[i] > (*maxIterator) ){
					(*maxIterator) = correlation[i];
				}
			}
		}
		++candidate;
		++maxIterator;
	}

	// keep only local maxima
	candidate = points.begin();
	maxIterator = maxVector.begin();
	corrs = correlation.begin();
	while ( candidate != points.end() ){
		if ( ( *corrs < this->correlationTh ) || ( *corrs < (*maxIterator) ) ){
			points.erase(candidate);
			correlation.erase(corrs);
		}
		else{
			++candidate;
			++corrs;
		}
		++maxIterator;
	}
}

template< class Pixel >
void nbfCutSubVolumes< Pixel > :: extractSubVolume( TinyVector< int, 3 > & position, Array< short, 3 > & volume, vtkTransform * transform, bool useFullResolution ){

	// set selected point as center of translation
	float magnificationFactor = 1;
	if ( ( useFullResolution == true ) && ( this->inputFullResolutionImageVtk != NULL ) ){
		int dimensions[3];
		this->inputImageVtk->GetDimensions(dimensions);
		int fullDimensions[3];
		this->inputFullResolutionImageVtk->GetDimensions(fullDimensions);
		magnificationFactor = (float)fullDimensions[0] / (float)dimensions[0];
		this->change->SetInput( this->inputFullResolutionImageVtk );		
		this->change->SetOriginTranslation( - position(0) * magnificationFactor, - position(1) * magnificationFactor, - position(2) * magnificationFactor );
		this->setNotToIgnoreValuesBelowSurface();
		this->change->Update();

		if ( this->templateImageFullResolutionVtk != NULL ){
			this->change2->SetInput( this->templateImageFullResolutionVtk );
			int size[3];
			this->templateImageFullResolutionVtk->GetDimensions(size);
			this->change2->SetOriginTranslation( - size[0] / 2.0, - size[1] / 2.0, - size[2] + this->offsetUnderSurfaceFullResolution );
			this->change2->Update();
			// set relative position of template wrt the image volume	
			this->reslice->SetInformationInput( this->change2->GetOutput() );
		}
	}
	else{
		this->change->SetInput( this->inputImageVtk );		
		this->change->SetOriginTranslation( - position(0), - position(1), - position(2) );
		this->change->Update();
	}

	this->changeSurface->SetOriginTranslation( - position(0), - position(1), - position(2) );

	// get normal from implicit surface representation
	double normal[3];
	
	// get normal as local average of normals in vecinity
	double normalAverage[3];
	normalAverage[0] = 0;
	normalAverage[1] = 0;
	normalAverage[2] = 0;
	int count = 0;
	for ( int i = -2; i <= 2; i++ ){
		for ( int j = -2; j <= 2; j++ ){
			for ( int k = -2; k <= 2; k++ ){
				// VTK gives the negative gradient for some reason
				this->surfaceImageVtk->GetPointGradient( position(0)+i, position(1)+j, position(2)+k, this->surfaceImageVtk->GetPointData()->GetScalars(), normal );
				vtkMath::Normalize(normal);
				count++;
				normalAverage[0] += normal[0];
				normalAverage[1] += normal[1];
				normalAverage[2] += normal[2];
			}
		}
	}
	normal[0] = normalAverage[0];
	normal[1] = normalAverage[1];
	normal[2] = normalAverage[2];
	vtkMath::Normalize(normal);

	//// VTK gives the negative gradient for some reason
	//this->surfaceImageVtk->GetPointGradient( position(0), position(1), position(2), this->surfaceImageVtk->GetPointData()->GetScalars(), normal );
	//vtkMath::Normalize(normal);

	// compute alignment angles with Z and Y axes
	double theta = vtkMath::DegreesFromRadians( - atan2( normal[1], -normal[0] ) );
	double phi = vtkMath::DegreesFromRadians( atan2( sqrt( pow2(normal[0]) + pow2(normal[1]) ), -normal[2] ) );
	
	// The transform object CANNOT be reused (I could not find the reason).
	bool createTransform = false;
	if ( transform == NULL ){
		createTransform = true;
		transform = vtkTransform::New();
	}
	transform->RotateZ(theta);
	transform->RotateY(phi);

	// assume the third angle is zero
	transform->RotateX(0.0);

	// cut subvolume from image
	this->reslice->SetResliceTransform( transform );
	this->reslice->Update();

	nbfVTKInterface::vtkToBlitz( this->reslice->GetOutput(), volume );

	vtkStructuredPointsWriter * writer = vtkStructuredPointsWriter::New();
	writer->SetFileName("tmp1.vtk");
	writer->SetInput( this->reslice->GetOutput() );
	writer->Write();
	writer->Delete();
	cout << "File written." << endl;

	if ( this->ignoreValuesBelowSurface == true ){
		// cut subvolume from surface
		this->resliceSurface->SetResliceTransform( transform );
		this->resliceSurface->Update();	

		Array< float, 3 > implicit;
		nbfVTKInterface::vtkToBlitz( this->resliceSurface->GetOutput(), implicit );

		// kill all components below the surface
		volume = where( implicit < 0, volume, 0 );

		//firstIndex i; secondIndex j;
		//implicit = sqrt( ( i - implicit.rows() / 2 ) + ( j - implicit.cols() / 2 ) );
		//volume = where( implicit < implicit.cols() / 2, volume, 0 );
	}

	//double range[2];
	//this->surfaceImageVtk->GetPointData()->GetScalars()->GetRange();
	//cout << range[0] << ", " << range[1] << endl;
	//cout << min(implicit) << ", " << max(implicit) << endl;

	if ( createTransform == true ){
		transform->Delete();
	}
	else{
		// store 3D position into transform object (as a way to track back the location of the particle)
		transform->Translate( position[0] * magnificationFactor, position[1] * magnificationFactor, position[2] * magnificationFactor );
	}
}
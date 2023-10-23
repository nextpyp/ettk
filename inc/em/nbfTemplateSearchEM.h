/** @file nbfTemplateSearchEM.h
*	Template search inside a volume of interest at specified point locations.
*	Given an image, a template and a list of points, this class will produce the
*   indexes (in the list of input points) where the cross-correlation locally maximizes. 
*	The metric used in the operations has to be provided to the constructor.
*/

#pragma once

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkTransform.h>
#include <vtkImageResample.h>
#include <vtkImageMedian3D.h>

#include <io/nbfVTKInterface.h>
#include <em/nbfWedgedSubImage3D.h>
#include <em/nbfImageMetric.h>

using namespace blitz;

/** Extract sub-volumes at given positions and orientations.
    
	Intended for particle extraction.
	
	@see 
	
	@todo Provide methods for handling correctly.	
*/

template< class Pixel >
class nbfTemplateSearchEM
{
public:

	nbfTemplateSearchEM( nbfImageMetric< Pixel, 3 > * );

	~nbfTemplateSearchEM();

	/// Evaluate CCC at specified locations (using cutting object) and return locations where CCC is locally maximized.
	void execute( vtkPolyData *, nbfWedgedSubImage3D< Pixel > *, vector< int > &, vector< TinyVector< Pixel, 3 > > & );

	/// Set template image to use in calculations.
	void setTemplate( Array< Pixel, 3 > & );
	void setTemplate( vtkImageData * );

	void setTemplates( Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	/// Set correlation threshold: CCC peaks will be considered only if above threshold. Defalut: -1.
	void setCorrelationTh( Pixel in ){ this->correlationTh = in; }

	/// Set local maxima separation range. Only peaks that are at least minimum sepration apart are considered. Default: 20.
	void setMinimumSeparation( int in ){ this->neighborhoodSize = in; }

	/// Set local maxima separation range. Only peaks that are at least minimum sepration apart are considered. Default: 20.
	void setSymmetryThreshold( Pixel in ){ this->symmetryThreshold = in; }

	/// Look for local maxima of CCC only in VOI.
	/// Region is specified as: dimension + lower + upper limits in that dimension.
	void setHeightRange( int a, int b, int c ){ this->dimBound = a; this->hlbound = b; this->hubound = c; }

protected:

	/// Find local maxima of correlation function.
	void nonLocalMaximumSuppression( vtkPolyData *, vector< int > & );

	/// Template reference
	Array< Pixel, 3 > templateImageBlitz;
	vtkImageData * templateImageVtk;
	bool destroyTemplateImage;

	Array< Pixel, 3 > templateStalkBlitz;
	Array< Pixel, 3 > templateHeadBlitz;

	/// Neighborhood size for the computation of local maxima of the correlation function
	int neighborhoodSize;

	/// 3-fold symmetry threshold
	Pixel symmetryThreshold;

	/// Correlation threshold, only keep local maxima above threshold value
	Pixel correlationTh;

	nbfImageMetric< Pixel, 3 > * metric;

	// geometry for constrained search of local maxima
	int dimBound, hlbound, hubound;
};

template< class Pixel >
nbfTemplateSearchEM< Pixel > :: nbfTemplateSearchEM( nbfImageMetric< Pixel, 3 > * m )
: destroyTemplateImage(false)
{
	this->metric = m;
	this->templateImageVtk = NULL;
	this->neighborhoodSize = 20;
	this->symmetryThreshold = 0.2;
	this->correlationTh = -1.0;
	this->dimBound = 0;
	this->hlbound = numeric_limits< int > :: min();
	this->hubound = numeric_limits< int > :: max();
}


template< class Pixel >
nbfTemplateSearchEM< Pixel > :: ~nbfTemplateSearchEM()
{
	if ( this->destroyTemplateImage != NULL ){
		this->templateImageVtk->Delete();
	}
}


template< class Pixel >
void nbfTemplateSearchEM< Pixel > :: setTemplate( Array< Pixel, 3 > & in )
{
	this->templateImageBlitz.reference( in );
	nbfVTKInterface converter;
	this->templateImageVtk = vtkImageData::New();
	converter.blitzToVtk( in, this->templateImageVtk );
	this->destroyTemplateImage = true;
}

template< class Pixel >
void nbfTemplateSearchEM< Pixel > :: setTemplate( vtkImageData * in )
{
	this->templateImageVtk = in;
	nbfVTKInterface converter;
	converter.vtkToBlitz( in, this->templateImageBlitz );
	this->destroyTemplateImage = false;
}

template< class Pixel >
void nbfTemplateSearchEM< Pixel > :: setTemplates( Array< Pixel, 3 > & stalk, Array< Pixel, 3 > & head )
{
	this->templateStalkBlitz.reference( stalk );
	this->templateHeadBlitz.reference( head );
}

template< class Pixel >
void nbfTemplateSearchEM< Pixel > :: execute( vtkPolyData * ccc3d, nbfWedgedSubImage3D< Pixel > * cutter, vector< int > & spikes, vector< TinyVector< Pixel, 3 > > & spikes_normals )
{
	spikes_normals.clear();

	ccc3d->GetPointData()->GetScalars()->SetNumberOfTuples( ccc3d->GetNumberOfPoints() );
	ccc3d->GetPointData()->GetScalars()->SetNumberOfComponents(1);

	vtkPolyData * sym3d = vtkPolyData :: New();
	sym3d->DeepCopy( ccc3d );

	nbfWedgedSubImage3D< Pixel > input1;
	input1.setFixedImage( this->templateImageVtk );
	this->metric->setInput1( &input1 );

	// set geometry of cut volumes to match template size
	TinyVector< int, 3 > d = this->templateImageBlitz.shape();
	cutter->setDimensions( d );

	vtkImageData * cutout = vtkImageData::New();

	cout << "Skipping sections below " << this->hlbound << " and above " << this->hubound << " in dimension " << this->dimBound << endl;

	// compute CCC manually
	//this->metric->imageFilter->execute( this->templateImageVtk );
	Array< float, 3 > B;
	nbfVTKInterface::vtkToBlitz(this->templateImageVtk, B );
	//B = where( B < mean(B), 0, 1 );
	B = B - mean(B);
	B = B / sqrt( sum( B * B ) );

	if ( this->templateStalkBlitz.size() > 0 ){
	
	this->templateStalkBlitz = this->templateStalkBlitz - mean(this->templateStalkBlitz);
	this->templateStalkBlitz = this->templateStalkBlitz / sqrt( sum( this->templateStalkBlitz * this->templateStalkBlitz ) );

	this->templateHeadBlitz = this->templateHeadBlitz - mean(this->templateHeadBlitz);
	this->templateHeadBlitz = this->templateHeadBlitz / sqrt( sum( this->templateHeadBlitz * this->templateHeadBlitz ) );
	}
	
	//Range I( fromStart, toEnd, 2 );
	//Array< float, 3 > Bbin( B(I,I,I) );
	//Bbin = Bbin - mean(Bbin);
	//Bbin = Bbin / sqrt( sum( Bbin * Bbin ) );

	//// ACF computation
	//this->metric->fourierFilter->initializeFFTWhalf( TinyVector< int, 3 >( B.shape() ) );

	//this->metric->fourierFilter->blitzFFTreal = B * this->metric->fourierFilter->shift;
	//fftw_execute( this->metric->fourierFilter->fftplanreal );

	//TinyVector< int, 3 > shape( B.rows(), B.cols(), B.depth() / 2 + 1 );
	//Array< complex< double >, 3 > view( reinterpret_cast< complex<double>* >( this->metric->fourierFilter->blitzFFT.data() ), shape, neverDeleteData );
	//reinterpret_cast< nbfCorrelationImageMetric<Pixel,3>* >(this->metric)->normalizeFourierHalf( view, false );
	//Array< complex< double >, 3 > FFT1( view.shape() );
	//FFT1 = view;

	//view = view * conj( view );
	//fftw_execute( this->metric->fourierFilter->ifftplanreal );

	//Array< double, 3 > acf( B.shape() );
	//acf = this->metric->fourierFilter->blitzFFTreal;
	//acf = acf - mean(acf);
	//acf = acf / sqrt( sum( acf * acf ) );

	nbfMatlabWriter w;
	w.setFileName("p.matlab");
	//w.write(acf);


	nbfWedgedAverageImage3D< Pixel > w1;

	for ( int i = 0; i < ccc3d->GetNumberOfPoints(); i+=1 ){
		
		// retrieve current 3d location
		double point[3];
		ccc3d->GetPoints()->GetPoint(i,point);

		Pixel ccc, sym;
		if ( ( point[this->dimBound] < this->hlbound ) || ( point[this->dimBound] > this->hubound ) ){
			ccc = -1;
			TinyVector< Pixel, 3 > norm(0,0,0);
			spikes_normals.push_back( norm );
		}
		else {

			// retrieve current normal
			double normal[3];	
			ccc3d->GetPointData()->GetNormals()->GetTuple(i,normal);

			//int dims[3];
			//cutter->getVolume()->GetDimensions(dims);

			/////////////
			// WARNING //
			/////////////

			// Set subvolume information (second dimension inverted to compensate ordering of full and half resolution input images)
			//cutter->setPosition( TinyVector< Pixel, 3 >( point[0], dims[1] - point[1], point[2] ) );
			TinyVector< Pixel, 3 > p( point[0], point[1], point[2] );
			cutter->setPosition( p );

			// VTK normal (second component is inverted to compensate as above)
			//cutter->setNormal( TinyVector< Pixel, 3 >( normal[0], -normal[1], normal[2] ) );
			//vtkTransform * t = vtkTransform :: New();
			//Pixel nnormal[3];
			//nnormal[0] = - normal[1];
			//nnormal[1] = normal[0];
			//nnormal[2] = normal[2];
			//t->RotateWXYZ(90,nnormal);
			//t->TransformVector(nnormal,normal);
			//t->Delete();

			// robustify normal computation
			vtkIdList * v = vtkIdList :: New();
			vtkIdList * cells = vtkIdList :: New();
			vtkCell * cell;
			ccc3d->GetPointCells( i, cells );
			for( int j = 0; j < cells->GetNumberOfIds(); j++ ){
				cell = ccc3d->GetCell( cells->GetId(j) );
				for( int k = 0; k < cell->GetNumberOfPoints(); k++ ){
					v->InsertUniqueId( cell->GetPointId(k) );
				}
			}
			cells->Delete();
			double rnormal[3];
			rnormal[0] = rnormal[1] = rnormal[2] = 0;
			double tmp[3];
			for ( int i = 0; i < v->GetNumberOfIds(); i++ ){
				ccc3d->GetPointData()->GetNormals()->GetTuple( v->GetId(i), tmp );
				rnormal[0] += tmp[0];
				rnormal[1] += tmp[1];
				rnormal[2] += tmp[2];
			}
			tmp[0] = rnormal[0] / ( 1.0 * v->GetNumberOfIds() );
			tmp[1] = rnormal[1] / ( 1.0 * v->GetNumberOfIds() );
			tmp[2] = rnormal[2] / ( 1.0 * v->GetNumberOfIds() );
			v->Delete();

			// NEW NORMALS
			//double tmp[3];
			//ccc3d->GetPointData()->GetNormals()->GetTuple( i, tmp );

			// convert normal orientation to Euler angles
			float rotX = vtkMath :: DegreesFromRadians( atan( fabs( tmp[0] / tmp[1] ) ) );
			if ( ( tmp[0] < 0 ) && ( tmp[1] > 0 ) ){
				rotX = - rotX;
			}
			if ( ( tmp[0] < 0 ) && ( tmp[1] < 0 ) ){
				rotX = - ( 180 - rotX );
			}
			if ( ( tmp[0] > 0 ) && ( tmp[1] < 0 ) ){
				rotX = 180 - rotX;
			}
			float rotZ = vtkMath :: DegreesFromRadians( atan( sqrt( tmp[0] * tmp[0] + tmp[1] * tmp[1] ) / fabs( tmp[2] ) ) );
			if ( ( tmp[1] < 0 ) && ( tmp[2] > 0 ) ){
				rotZ = - rotZ;
			}
			if ( ( tmp[1] < 0 ) && ( tmp[2] < 0 ) ){
				rotZ = - ( 180 - rotZ );
			}
			if ( ( tmp[1] > 0 ) && ( tmp[2] < 0 ) ){
				rotZ = 180 - rotZ;
			}
			rotZ = fabs(rotZ);
			normal[0] = rotZ;
			normal[1] = 0;
			normal[2] = rotX;
			/// END NEW NORMALS

			TinyVector< Pixel, 3 > norm( normal[0], normal[1], normal[2] );
			spikes_normals.push_back( norm );
			cutter->setNormal( norm );

			// retrieve sub volume
			cutter->getImage( cutout, (vtkTransform*)NULL, true );
			
			Array< float, 3 > A;

			//nbfVTKInterface::vtkToBlitz(cutout,A);
			//w.write(A);
			//w.write(this->templateImageBlitz);

 		//	this->metric->setInput1( cutter );
 			//this->metric->setInput2( cutter );
			//this->metric->execute();
			////ccc = - sum( A );
			//////ccc = - this->metric->getCorrelationPeak();
			//ccc = this->metric->getCorrelationPeak();

			this->metric->imageFilter->execute( cutout );
			nbfVTKInterface::vtkToBlitz(cutout,A);

			A = A - mean(A);
			A = A / sqrt( sum( A * A ) );

			if ( this->templateStalkBlitz.size() > 0 ){
			
			float d_spike = sum( A * B );
			float d_stalk = sum( A * this->templateStalkBlitz );
			float d_head = sum( A * this->templateHeadBlitz );

#if 0
	w.write(A);
	w.write(B);
	w.write(this->templateStalkBlitz);
	w.write(this->templateHeadBlitz);
#endif


			if ( ( d_spike > d_stalk ) && ( d_spike > d_head ) ){
				ccc = d_spike;
			} else {
				ccc = - numeric_limits< Pixel > :: max();
			}
			} else {
				ccc = sum( A * B );
			}
			
			//float cc1 = sum( A * B );
			//ccc = cc1;
			//cout << cc1 << endl;

			//w.write(A);
			//w.write(B);

			//// signature style particle picking

			//// LCF
			//this->metric->fourierFilter->blitzFFTreal = A * this->metric->fourierFilter->shift;
			//fftw_execute( this->metric->fourierFilter->fftplanreal );

			//reinterpret_cast< nbfCorrelationImageMetric<Pixel,3>* >(this->metric)->normalizeFourierHalf( view, false );

			//view = FFT1 * conj( view );
			//fftw_execute( this->metric->fourierFilter->ifftplanreal );

			//this->metric->fourierFilter->blitzFFTreal = this->metric->fourierFilter->blitzFFTreal - mean(this->metric->fourierFilter->blitzFFTreal);
			//this->metric->fourierFilter->blitzFFTreal = this->metric->fourierFilter->blitzFFTreal / sqrt( sum( this->metric->fourierFilter->blitzFFTreal * this->metric->fourierFilter->blitzFFTreal ) );

			//ccc = ( 1 + cc1 ) / 2.0 * ( 1 + sum( acf * this->metric->fourierFilter->blitzFFTreal ) ) / 2.0;

			//cout << sum( acf * this->metric->fourierFilter->blitzFFTreal ) << endl;

			// normalize FFT2


			//Array< float, 3 > Abin( A(I,I,I) );
			//Abin = Abin - mean(Abin);
			//Abin = Abin / sqrt( sum( Abin * Abin ) );

			//A = where( A < mean(A), 0, 1 );
			//A = A - mean(A);
			//A = A / sqrt( sum( A * A ) );

			//Array< Pixel, 3 > cropA( B.shape() );
			//Range I( A.rows() / 4, A.rows() / 4 + B.rows() );
			//cropA = A( I, I, Range( A.depth() / 2, toEnd ) );
			////w.write(Abin);
			////w.write(Bbin);
			////w.write(A);
			////w.write(cropA);

			//Range I( A.rows() / 2 - B.rows() / 2 + 1, A.rows() / 2 - B.rows() / 2 + B.rows() );
			//Array< Pixel, 3 > cropA( A( I, I, Range( A.depth() / 2, toEnd ) ) );

			//// head
			//Range IsA( A.rows() / 2 - 15, A.rows() / 2 + 15 );
			//Array< Pixel, 3 > stalkA = A( IsA, IsA, Range( fromStart, 6 ) );
			////stalkA = stalkA - mean(stalkA);
			////stalkA = stalkA / sqrt( sum( stalkA * stalkA ) );

			//Range IsB( B.rows() / 2 - 15, B.rows() / 2 + 15 );
			//Array< Pixel, 3 > stalkB = B( IsB, IsB, Range( fromStart, 6 ) );
			////stalkB = stalkB - mean(stalkB);
			////stalkB = stalkB / sqrt( sum( stalkB * stalkB ) );

			//Pixel cccStalk = sum( stalkA * stalkB );
			//
			////w.write(stalkA);
			////w.write(stalkB);

			//Range IhA( A.rows() / 2 - B.rows() / 2 + 1, A.rows() / 2 - B.rows() / 2 + 1 + B.rows() - 1 );
			//Array< Pixel, 3 > headA = A( Range :: all(), Range :: all(), Range( 6 + 1, toEnd ) );
			////headA = headA - mean(headA);
			////headA = headA / sqrt( sum( headA * headA ) );

			//Array< Pixel, 3 > headB = B( Range::all(), Range::all(), Range( 6 + 1, toEnd ) );
			////headB = headB - mean(headB);
			////headB = headB / sqrt( sum( headB * headB ) );

			//Pixel cccHead = sum( headA * headB );
			//
			////ccc = -ccc;

			//if ( ( cccHead > this->correlationTh ) && ( cccStalk > this->correlationTh ) ){
			//	ccc = ( cccHead + cccStalk ) / 2.0;
			//} else {
			//	ccc = -1;
			//}

			//ccc = ( .75 * sum( stalkA * stalkB ) + 1.25 * sum( headA * headB ) ) / 2;
			//w.write(headA);
			//w.write(headB);

			// stalk

			//ccc = sum( cropA * B );
			//ccc = - sum( abs( cropA - B ) ) / B.size();
			//cout << ccc << endl;

			if ( false ){
				// compute 3-foldedness
				norm[1] = 60;
				cutter->setNormal( norm );
				cutter->getImage( cutout, (vtkTransform*)NULL, true );
				this->metric->imageFilter->execute( cutout );

				Array< Pixel, 3 > C;
				nbfVTKInterface::vtkToBlitz( cutout, C );
				C = C - mean(C);
				C = C / sqrt( sum( C * C ) );

				A = where( B < 0, A, mean(A) );
				C = where( B < 0, C, mean(C) );
				//w.write(A);
				//w.write(C);
				sym = - sum( A * C );
			}

			// ccc += sym;
			//ccc = - sum( A );
			//if ( ccc > - numeric_limits< Pixel > :: max() ){
			//	cout << "[ccc] = [" << ccc << "]" << endl;
			//}
		}

		// store CCC value
		ccc3d->GetPointData()->GetScalars()->SetTuple1( i, ccc );
		//sym3d->GetPointData()->GetScalars()->SetTuple1( i, ccc );

		// DEBUG
		if ( ( i % 1000 == 0 ) && ( i > 0 ) && ( ccc != -1 ) ){
		//if ( ( ccc != -1 ) ){
			//cout << i << " : " << ccc << endl;
			cout << i << " of " << ccc3d->GetNumberOfPoints() << endl;
		}
	}

	cutout->Delete();

	// first compute all local maxima of correlation
	vector< int > local_maximas;

	// only keep local maxima of correlation
	Pixel saveSize = this->neighborhoodSize;
	this->neighborhoodSize = 2;
	this->nonLocalMaximumSuppression( ccc3d, local_maximas );
	cout << "Finding local maxima above " << this->correlationTh << " ... ";
	cout << local_maximas.size() << " points extracted." << endl;
	this->neighborhoodSize = saveSize;

	cout << "Eliminating peaks by proximity " << this->neighborhoodSize << "... ";
	// now eliminate spikes that are too close to each other
	spikes.clear();
	Pixel minimum = - numeric_limits< Pixel > :: max();
	for ( int i = 0; i < local_maximas.size(); i++ ){
		// get point coordinates
		double point[3];
		ccc3d->GetPoint( local_maximas[i], point );

		// get correlation value
		Pixel ccc = ccc3d->GetPointData()->GetScalars()->GetTuple1( local_maximas[i] );

		if ( ccc > minimum ){
			minimum = ccc;
		}

		// compute local maxima for each location
		Pixel maxima = - numeric_limits< Pixel > :: max();

		// traverse all other points
		for ( int j = 0; j < local_maximas.size(); j++ ){
			// if point less than given distance away
			if ( vtkMath::Distance2BetweenPoints( point, ccc3d->GetPoint( local_maximas[j] ) ) < this->neighborhoodSize*this->neighborhoodSize ){
				Pixel current = ccc3d->GetPointData()->GetScalars()->GetTuple1( local_maximas[j] );

				if ( current > maxima ){
					maxima = current;
					// if maxima exceed current I'm done (it will not be local maxima)
					if ( maxima > ccc ){
						continue;
					}
				}
			}
		}

		// do not consider points that are too close to the boundary of the search region
		if ( ( ccc >= maxima ) && ( point[this->dimBound] - this->hlbound > 4 ) && ( this->hubound - point[this->dimBound] > 4 ) ){
			spikes.push_back( local_maximas[i] );
		}
	}
	cout << spikes.size() << " spikes left." << endl;

	//cout << "Eliminating peaks by 3-fold symmetry threshold " << this->symmetryThreshold << "... ";
	//vector< int > :: iterator iter = spikes.begin();
	//while ( iter != spikes.end() ){
	//	cout << "sym = " << sym3d->GetPointData()->GetScalars()->GetTuple1( *iter ) << endl;
	//	if ( sym3d->GetPointData()->GetScalars()->GetTuple1( *iter ) < this->symmetryThreshold ){
	//		spikes.erase(iter);
	//	} else {
	//		iter++;
	//	}
	//}
	//cout << spikes.size() << " 3-fold spikes left." << endl;

}

template< class Pixel >
void nbfTemplateSearchEM< Pixel > :: nonLocalMaximumSuppression( vtkPolyData * points, vector< int > & spikes ){

	spikes.clear();

	Pixel minimum = - numeric_limits< Pixel > :: max();

	for ( int i = 0; i < points->GetNumberOfPoints(); i++ ){

		// get point coordinates
		double point[3];
		points->GetPoint(i,point);

		// get correlation value
		Pixel ccc = points->GetPointData()->GetScalars()->GetTuple1(i);

		if ( ccc > minimum ){
			minimum = ccc;
		}

		// if ccc > th and point inside computational domain
		if ( ( ccc > this->correlationTh ) &&
			 ( point[this->dimBound] > this->hlbound ) && ( point[this->dimBound] < this->hubound ) ){

			// compute local maxima for each location
			Pixel maxima = - numeric_limits< Pixel > :: max();
			
			// traverse all other points
			for ( int j = 0; j < points->GetNumberOfPoints(); j++ ){
				// if point less than given distance away
				if ( vtkMath::Distance2BetweenPoints( point, points->GetPoint(j) ) < this->neighborhoodSize*this->neighborhoodSize ){
					Pixel current = points->GetPointData()->GetScalars()->GetTuple1(j);

					if ( current > maxima ){
						maxima = current;
						// if maxima exceed current I'm done (it will not be local maxima)
						if ( maxima > ccc ){
							continue;
						}
					}
				}
			}
			if ( ccc >= maxima ){
				spikes.push_back(i);
			}
		}
	}

	cout << "Global correlation maximum = " << minimum << endl;
}
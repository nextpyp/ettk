#pragma once

#include <io/nbfMrcWriter.h>
#include <em/nbfRadonStructure.h>

#include <em/nbfCorrelationImageMetric.h>
#include <em/nbfFourierImageMetric.h>

/** Top level data structure to handle Radon projections and backprojections.
	This class provides a framework so all iterative algorithms ART-type, 
	SIRT-type, etc. have a trivial implementation.mpi
	A data structure is assembled that contains coeficients (matrix A) and references
	to pixel values (vector x) in which all computations are based. This structure
	has to be computed only once, since A remains fixed throughout the process.
*/
template< class Pixel >
class nbfReconstruction3D : public nbfRadonStructure< Pixel >
{
public:

	// Default constructor
	nbfReconstruction3D(){};

	~nbfReconstruction3D(){};

	// Set output image (where resulting image is stored)
	void setImage( Array< Pixel, 3 > & );

	// Set target projections.
	void setProjections( Array< Pixel, 3 > & );
	void setProjections( Array< short, 3 > & );

	void getProjections( Array< Pixel, 3 > & );
	void getProjections( Array< Pixel, 3 > &, vector< int > & );

	// SIRT iterative technique.
	// Input image is used as initial state and number of iterations is specified in argument.
	void sirt(int);
	void art(int);

	// MPI
	void sirtMPI(int,int); 
	void finalizeMPI(); 
	void slaveMPI();
	static const int tag_done = 1;
	static const int tag_processing = 2;
	static const int parameterSize = 3; // jobNumber, sizeOfObject

	void sirtMedian(int,int=1);

	// SIRT iterative refinement.
	void iterativeRefinement(int,int);

	void progresiveIterativeRefinement( int );

protected:

	// Project volume at tilt angle i
	void project( int, Array< Pixel, 2 > & );

	// Project volume at tilt angle i
	void realignProjections( Array< Pixel, 3 > &, Array< Pixel, 3 > &, vector< int > & );
	void realignProjections( Array< Pixel, 3 > &, Array< Pixel, 3 > & );

	Array< Pixel, 3 > tiltSeries;
	Array< Pixel, 3 > reconstruction;
};

template< class Pixel >
void nbfReconstruction3D< Pixel > :: setImage( Array< Pixel, 3 > & rec )
{
	this->reconstruction.reference( rec );
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: setProjections( Array< short, 3 > & in )
{
	this->tiltSeries.resize( in.shape() );
	this->tiltSeries = cast< Pixel >( in );
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: setProjections( Array< Pixel, 3 > & in )
{
	this->tiltSeries.resize( in.shape() );
	this->tiltSeries = in;
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: project( int angle, Array< Pixel, 2 > & P )
{
	P.resize( this->input.rows(), this->input.cols() );
	for ( int currentSlice = 0; currentSlice < this->input.cols(); currentSlice++ ){
		Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
		nbfRadonStructure< Pixel > :: setImage( rec2D );
		nbfRadonStructure< Pixel > :: project( angle );
		P( Range::all(), currentSlice ) = this->outImage( Range::all(), angle );
	}
}


template< class Pixel >
void nbfReconstruction3D< Pixel > :: sirt( int maxIters )
{
	//nbfTimer t;
	//t.start();
	
	this->inImage.resize( this->tiltSeries.rows(), this->tiltSeries.rows() );
	this->outImage.resize( this->inImage.shape() );
	
	for ( int currentSlice = 0; currentSlice < this->tiltSeries.cols(); currentSlice++ ){
		Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
		this->inImage = rec2D.transpose(secondDim,firstDim);
		Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), currentSlice, Range::all() ) );
		nbfRadonStructure< Pixel > :: setProjections( proj2D );
		nbfRadonStructure< Pixel > :: sirt( maxIters );
		this->reconstruction( Range::all(), currentSlice, Range::all() ) = this->inImage.transpose(secondDim,firstDim);
		//cout << "slice " << currentSlice << " of " << this->tiltSeries.cols() << " done." << endl;
	}

	//t.stop();
	//cout << "Time elapsed (SIRT) = " << t.elapsedSeconds() << endl;
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: art( int maxIters )
{
	//nbfMatlabWriter w;
	//w.setFileName("p.matlab");

	//nbfTimer t;
	//t.start();
	
	this->inImage.resize( this->tiltSeries.rows(), this->tiltSeries.rows() );
	this->outImage.resize( this->inImage.shape() );
	
	for ( int currentSlice = 0; currentSlice < this->tiltSeries.cols(); currentSlice++ ){
		Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
		this->inImage = rec2D.transpose(secondDim,firstDim);
		//w.write( this->inImage );
		Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), currentSlice, Range::all() ) );
		//w.write( proj2D );
		nbfRadonStructure< Pixel > :: setProjections( proj2D );
		nbfRadonStructure< Pixel > :: art( maxIters );
		this->reconstruction( Range::all(), currentSlice, Range::all() ) = this->inImage.transpose(secondDim,firstDim);
		//w.write(this->reconstruction( Range::all(), currentSlice, Range::all() ));
		//cout << "slice " << currentSlice << " of " << this->tiltSeries.cols() << " done." << endl;
	}

	//t.stop();
	//cout << "Time elapsed (ART) = " << t.elapsedSeconds() << endl;
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: sirtMedian( int maxIters, int every )
{
	//nbfTimer t;
	//t.start();
	
	this->inImage.resize( this->tiltSeries.rows(), this->tiltSeries.rows() );
	this->outImage.resize( this->inImage.shape() );

	vtkImageMedian3D * median = vtkImageMedian3D::New();
	vtkImageData * data = vtkImageData::New();
	median->SetKernelSize(3,3,3);

	vtkImageGaussianSmooth * smooth = vtkImageGaussianSmooth::New();
	smooth->SetDimensionality(3);
	smooth->SetRadiusFactor(1);

	for ( int iter = 0; iter < maxIters; iter++ ){
		nbfVTKInterface::blitzToVtk( this->reconstruction, data );
		median->SetInput( data );
		median->Modified();
		median->Update();
		//smooth->SetInput( data );
		//smooth->Modified();
		//smooth->Update();
		Array< Pixel, 3 > C;
		nbfVTKInterface::vtkToBlitzReference( median->GetOutput(), C );
		//nbfVTKInterface::vtkToBlitzReference( smooth->GetOutput(), C );
		this->reconstruction = C;
		for ( int currentSlice = 0; currentSlice < this->tiltSeries.cols(); currentSlice++ ){
			Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
			this->inImage = rec2D.transpose(secondDim,firstDim);
			Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), currentSlice, Range::all() ) );
			nbfRadonStructure< Pixel > :: setProjections( proj2D );
			nbfRadonStructure< Pixel > :: art( every );
			this->reconstruction( Range::all(), currentSlice, Range::all() ) = this->inImage.transpose(secondDim,firstDim);
		}
		cout << "iter " << iter << " done." << endl;
	}

	data->Delete();
	median->Delete();
	smooth->Delete();

	//t.stop();
	//cout << "sirt = " << t.elapsedSeconds() << endl;
}


template< class Pixel >
void nbfReconstruction3D< Pixel > :: getProjections( Array< Pixel, 3 > & projs )
{
	vector< int > indexes;
	for ( int i = 0; i < this->angles.size(); i++ ){
		indexes.push_back(i);
	}
	this->getProjections( projs, indexes );
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: getProjections( Array< Pixel, 3 > & projs, vector< int > & indexes )
{
	projs.resize( this->tiltSeries.shape() );

	this->inImage.resize( this->tiltSeries.rows(), this->tiltSeries.rows() );
	this->outImage.resize( this->inImage.shape() );

	for ( int currentSlice = 0; currentSlice < this->tiltSeries.cols(); currentSlice++ ){
		Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
		this->inImage = rec2D.transpose(secondDim,firstDim);
		Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), currentSlice, Range::all() ) );
		nbfRadonStructure< Pixel > :: setProjections( proj2D );
		if ( this->geometryReady != true ){
			this->buildRadonFramework();
			this->geometryReady = true;
		}
		for ( int i = 0; i < indexes.size(); i++ ){
			nbfRadonStructure< Pixel > :: project( indexes[i] );
		}
		Array< Pixel, 2 > dummy( projs( Range::all(), currentSlice, Range::all() ) );
		nbfRadonStructure< Pixel > :: getProjections( dummy );
		projs( Range::all(), currentSlice, Range::all() ) /= this->circle;
	}
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: iterativeRefinement( int levels, int maxIters )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	Array< Pixel, 3 > projs;

	this->inImage.resize( this->tiltSeries.rows(), this->tiltSeries.rows() );
	// merge input and output image in single reference (for computational efficiency)

	//this->outImage.resize( this->inImage.shape() );
	this->outImage.reference( this->inImage );

	for ( int i = 0; i < levels; i++ ){

		//for ( int tilt = 0; tilt < this->tiltSeries.depth(); tilt++ ){

			this->reconstruction = 0;

			for ( int currentSlice = 0; currentSlice < this->tiltSeries.cols(); currentSlice++ ){

				Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
				this->inImage = rec2D.transpose(secondDim,firstDim);
				Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), currentSlice, Range::all() ) );
				nbfRadonStructure< Pixel > :: setProjections( proj2D );

				// build geometry if neccesary
				if ( this->geometryReady != true ){
					this->buildRadonFramework();
					this->geometryReady = true;
				}


				for ( int iter = 0; iter < maxIters; iter++ ){
					for ( int block = 0; block < this->tiltSeries.depth(); block++ ){		
						// compute forward projection at current angle
						//if ( block != tilt ){
							this->blocks[ block ].project();
							this->blocks[ block ].backProject( 1.0 );
						//}
					}
				}			

				this->reconstruction( Range::all(), currentSlice, Range::all() ) = this->inImage.transpose(secondDim,firstDim);
			}

			//w.write(this->reconstruction);
			//cout << "Reconstruction written. Press any key to continue." << endl;
			//int jolla;
			//cin >> jolla;

			vector< int > positions;
			for ( int ti = 0; ti < this->tiltSeries.depth(); ti++ ){
				positions.push_back(ti);
			}

			//positions.push_back(tilt);

			//this->art( maxIters );
			this->getProjections( projs, positions );

			//w.write( projs );
			//w.write( this->tiltSeries );

			// overwrite tiltSeries with realigned series 
			this->realignProjections( projs, this->tiltSeries, positions );

			//positions.clear();
			//for ( int ti = 0; ti < this->tiltSeries.depth(); ti++ ){
			//	positions.push_back(ti);
			//}
			//this->getProjections( projs, positions );
			//w.write( projs );

		//}

		//w.write(this->tiltSeries);
		//cout << "Realigned tilt series written. Press any key to continue." << endl;
		//cin >> jolla;

		if ( i == levels - 1 ){
		
			Array< Pixel, 3 > A( this->tiltSeries.shape() );
			A = cast< Pixel >( this->tiltSeries );

			vtkImageData * data = vtkImageData::New();
			nbfVTKInterface::blitzToVtk( A, data );
			////vtkStructuredPointsWriter * w = vtkStructuredPointsWriter::New();
			////w->SetFileName("aligned.vtk");
			////w->SetInput( data );
			////w->Write();
			////w->Delete();

			nbfMrcWriter mrcw;
			mrcw.setFileName("aligned.ali");
			mrcw.write(data);
			cout << "realigned stack written" << endl;
			data->Delete();

			w.write(this->tiltSeries );

			this->reconstruction = 0;
			this->art( 3 );

		}

		cout << "step " << i << " done." << endl;
	}
}


template< class Pixel >
void nbfReconstruction3D< Pixel > :: realignProjections( Array< Pixel, 3 > & in1, Array< Pixel, 3 > & in2 )
{
	vector< int > indexes;
	for ( int i = 0; i < this->angles.size(); i++ ){
		indexes.push_back(i);
	}
	this->realignProjections( in1, in2, indexes );
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: realignProjections( Array< Pixel, 3 > & in1, Array< Pixel, 3 > & in2, vector< int > & indexes )
{
	vtkImageData * data2D = vtkImageData::New();
	Array< Pixel, 2 > dummy( in1( Range::all(), Range::all(), 0 ) );
	nbfVTKInterface::blitzToVtk( dummy, data2D );

	nbfImageFilter< double, 2 > imFilter;
	imFilter.setMaskSize( in1.rows() / 2.0 - 4, in1.cols() /2.0 - 4, 0, 2 );
	//imFilter.edgeFilterOn();
	//corr.windowOff();

	nbfFourierFilter< double, 2 > ffFilter;
	//filter.setDimensions( data2D );
	ffFilter.bandPassOn(.05,.25,.025,.01);
	//ffFilter.bandPassOn(.01,.2,.01,.1);
	
	//nbfCorrelationImageMetric< double, 2 > corr( &imFilter, &ffFilter );
	nbfFourierImageMetric< double, 2 > corr( &imFilter, &ffFilter );

	Pixel tolerance = 5;
	corr.setTranslationSearchRestriction( tolerance );
	corr.setToUseMutualCorrelation( true );

	Array< double, 2 > A( in1.rows(), in1.cols() );
	Array< double, 2 > B( A.shape() );

	// transform image if ccc succesful
	vtkImageReslice * reslice = vtkImageReslice::New();
	reslice->SetInterpolationModeToCubic();

	vtkImageData * input1 = vtkImageData::New();
	vtkImageData * input2 = vtkImageData::New();
	//vtkImageData * aligned = vtkImageData::New();

	Pixel globalScore = 0;
	//for ( int i = 0; i < in1.depth(); i++ ){
	for ( int i = 0; i < indexes.size(); i++ ){

		A = cast< double >( in1( Range::all(), Range::all(), indexes[i] ) );
		nbfVTKInterface::blitzToVtk( A, input1 );
		corr.setInput1( input1 );
		
		nbfMatlabWriter writer;
		writer.setFileName("p.matlab");
		//writer.write(A);

		B = cast< double >( in2( Range::all(), Range::all(), indexes[i] ) );
		nbfVTKInterface::blitzToVtk( B, input2 );
		corr.setInput2( input2 );

		reslice->SetBackgroundLevel( mean(B) );
		
		//writer.write(B);
		
		vtkTransform * t1 = vtkTransform :: New();
		corr.executeFourierNewHalf2D(t1);
		double t[3];
		t1->GetPosition(t);
		cout << " " << i << ", t = [" << t[0] << ", " << t[1] << "], ccc = " << corr.getCorrelationPeak() << endl;
		t1->Delete();

		globalScore += corr.getCorrelationPeak();

		vtkTransform * transform = vtkTransform::New();
		transform->Translate( t[0], t[1], 0 );

		// only correct in Y axis
		//transform->Translate( 0, -t[1], 0 );

		if ( sqrt( t[0]*t[0] + t[1]*t[1] ) < tolerance ){
		//if ( fabs( t[1] ) < tolerance ){
			reslice->SetInput( input2 );
			reslice->SetResliceTransform( transform );
			reslice->Update();
			//aligned->DeepCopy( reslice->GetOutput() );
			// write result
			Array< Pixel, 2 > aligned;
			nbfVTKInterface::vtkToBlitz( reslice->GetOutput(), aligned );
			in2( Range::all(), Range::all(), indexes[i] ) = aligned;
			//nbfVTKInterface::vtkToBlitz( aligned,  );
		}
		else{
			cout << "tolerance exceeeded" << endl;
		}
		transform->Delete();

	}

	cout << "Total score = " << globalScore << endl;

	//aligned->Delete();
	input2->Delete();
	input1->Delete();
	reslice->Delete();
	data2D->Delete();
}

template< class Pixel >
void nbfReconstruction3D< Pixel > :: progresiveIterativeRefinement( int maxIters )
{
	nbfMatlabWriter w;
	w.setFileName("p.matlab");

	int zeroTilt = 0;
	for ( int i = 0; i < this->angles.numElements(); i++ ){
		if ( fabs( this->angles(i) ) < fabs( this->angles( zeroTilt ) ) ){
			zeroTilt = i;
		}
	}

	vector< int > currentTilts;
	currentTilts.push_back( zeroTilt );

	int tiltIncrement = 1;

	this->inImage.resize( this->tiltSeries.rows(), this->tiltSeries.rows() );
	// merge input and output image in single reference (for computational efficiency)

	//this->outImage.resize( this->inImage.shape() );
	this->outImage.reference( this->inImage );

	Array< Pixel, 3 > projs;

	// for all tilt angles
	while ( currentTilts.size() < this->angles.size() ){

		// do ART reconstruction
		for ( int currentSlice = 0; currentSlice < this->tiltSeries.cols(); currentSlice++ ){
			Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
			this->inImage = rec2D.transpose(secondDim,firstDim);
			Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), currentSlice, Range::all() ) );
			nbfRadonStructure< Pixel > :: setProjections( proj2D );
			
			// nbfRadonStructure< Pixel > :: art( maxIters );

			// build geometry if neccesary
			if ( this->geometryReady != true ){
				this->buildRadonFramework();
				this->geometryReady = true;
			}

			for ( int iter = 0; iter < maxIters; iter++ ){
				for ( int block = 0; block < currentTilts.size(); block++ ){		
					// compute forward projection at current angle
					this->blocks[ currentTilts[block] ].project();
					this->blocks[ currentTilts[block] ].backProject( 1.0 );
				}
			}			
			
			this->reconstruction( Range::all(), currentSlice, Range::all() ) = this->inImage.transpose(secondDim,firstDim);
			// cout << "slice " << currentSlice << " of " << this->tiltSeries.cols() << " done." << endl;
		}

		w.write( this->reconstruction );

		vector< int > update;
		if ( zeroTilt + tiltIncrement < this->angles.size() ){
			currentTilts.push_back( zeroTilt + tiltIncrement );
			update.push_back( zeroTilt + tiltIncrement );
		}
		if ( zeroTilt - tiltIncrement > -1 ){
			currentTilts.push_back( zeroTilt - tiltIncrement );
			update.push_back( zeroTilt - tiltIncrement );
		}

		tiltIncrement++;

		this->getProjections( projs, update );

		w.write( projs );
		w.write( this->tiltSeries );

		// overwrite tiltSeries with realigned series 
		this->realignProjections( projs, this->tiltSeries, update );

		//cout << "step " << i << " done." << endl;
	}
	w.write( this->tiltSeries );
}


template< class Pixel >
void nbfReconstruction3D< Pixel > :: sirtMPI( int maxIters, int algorithm )
{
	int jobNotSubmitted = 0;
	int jobSubmitted = 1;
	int jobDone = 2;

	Array< int, 1 > queueStatus( this->reconstruction.cols() );
	queueStatus = jobNotSubmitted;

	int distancesToCompute = sum( where( queueStatus == jobNotSubmitted, 1, 0 ) );
	cout << "Total slices to reconstruct = " << distancesToCompute << endl;

	int num_procs;
	MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

	// switch to uni-processor
	if ( num_procs == 1 ){
		for ( int currentSlice = 0; currentSlice < this->reconstruction.cols(); currentSlice++ ){
			this->inImage.resize( this->tiltSeries.rows(), this->tiltSeries.rows() );
			this->outImage.resize( this->inImage.shape() );

			Array< Pixel, 2 > rec2D( this->reconstruction( Range::all(), currentSlice, Range::all() ) );
			this->inImage = rec2D.transpose(secondDim,firstDim);
			Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), currentSlice, Range::all() ) );
			nbfRadonStructure< Pixel > :: setProjections( proj2D );
			switch ( algorithm ){
				case 0:
					nbfRadonStructure< Pixel > :: art( maxIters );
					break;
				case 1:
					nbfRadonStructure< Pixel > :: sirt( maxIters );
			}
			//cout << "[" << min(this->inImage) << "," << max(this->inImage) << "]" << endl;
			cout << "Reconstructing slice " << currentSlice << " of " << this->reconstruction.cols() << endl;
			this->reconstruction( Range::all(), currentSlice, Range::all() ) = this->inImage.transpose(secondDim,firstDim);
		}
	} else {
		// multiprocessor
		int * params = new int[ this->parameterSize ];

		// start all slave processes
		for ( int i = 1; i < num_procs; i++ ){
			if ( sum( where( queueStatus == jobNotSubmitted, 1, 0 ) ) > 0 ){

				// search for next unsent job
				int nextJob = 0;
				while ( nextJob < queueStatus.size() ){
					if ( queueStatus(nextJob) == jobNotSubmitted ){
						break;
					}
					nextJob++;
				}

				cout << "Process " << i << "\treconstructing slice " << nextJob + 1 << "\tof " << queueStatus.size() << endl;

				// update queue status
				queueStatus( nextJob ) = jobSubmitted;

				// build and send first argument
				stringstream objectDataStream;

				Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), nextJob, Range::all() ) );
				//Array< Pixel, 2 > proj2D( this->tiltSeries.rows(), this->tiltSeries.depth() );
				//proj2D = this->tiltSeries( Range::all(), nextJob, Range::all() );

				objectDataStream << maxIters << endl;
				objectDataStream << this->angles << endl;
				objectDataStream << proj2D << endl;

				// send current slice and size of object
				params[0] = nextJob;
				params[1] = objectDataStream.str().size();
				params[2] = algorithm;
				
				// send data size and tag information
				MPI_Send ( params, this->parameterSize, MPI_INT, i, this->tag_processing, MPI_COMM_WORLD );
				
				// send object data and object type information
				MPI_Send ( (char*)objectDataStream.str().c_str(), params[1], MPI_CHAR, i, this->tag_processing, MPI_COMM_WORLD );
			}
		}

		// loop while not all jobs completed
		while ( sum( where( queueStatus != jobDone, 1, 0 ) ) > 0 ){

			MPI_Status status;

			//cout << "MASTER: about to receive control data" << endl;

			// receive notification from slaves
			MPI_Recv( params, this->parameterSize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			// cout << "MASTER: received control data for slice " << result[0] << endl;

			int source = status.MPI_SOURCE;
			int tag = status.MPI_SOURCE;

			// retrieve current index
			int index = params[0];

			//cout << "MASTER: about to receive data of size " << result[1] << " for index " << index << endl;

			stringstream processFile;
			processFile << "mpi_" << source << ".tmp";
			nbfMatlabReader r;
			r.setFileName( processFile.str().c_str() );
			Array< Pixel, 2 > rec2D;
			r.read( rec2D );

			//// allocate array
			//char * objectData = new char[ params[1] ];

			//// receive notification from slaves
			//MPI_Recv( objectData, params[1], MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			//stringstream inputObjectDataStream( objectData );

			//Array< Pixel, 2 > rec2D;
			//inputObjectDataStream >> rec2D;

			////cout << "result =[" << min(rec2D) << "," << max(rec2D) << "]" << endl;

			// assign result
			this->reconstruction( Range::all(), index, Range::all() ) = rec2D.transpose(secondDim,firstDim);

			// update queue status
			queueStatus( index ) = jobDone;

			// send remaining jobs if not already submitted
			if ( sum( where( queueStatus == jobNotSubmitted, 1, 0 ) ) > 0 ){

				// first search for next unsent job
				int nextJob = 0;
				while ( nextJob < queueStatus.size() ){
					if ( queueStatus(nextJob) == jobNotSubmitted ){
						break;
					}
					nextJob++;
				}

				if ( queueStatus( nextJob ) != jobNotSubmitted ){
					cerr << "ERROR - could not find unsent job.\nIn " << __FILE__ << ", line " << __LINE__ << endl;
				} else {
					queueStatus( nextJob ) = jobSubmitted;
				}

				cout << "Process " << source << "\treconstructing slice " << nextJob + 1 << "\tof " << queueStatus.size() << endl;

				stringstream objectDataStream;
				Array< Pixel, 2 > proj2D( this->tiltSeries( Range::all(), nextJob, Range::all() ) );

				objectDataStream << maxIters << endl;
				objectDataStream << this->angles << endl;
				objectDataStream << proj2D << endl;

				// send current (i,j) and size of object
				params[0] = nextJob;
				params[1] = objectDataStream.str().size();
				params[2] = algorithm;

				// send data size and tag information
				MPI_Send ( params, this->parameterSize, MPI_INT, source, this->tag_processing, MPI_COMM_WORLD );
				
				// send object data and object type information
				MPI_Send ( (char*)objectDataStream.str().c_str(), params[1], MPI_CHAR, source, this->tag_processing, MPI_COMM_WORLD );

			}
			//delete [] objectData;
		}
		delete [] params;
	}
}


template< class Pixel >
void nbfReconstruction3D< Pixel > :: finalizeMPI()
{
	int num_procs;
	MPI_Comm_size( MPI_COMM_WORLD, &num_procs );

	int * voidParams = new int[ this->parameterSize ];

	// send termination message to all nodes
	for ( int i = 1; i < num_procs; i++ ){
		MPI_Send ( voidParams, this->parameterSize, MPI_INT, i, this->tag_done, MPI_COMM_WORLD );
	}
	delete [] voidParams;
}

template< class Pixel>
void nbfReconstruction3D< Pixel > :: slaveMPI()
{
	int source = 0;
	int tag = this->tag_processing;

	char * objectData;

	MPI_Status status;

	int * params = new int[ this->parameterSize ];

	while ( tag != this->tag_done ){

		//cout << "SLAVE: about to receive control data" << endl;

		// retrieve object size and tag information
		MPI_Recv( params, this->parameterSize, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		
		tag = status.MPI_TAG;

		//cout << "SLAVE: about to receive data of size " << params[1] << " for slice " << params[0] + 1 << endl;

		if ( tag != this->tag_done ){
			
			// allocate array
			objectData = new char[ params[1] ];

			// retrieve object data
			MPI_Recv( objectData, params[1], MPI_CHAR, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			//cout << "SLAVE: data received succesfully for slice " << params[0] + 1 << endl;

			stringstream inputString( objectData );

			Array< Pixel, 1 > anglesParam;
			Array< Pixel, 2 > proj2D;
			int maxIters;
			inputString >> maxIters;
			inputString >> anglesParam;
			inputString >> proj2D;
			
			//Array< Pixel, 2 > rec2D( proj2D.rows(), proj2D.rows() );
			//rec2D = 0;
			//cout << "angles = " << angles.shape() << endl;
			//cout << "rec2D = " << rec2D.shape() << endl;
			//cout << "proj2D = " << proj2D.shape() << endl;
			//cout << "maxIters = " << maxIters << endl;

			// do processing
			//cout << "SLAVE: processing entry " << params[0] + 1 << endl;
			if ( this->angles.size() == 0 ){
				// cout << "Setting angles for the first time" << endl;
				this->setAngles( anglesParam );
				this->inImage.resize( proj2D.rows(), proj2D.rows() );
				this->outImage.resize( this->inImage.shape() );
			}

			//this->inImage = rec2D.transpose(secondDim,firstDim);
			this->inImage = 0;

			//cout << "this->inImage: " << this->inImage.shape() << ",[" << min(this->inImage) << "," << max(this->inImage) << "]" << endl;
			//cout << "proj2D: " << proj2D.shape() << ",[" << min(proj2D) << "," << max(proj2D) << "]" << endl;
			//cout << "SLAVE: processing entry before exec " << params[0] + 1 << endl;
			nbfRadonStructure< Pixel > :: setProjections( proj2D );
			//cout << "SLAVE: processing entry between exec " << params[0] + 1 << endl;
			
			int algorithm = params[2];

			switch ( algorithm ){
				case 0:
					nbfRadonStructure< Pixel > :: art( maxIters );
					break;
				case 1:
					nbfRadonStructure< Pixel > :: sirt( maxIters );
					break;
			}

			//cout << "SLAVE: processing entry after exec " << params[0] + 1 << endl;
			
			//  Get this processes's rank.
			int my_id;
			MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

			// we save the result to disk
			stringstream currentFile;
			currentFile << "mpi_" << my_id << ".tmp";
			nbfMatlabWriter w;
			w.setFileName( currentFile.str().c_str() );
			w.write( this->inImage );

			//stringstream outputString;
			//outputString << this->inImage << endl;
			////outputString << rec2D << endl;

			////cout << "SLAVE: this->inImage =[" << min(this->inImage) << "," << max(this->inImage) << "], " << this->inImage.shape() << endl;
			//params[0] = params[0];
			//params[1] = outputString.str().size();
			////cout << "SLAVE: about to send control result of size " << result[1] << "(" << this->inImage.size() << ")," << result[0] + 1 << endl;

			// send result and notify master node we have finished processing
			MPI_Send ( params, this->parameterSize, MPI_INT, 0, this->tag_processing, MPI_COMM_WORLD );

			//cout << "SLAVE: about to send data result of size " << result[1] << "(" << this->inImage.size() << "), " << result[0] + 1 << endl;

			// send reconstructed slice
			//cout << "sending result back to source... size " << params[1] << endl;
			//MPI_Send ( (char*)outputString.str().c_str(), params[1], MPI_CHAR, 0, this->tag_processing, MPI_COMM_WORLD );

			//cout << "SLAVE: finished processing entry " << result[0] + 1 << endl;
			delete [] objectData;
		}
	}
	delete [] params;
}
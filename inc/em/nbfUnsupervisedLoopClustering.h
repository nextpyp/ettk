#pragma once

//using namespace blitz;

#include <em/nbfLoopClustering.h>
#include <em/nbfHierarchicalClustering.h>
#include <em/nbfFourierImageMetricCore.h>

/** Iterative Clustering Method "Loop".
*/
template< class Pixel >
class nbfUnsupervisedLoopClustering : public nbfLoopClustering< Pixel >
{
public:

	nbfUnsupervisedLoopClustering();

	~nbfUnsupervisedLoopClustering(){};

protected:

	void doClustering( Array< Pixel, 2 > & );

	void normalDistribution( int, Pixel, Array< Pixel, 1 > & );

	void alignVolumesToReference(int);

	// store variances
	vector< Pixel > variances;

	vector< TinyVector< Pixel, 2 > > referencey;
};

template< class Pixel >
nbfUnsupervisedLoopClustering< Pixel > :: nbfUnsupervisedLoopClustering()
: nbfLoopClustering< Pixel >()
{
}

template< class Pixel >
void nbfUnsupervisedLoopClustering< Pixel > :: doClustering( Array< Pixel, 2 > & classification )
{
	//nbfMatlabReader reader;
	//reader.setFileName("points.matlab");
	//Array< Pixel, 2 > y;
	//reader.read(y);

	// retrieve saved state
	stringstream currentFile;
	currentFile << this->fileHeader;
	this->loadState( currentFile );

	// if empty, run initial classification
	if ( this->classes.rows() == 0 ){

		// compute initial references by thresholding dendrogram
		// and enforcing self-consistency
		cout << "Computing initial refernces..." << endl;
		this->computeInitialReferences();

		// bundle align all references with MPI
		this->bundleAlignment( this->references, 5 );

		// save state
		currentFile.clear();
		currentFile << this->fileHeader << "_initial";
		this->saveState( currentFile );

		return;
	}
	

	// Align volumes to new references
	this->alignVolumesToReferences();

	//reader.setFileName("randindex.matlab");
	//Array< Pixel, 2 > randindex;
	//reader.read(randindex);
	//for ( int i = 0; i < randindex.cols(); i++ ){
	//	this->referencey.push_back( TinyVector< Pixel, 2 >( y( 0, (int)randindex(0,i) - 1 ), y( 1, (int)randindex(0,i) - 1 ) ) );
	//}

	// by now we have initial references stored in this->references

	// LOOP INITIALIZATION

	// minimization tolerance
	Pixel th = 1e-4;

	// minimum acceptable number of classes
	int kmin = 1;

	// store cost function evolution
	vector< Pixel > evolutionLikelihood;

	// initialize variances
	this->variances.resize( this->references.size() );
	//this->variances.resize( this->referencey.size() );

	//reader.setFileName("estcov.matlab");
	//reader.read(randindex);

	// ... and initialize to some arbitrary value
	for ( int i = 0; i < this->variances.size(); i++ ){
		this->variances[i] = .01;
		//this->variances[i] = randindex(0,i);
	}

	// get dimensionality of data (number of voxels in sub-volume)
	Pixel data_dims = product( this->input[0].getDimensions() );
	//data_dims = y.rows();
	data_dims = 0.0;
	
	// dimension of parameter space THETA ( mean + variance )
	Pixel theta_dims = data_dims + 1; 

	// kmax is the initial number of mixture components
	Pixel current_num_references = this->references.size();
	//Pixel current_num_references = this->referencey.size();

	// the initial estimates of the alpha probabilities are set to: 1.0 / current_num_references
	Array< Pixel, 1 > alpha( this->references.size() );
	alpha = 1.0 / this->references.size();
	//Array< Pixel, 1 > alpha( this->referencey.size() );
	//alpha = 1.0 / this->referencey.size();

	// w_m^i - memberships will contain the assignments of each data point to the mixture components, as result of the E-step
	Array< Pixel, 2 > w_mi( this->references.size(), this->input.size() );
	//Array< Pixel, 2 > w_mi( this->referencey.size(), y.cols() );
	w_mi = 0;

	// u_m^i - probabilities will contain the normal distributions
	Array< Pixel, 2 > u_mi( w_mi.shape() );
	u_mi = 0;

	// having the initial means, covariances, and probabilities, we can initialize the indicator functions following the standard EM equation
	// Notice that these are unnormalized values
	for ( int i = 0; i < current_num_references; i++ ){
		Array< Pixel, 1 > dummy( u_mi( i, Range::all() ) );
		this->normalDistribution( i, data_dims, dummy );
		w_mi( i, Range::all() ) = u_mi( i, Range::all() ) * alpha(i);
		//cout << w_mi( i, Range::all() ) << endl;
	}

	// we can use the memberships variables (unnormalized) to compute the log-likelihood and store it for later plotting its evolution
	// we also compute and store the number of components
	int iteration_counter = 0;

	// store log-likelihood evolution
	vector< Pixel > loglike;

	firstIndex i; secondIndex j;
	Pixel aux = sum( log( numeric_limits< Pixel > :: min() + sum( w_mi(j,i), j ) ) );
	loglike.push_back(aux);

	// compute and store total likelihood
	aux = - loglike[iteration_counter] + ( theta_dims / 2.0 * sum( log(alpha) ) ) + ( theta_dims / 2.0 + 0.5 ) * current_num_references * log( 1.0 * this->input.size() );
	//aux = - loglike[iteration_counter] + ( theta_dims / 2.0 * sum( log(alpha) ) ) + ( theta_dims / 2.0 + 0.5 ) * current_num_references * log( 1.0 * y.cols() );
	evolutionLikelihood.push_back( aux );
	
	// store best configuration
	Array< Pixel, 1 > bestAlpha;
	vector< nbfWedgedAverageImage3D< Pixel > >  bestReferences;
	//vector< TinyVector< Pixel, 2 > > bestReferencey;
	vector< Pixel > bestVariances;
	Array< Pixel, 3 > bestAlignmentToReferences;

	// minimum description length seen so far, and corresponding parameter estimates
	Pixel bestEvolutionLikelihood = evolutionLikelihood[iteration_counter];

	bestAlpha.resize( alpha.shape() );
	bestAlpha = alpha;
	bestReferences = this->references;
	//bestReferencey = this->referencey;
	bestVariances = this->variances;
	bestAlignmentToReferences.resize( this->alignmentToReferences.shape() );
	bestAlignmentToReferences = this->alignmentToReferences;

	// variable for the outer loop
	bool reached_minimum_number_of_classes = false;

	// OUTER LOOP (ITERATE FROM KMAX TO KMIN NUMBER OF REFERENCES)

	// the outer loop will take us down from kmax to kmin components
	while ( reached_minimum_number_of_classes == false ){

		// INNER LOOP (ITERATE UNTIL TOLERANCE REACHED)

		// auxiliary variable of the inner loop
		bool reached_tolerance = false;

		// this inner loop is the component-wise EM algorithm with the modified M-step that can kill components
		while ( reached_tolerance == false ){

			// traverse for all references: [0,current_num_references]
			int current_reference = 0;

			while ( current_reference < current_num_references ){

				w_mi.resize( current_num_references, this->input.size() );
				//w_mi.resize( current_num_references, y.cols() );

				// we start with the M step
				// first, we compute the indicator function: w_m^i = u_m^i * \alpha_m
				for ( int i = 0; i < current_num_references; i++ ){
					w_mi( i, Range::all() ) = u_mi( i, Range::all() ) * alpha(i);
					// truncate weights if too small wrt principal component
					// w_mi( i, Range::all() ) = where( w_mi( i, Range::all() ) > max( w_mi( i, Range::all() ) ) / 20.0, w_mi( i, Range::all() ), 0 );
				}

				// compute normalized memberships (so contributions from each volume add-up to 1)
				Array< Pixel, 2 > normalized_w_mi( w_mi.shape() );

				for ( int i = 0; i < w_mi.cols(); i++ ){
					normalized_w_mi( Range::all(), i ) = w_mi( Range::all(), i ) / ( numeric_limits< Pixel > :: min() + sum( w_mi( Range::all(), i ) ) );
				}

				Pixel normalize = 1.0 / sum( normalized_w_mi( current_reference, Range::all() ) );

				// now we perform the standard M-step for mean, covariance and alignments
#if 1
				// estimate new mean
				Array< Pixel, 2 > alignments2( this->input.size(), 17 );
				// set weights as first components (normalization is done inside nbfWedgedAverageImage3D)
				alignments2( Range::all(), 0 ) = normalized_w_mi( current_reference, Range::all() ); // this->classes( i, Range::all() );
				// truncate weights if too small wrt principal component
				alignments2( Range::all(), 0 ) = where( alignments2( Range::all(), 0 ) > max( alignments2( Range::all(), 0 ) ) / 20.0, alignments2( Range::all(), 0 ), 0 );
				// set alignments
				alignments2( Range::all(), Range(1,toEnd) ) = this->alignmentToReferences( current_reference, Range::all(), Range(3,toEnd) );				
				// set to new weights
				this->references[ current_reference ].setAlignments( alignments2 );

				// cout << "New weights for reference " << current_reference << " = " << alignments2( Range::all(), 0 ) << endl;

				// estimate exponential decay
				//this->variances[ current_reference ] = normalize * sum( normalized_w_mi( current_reference, Range::all() ) * pow2( this->alignmentToReferences( current_reference, Range::all(), 0 ) ) ) / data_dims;
				this->variances[ current_reference ] = normalize * sum( normalized_w_mi( current_reference, Range::all() ) * pow2( this->alignmentToReferences( current_reference, Range::all(), 0 ) ) );

				//cout << "New variance for reference " << current_reference << " = " << this->variances[ current_reference ] << endl;

				// compute new alignments
				this->alignVolumesToReference( current_reference );
#else
				// 1. estimate mean
				Array< Pixel, 2 > aux( data_dims, y.cols() );
				// update this->references[current_reference]
				for ( int i = 0; i < data_dims; i++ ){
					aux( i, Range::all() ) = normalized_w_mi( current_reference, Range::all() ) * y( i, Range::all() );
				}
				// estmu(:,current_reference) = normalize * sum(aux,2);
				firstIndex i; secondIndex j;
				Array< Pixel, 1 > A(2);
				A = normalize * sum( aux, j );
				this->referencey[ current_reference ] = TinyVector< Pixel, 2 >( A(0), A(1) );

				// 2. estimate variance
				aux = aux * y;
				Array< Pixel, 1 > M(2);
				M(0) = this->referencey[ current_reference ][0];
				M(1) = this->referencey[ current_reference ][1];
				this->variances[ current_reference ] = mean( normalize * sum( aux, j ) - pow2( M ) );
				//variances[ current_reference ] = mean( diag( normalize * diag(sum(aux.*y,2)) - diag(estmu(:,current_reference).^2) ) );
				
				//cout << "m=" << this->referencey[ current_reference ] << ",v=" << this->variances[ current_reference ] << endl;
#endif
				// 3. update alignment of input volumes to references[i]

				// this is the special part of the M step that is able to kill components
				// cout << sum( normalized_w_mi( current_reference, Range::all() ) ) << ", " << theta_dims / 2.0 << endl;

				alpha(current_reference) = max( sum( normalized_w_mi( current_reference, Range::all() ) ) - theta_dims / 2.0, 0.0 ) / this->input.size();
				//alpha(current_reference) = max( sum( normalized_w_mi( current_reference, Range::all() ) ) - theta_dims / 2.0, 0.0 ) / y.cols();
				alpha = alpha / sum(alpha);

				// this is an auxiliary variable that will be used to signal the killing of the current component being updated
				bool reference_killed = false;

				// we now have to do some book-keeping if the current component was killed that is, we have to rearrange the vectors and matrices that store the parameter estimates
				if ( alpha(current_reference) == 0 ){

					cout << "Reference " << current_reference << " killed." << endl;
					reference_killed = true;

					// since we've just killed a component, current_num_references must decrease
					current_num_references = current_num_references - 1;

					// re-assign references, variances, alpha, u_mi and alignments

					this->references.erase( this->references.begin() + current_reference );
					//this->referencey.erase( this->referencey.begin() + current_reference );
					this->variances.erase( this->variances.begin() + current_reference );

					Array< Pixel, 1 > newAlpha( current_num_references );
					Array< Pixel, 2 > new_u_mi( current_num_references, this->input.size() );
					//Array< Pixel, 2 > new_u_mi( current_num_references, y.cols() );

					Array< Pixel, 3 > newAlignmentToReferences( current_num_references, this->input.size(), 19 );
					//Array< Pixel, 3 > newAlignmentToReferences( current_num_references, y.cols(), 19 );

					int count = 0;
					for ( int i = 0; i < alpha.rows(); i++ ){
						if ( i != current_reference ){
							newAlpha(count) = alpha(i);
							new_u_mi( count, Range::all() ) = u_mi( i, Range::all() );
							newAlignmentToReferences( count, Range::all(), Range::all() ) = this->alignmentToReferences( i, Range::all(), Range::all() );
							count++;
						}
					}

					alpha.resize( newAlpha.shape() );
					alpha = newAlpha;

					u_mi.resize( new_u_mi.shape() );
					u_mi = new_u_mi;

					this->alignmentToReferences.resize( newAlignmentToReferences.shape() );
					this->alignmentToReferences = newAlignmentToReferences;
				}

				if ( reference_killed == false ){
					// if the component was not killed, we update the corresponding indicator variables...
					Array< Pixel, 1 > dummy( u_mi( current_reference, Range::all() ) );
					this->normalDistribution( current_reference, data_dims, dummy );
					// ...and go on to the next component
					current_reference = current_reference + 1;  
				}
				// if ( reference_killed == true ), it means the in the position "current_reference", we now have a component that was not yet visited in this sweep, and
				// so all we have to do is go back to the M setp without increasing "current_reference"

			} // this is the end of the innermost "while current_reference < current_num_references" loop which cycles through the components

			cout << "alpha(" << iteration_counter << ")= " << alpha << endl;
			cout << "variances(" << iteration_counter << ")= ";
			for ( int i = 0; i < this->variances.size(); i++ ){
				cout << this->variances[i] << " ";
			}
			cout << endl;

			// increment the iterations counter            
			iteration_counter = iteration_counter + 1;

			w_mi.resize( current_num_references, this->input.size() );
			//w_mi.resize( current_num_references, y.cols() );
			//u_mi.resize( w_mi.shape() );

			for ( int i = 0; i < current_num_references; i++ ){
				Array< Pixel, 1 > dummy( u_mi( i, Range::all() ) );
				this->normalDistribution( i, data_dims, dummy );
				w_mi( i, Range::all() ) = u_mi( i, Range::all() ) * alpha(i);
			}

			if ( current_num_references != 1 ){
				// if the number of surviving components is not just one, we compute the loglikelihood from the un-normalized assignment variables
				firstIndex i; secondIndex j;
				loglike.push_back( sum( log( numeric_limits< Pixel > :: min() + sum( w_mi(j,i), j ) ) ) );
			} else {
				// if it is just one component, it is even simpler
				loglike.push_back( sum( log ( numeric_limits< Pixel > :: min() + w_mi ) ) );
			}

			// compute and store the description length and the current number of components
			evolutionLikelihood.push_back( - loglike[iteration_counter] + ( theta_dims / 2.0 * sum( log(alpha) ) ) + ( theta_dims / 2.0 + 0.5 ) * current_num_references * log( 1.0 * this->input.size() ) );
			//evolutionLikelihood.push_back( - loglike[iteration_counter] + ( theta_dims / 2.0 * sum( log(alpha) ) ) + ( theta_dims / 2.0 + 0.5 ) * current_num_references * log( 1.0 * y.cols() ) );

			cout << "ML(" << iteration_counter << ") = " << evolutionLikelihood[iteration_counter] << endl;

			// compute the change in loglikelihood to check if we should stop
			Pixel deltlike = loglike[ iteration_counter ] - loglike[ iteration_counter - 1 ];

			// if the relative change in loglikelihood is below the threshold, we stop CEM2
			if ( fabs( deltlike / loglike[ iteration_counter - 1 ] ) < th ){
				reached_tolerance = true;
				cout << "Iteration " << iteration_counter << ", tolerance reached." << endl;
			}

		} // this end is of the inner loop: "while(cont)"

		// now check if the latest description length is the best; 
		// if it is, we store its value and the corresponding estimates 
		if ( evolutionLikelihood[iteration_counter] < bestEvolutionLikelihood ){
			bestAlpha.resize( alpha.shape() );
			bestAlpha = alpha;
			bestReferences = this->references;
			//bestReferencey = this->referencey;
			bestVariances = this->variances;
			bestEvolutionLikelihood = evolutionLikelihood[ iteration_counter ];
			bestAlignmentToReferences.resize( this->alignmentToReferences.shape() );
			bestAlignmentToReferences = this->alignmentToReferences;
		}

		// at this point, we may try smaller mixtures by killing the component with the smallest alpha probability and then restarting CEM2,
		// as long as current_num_references is not yet at kmin
		if ( current_num_references > kmin ){

			TinyVector< int, 1 > minAlphaIndex = minIndex(alpha);

			// what follows is the book-keeping associated with removing one component

			cout << "Reference " << minAlphaIndex[0] << " killed by least probability criteria." << endl;

			current_num_references = current_num_references - 1;

			this->references.erase( this->references.begin() + minAlphaIndex[0] );
			//this->referencey.erase( this->referencey.begin() + minAlphaIndex[0] );
			this->variances.erase( this->variances.begin() + minAlphaIndex[0] );

			Array< Pixel, 1 > newAlpha( current_num_references );
			Array< Pixel, 2 > new_u_mi( current_num_references, this->input.size() );
			//Array< Pixel, 2 > new_u_mi( current_num_references, y.cols() );
			Array< Pixel, 3 > newAlignmentToReferences( current_num_references, this->input.size(), 19 );
			//Array< Pixel, 3 > newAlignmentToReferences( current_num_references, y.cols(), 19 );

			int count = 0;
			for ( int i = 0; i < alpha.rows(); i++ ){
				if ( i != minAlphaIndex[0] ){
					newAlpha(count) = alpha(i);
					new_u_mi( count, Range::all() ) = u_mi( i, Range::all() );
					newAlignmentToReferences( count, Range::all(), Range::all() ) = this->alignmentToReferences( i, Range::all(), Range::all() );
					count++;
				}
			}

			alpha.resize( newAlpha.shape() );
			alpha = newAlpha;

			u_mi.resize( new_u_mi.shape() );
			u_mi = new_u_mi;

			this->alignmentToReferences.resize( newAlignmentToReferences.shape() );
			this->alignmentToReferences = newAlignmentToReferences;

			// we re-normalize the alpha probabilities after killing the component
			alpha = alpha / sum(alpha);

			// increment the iterations counter
			iteration_counter = iteration_counter + 1;

			// ...and compute the log-likelihhod function and the description length
			//clear memberships;
			//clear probabilities;

			w_mi.resize( current_num_references, this->input.size() );
			//w_mi.resize( current_num_references, y.cols() );

			for ( int i = 0; i < current_num_references; i++ ){
				Array< Pixel, 1 > dummy( u_mi( i, Range::all() ) );
				this->normalDistribution( i, data_dims, dummy );
				w_mi( i, Range::all() ) = u_mi( i, Range::all() ) * alpha(i);
			}
			
			if ( current_num_references != 1 ){
				firstIndex i; secondIndex j;
				loglike.push_back( sum( log( numeric_limits< Pixel > :: min() + sum(w_mi(j,i),j) ) ) );
			} else {
				loglike.push_back( sum( log( numeric_limits< Pixel > :: min() + w_mi ) ) );
			}
			evolutionLikelihood.push_back( - loglike[iteration_counter] + ( theta_dims / 2.0 * sum( log(alpha) ) ) + ( theta_dims / 2.0 + 0.5 ) * current_num_references * log( 1.0 * this->input.size() ) );
			//evolutionLikelihood.push_back( - loglike[iteration_counter] + ( theta_dims / 2.0 * sum( log(alpha) ) ) + ( theta_dims / 2.0 + 0.5 ) * current_num_references * log( 1.0 * y.cols() ) );

			cout << "ML(" << iteration_counter << ") = " << evolutionLikelihood[iteration_counter] << endl;

		} else { // this else corresponds to "if current_num_references > kmin"
			// of course, if current_num_references is not larger than kmin, we must stop      
			reached_minimum_number_of_classes = true;
		}
	} // this is the end of the outer loop "while(k_cont)" 

	//for ( int i = 0; i < bestReferencey.size(); i++ ){
	//	cout << "Mean " << i << " = " << bestReferencey[i] << endl;
	//}
	
	cout << "Best k = " << bestReferences.size() << endl;

	nbfMatlabWriter w;
	w.setFileName("evolution.matlab");
	Array< Pixel, 1 > evolution( evolutionLikelihood.size() );
	for ( int i = 0; i < evolution.rows(); i++ ){
		evolution(i) = evolutionLikelihood[i];
	}
	w.write(evolution);

	cout << evolution << endl;
	cout << "Evolution likelihood saved in file evolution.matlab" << endl;

	//this->referencey = bestReferencey;
	this->references = bestReferences;
	this->alignmentToReferences.resize( bestAlignmentToReferences.shape() );
	this->alignmentToReferences = bestAlignmentToReferences;

	// save state
	currentFile.clear();
	currentFile << this->fileHeader << "_iteration_" << iteration_counter;
	this->saveState( currentFile );

	cout << "Best classification state saved with prefix: " << currentFile.str() << endl;
}

template< class Pixel >
void nbfUnsupervisedLoopClustering< Pixel > :: normalDistribution( int i, Pixel dimensions, Array< Pixel, 1 > & out )
{
#if 0
	Array< Pixel, 2 > y;
	nbfMatlabReader reader;
	reader.setFileName("points.matlab");
	reader.read(y);

	//int dim = 2;
	//[dim npoints] = size(x);
	Pixel dd = pow2( this->variances[i] );
	//in = inv(covar+ realmin*eye(dim));
	Pixel ff = pow( 2.0f * vtkMath::Pi() , - dimensions / 2.0f ) / this->variances[i];
	//quadform = zeros(1,npoints);
	//centered = (x-m*ones(1,npoints));
	//if dim ~= 1
	//	y = ff * exp(-0.5*sum(centered.*(in*centered)));
	//else
	//	y = ff * exp(-0.5*in*centered.^2 );
	//end
	for ( int j = 0; j < y.cols(); j++ ){
		//cout << y( Range::all(), j ) << endl;
		//cout << this->referencey[i] << endl;
		out(j) = ff * exp( -.5 * ( pow2( y( 0, j ) - this->referencey[i][0] ) + pow2( y( 1, j ) - this->referencey[i][1] ) ) / this->variances[i] );
		//cout << out(j) << endl;
	}
#else
	//out = pow( 2.0 * vtkMath::Pi(), - dimensions / 2.0 ) / sqrt( pow( this->variances[i], dimensions ) ) * exp( - .5 / this->variances[i] * this->alignmentToReferences( i, Range::all(), 0 ) * dimensions );
	out = 1.0 / sqrt( 2.0 * vtkMath::Pi() * this->variances[i] ) * exp( - .5 / this->variances[i] * pow2( this->alignmentToReferences( i, Range::all(), 0 ) ) );
	//cout << this->alignmentToReferences( i, Range::all(), 0 ) << endl;
	//cout << out << endl;
#endif
}
		

template< class Pixel >
void nbfUnsupervisedLoopClustering< Pixel > :: alignVolumesToReference( int reference )
{
	// reset transforms before computing distances
	for ( int i = 0; i < this->input.size(); i++ ){
		this->input[i].setTransform( (vtkTransform*)NULL );
	}

	// initialize new alignment matrix of references to volumes
	Array< Pixel, 3 > aligs( 1, this->input.size(), 19 );
	aligs = -1;

	//this->alignmentToReferences.resize( this->classes.rows(), this->input.size(), 19 );
	//this->alignmentToReferences = -1;

	vector< nbfWedgedAverageImage3D< Pixel > > dummy;
	dummy.push_back( this->references[reference] );

	this->metric->getDistances( dummy, this->input, aligs );

	// compute distance to references as weighted-mean distance to elements

	vector< nbfWedgedSubImage3D< Pixel > > list1, list2;
	vector< TinyVector< int, 2 > > positions;

	// for all reference components
	for ( int j = 0; j < aligs.cols(); j++ ){

		if ( this->references[reference].weights(j) > 0 ){

			list1.push_back( this->references[reference].getVolumes()[j] );

			// for all input volumes
			for ( int k = 0; k < aligs.cols(); k++ ){

				double matrix[16];
				for ( int m = 0; m < 16; m++ ){
					matrix[m] = aligs( 0, k, 3 + m );
				}
				vtkMatrix4x4 * mat = vtkMatrix4x4::New();
				mat->DeepCopy( matrix );

				this->input[k].setTransform( mat );
				mat->Delete();

				list2.push_back( this->input[k] );

				TinyVector< int, 2 > currentPosition( list1.size() - 1, list2.size() - 1 );
				positions.push_back( currentPosition );
			}
		}
	}

	vector< nbfWedgedImage3D< Pixel > * > plist1, plist2;

	for ( int j = 0; j < list1.size(); j++ ){
		plist1.push_back( &(list1[j]) );
	}

	for ( int j = 0; j < list2.size(); j++ ){
		plist2.push_back( &(list2[j]) );
	}

	Array< Pixel, 2 > allDistances( positions.size(), 19 );
	allDistances = -1;
	this->metric->getDistances( plist1, plist2, positions, allDistances, 2 ); // use executeFourierNew

	int count = 0;

	Array< Pixel, 2 > distances( aligs.cols(), aligs.cols() );
	distances = 0;
	for ( int j = 0; j < aligs.cols(); j++ ){
		if ( this->references[reference].weights(j) > 0 ){
			for ( int k = 0; k < aligs.cols(); k++ ){
				distances(j,k) = allDistances(count,0);
				count++;
			}
		}
	}
	for ( int j = 0; j < aligs.cols(); j++ ){
		this->alignmentToReferences( reference, j, 0 ) = sum( this->references[reference].weights * distances( Range::all(), j ) ) / sum( this->references[reference].weights );
	}
}

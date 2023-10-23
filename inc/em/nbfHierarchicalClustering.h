#pragma once

using namespace blitz;

#include <em/nbfClustering.h>
#include <algorithm>

/** Iterative Clustering Method "Loop".
*/
template< class Pixel >
class nbfHierarchicalClustering : public nbfClustering< Pixel >
{
public:

	nbfHierarchicalClustering();

	~nbfHierarchicalClustering(){};

	void setMaxIterations( int max ){ this->maxIterations = max; }

	void setMaxDistance( double d ){ this->maxDistance = d;};
	void setMaxClusters( int d ){ this->maxClusters = d;};
	void setMinOverlap( double o ){ this->minOverlap = o;};

	void setMinElementNumber( double m ){ this->minElementNumber = m;};

	Array< Pixel, 2 > currentTree;
	vector< vector< int > > cumulativeClasses;
	vector< TinyVector< int, 4 > > cumulativeIndexes;
	vector< int > classIndexes;

protected:

	void doClustering( Array< Pixel, 3 > & );

	Pixel classDistance( vector< int > &, vector< int > & );

	Pixel classOverlap( vector< int > &, vector< int > & );

	Array< Pixel, 2 > correlationVector;
	Pixel distanceRadius;

	int maxIterations;

	double maxDistance;
	double minOverlap;
	double minElementNumber;

	bool convergence;
	int maxClusters;
};

template< class Pixel >
nbfHierarchicalClustering< Pixel > :: nbfHierarchicalClustering()
: nbfClustering< Pixel >()
{
	this->convergence = false;
	this->maxIterations = 5;
	
	this->minElementNumber = 1;
	this->maxDistance = numeric_limits< double > :: max();
	this->maxClusters = numeric_limits< int > :: max();
	this->minOverlap = 0.5;
}

template< class Pixel >
void nbfHierarchicalClustering< Pixel > :: doClustering( Array< Pixel, 3 > & classification )
{	
	int N = this->alignments.rows();

	vector< vector< int > > finalClasses;

	vector< vector< int > > classesAux( N );

	// Each image is a class
	for ( int i = 0; i < this->alignments.rows(); i++ ){
		vector< int > aux;
		aux.push_back( i );
		classesAux[ i ] = aux;
	}

	vector< int > globalIndexes;
	for ( int i = 0; i < N; i++ ){
		globalIndexes.push_back(i);
	}

	this->cumulativeClasses.clear();
	this->cumulativeIndexes.clear();

	this->classIndexes.clear();

	bool doneWithClasses = false;

	// if classification tree already computed, only need to build classes
	if ( this->currentTree.size() > 0 ){
		int index = 0;
		while ( index < N - 1 ){

			//// determine number of classes with at least the minimum number of elements
			//int bigEnoughClasses = 0;
			//for ( int h = 0; h < classesAux.size(); h++ ){
			//	if ( classesAux[h].size() >= this->minElementNumber ){
			//		bigEnoughClasses++;
			//	}
			//}

			//bool stopping;
			//if ( this->minElementNumber > 0 ){
			//	stopping = bigEnoughClasses == this->maxClusters;
			//} else {
			//	stopping = classesAux.size() == this->maxClusters;
			//}

			// check if we are done
			if ( ( this->currentTree( index, 2 ) > this->maxDistance ) || ( classesAux.size() == this->maxClusters ) ){
			//if ( ( this->currentTree( index, 2 ) > this->maxDistance ) || ( bigEnoughClasses == this->maxClusters ) ){
			//if ( ( this->currentTree( index, 2 ) > this->maxDistance ) || ( stopping == true ) ){
				if ( doneWithClasses == false ){
					finalClasses.clear();
					finalClasses = classesAux;
					doneWithClasses = true;

					// store global class indexes
					this->classIndexes = globalIndexes;
				}
			}

			// merge second class with first one
			vector< int > :: iterator next = find( globalIndexes.begin(), globalIndexes.end(), this->currentTree( index, 0 ) - 1 );
			int indexI = next - globalIndexes.begin();
			next = find( globalIndexes.begin(), globalIndexes.end(), this->currentTree( index, 1 ) - 1 );
			int indexJ = next - globalIndexes.begin();
			for ( int k = 0; k < classesAux[ indexJ ].size(); k++ ){ 
				classesAux[ indexI ].push_back( classesAux[ indexJ ][k] );
			}
			// eliminate second class
			classesAux.erase( classesAux.begin() + indexJ );
			// update class index
			globalIndexes[ indexI ] = N + index;
			globalIndexes.erase( globalIndexes.begin() + indexJ );

			// if threshold reached, build top-portion of tree
			if ( doneWithClasses == true ){
				this->cumulativeClasses.push_back( classesAux[indexI] );
				// class number, by combining indexI and indexJ
				TinyVector< int, 4 > t( globalIndexes[ indexI ], this->currentTree( index, 0 ) - 1, this->currentTree( index, 1 ) - 1, this->currentTree( index, 2 ) );
				this->cumulativeIndexes.push_back( t );
			}
			index++;
		}
	} else { // otherwise, fully compute tree and classes

		this->currentTree.resize( N - 1, 3 );

		// D stores distances and overlaps
		Array< Pixel, 3 > D( this->alignments.rows(), this->alignments.cols() , 2 );
		Array< Pixel, 3 > Daux;

		D = this->alignments( Range::all(), Range::all(), Range::all(), 0 ); 

		D = D * D / 2;

		// Put \infty in diagonal elements - Distance and Overlap matrix
		for ( int i = 0; i < D.rows(); i++ ){
			D( i, i, Range::all() ) = numeric_limits< Pixel > :: max();
		}

		int indexI, indexJ;

		int aggregateClasses = N;

		for ( int i = 0; i < N - 1; i++ ){

			// find the position of the minimum distance
			Array< Pixel, 2 > subD( D( Range::all(), Range::all(), 0 ) );
			firstIndex h; secondIndex k;
			//indexes = minIndex( subD(k,h) , k );

			TinyVector< int, 2 > new_min_ind;

			if ( this->minOverlap > 0 ){

				Array< Pixel, 1 > values( D.cols() );
				Array< int, 1 > indexes( D.cols() );
				TinyVector< int, 1 > min_ind;

				// Check both conditions
				for ( int j = 0; j < values.size(); j++ ){

					indexes = minIndex( subD(k,h) , k );
					values = min( subD(k,h) , k );

					min_ind = minIndex( values );

					if ( D( min_ind(0) , indexes( min_ind(0) ) , 1 ) < this->minOverlap ){
						//cout << "overlap en iteracion " << i << " es, descartado " << D( min_ind(0) , indexes( min_ind(0) ) , 1 ) <<  endl;
						subD( min_ind(0) , indexes( min_ind(0) ) ) = numeric_limits< Pixel > :: max();
						subD( indexes( min_ind(0)) , min_ind(0)  ) = numeric_limits< Pixel > :: max();
						values( min_ind(0) ) = numeric_limits< Pixel > :: max();
					} else {
						//cout << "overlap en iteracion " << i << " es " << D( min_ind(0) , indexes( min_ind(0) ) , 1 ) <<  endl;
						break;
					}
				}
			} else {
				new_min_ind = minIndex( subD );
			}

			if ( new_min_ind(0) < new_min_ind(1) ){
				indexI = new_min_ind(0);
				indexJ = new_min_ind(1);
			} else {
				indexI = new_min_ind(1);
				indexJ = new_min_ind(0);
			}

			// determine number of classes with at least the minimum number of elements
			int bigEnoughClasses = 0;
			for ( int h = 0; h < classesAux.size(); h++ ){
				if ( classesAux[h].size() >= this->minElementNumber ){
					bigEnoughClasses++;
				}
			}

			// check if we are done
			if ( ( D( indexI, indexJ, 0 ) > this->maxDistance ) || ( classesAux.size() == this->maxClusters ) ){
			//if ( ( D( indexI, indexJ, 0 ) > this->maxDistance ) || ( bigEnoughClasses == this->maxClusters ) ){
				if ( doneWithClasses == false ){
					finalClasses.clear();
					finalClasses = classesAux;
					doneWithClasses = true;

					// store global class indexes
					this->classIndexes = globalIndexes;
				}
			}

			this->currentTree( i, 0 ) = globalIndexes[ indexI ] + 1;
			this->currentTree( i, 1 ) = globalIndexes[ indexJ ] + 1;
			this->currentTree( i, 2 ) = D( indexI, indexJ, 0 );

			vector< vector< int > > saveclasses;
			saveclasses = classesAux;

			// add indexJ part
			for ( int j = 0; j < classesAux[indexJ].size(); j++ ){
				//aux.push_back( classesAux[indexJ][j] );
				classesAux[indexI].push_back( classesAux[indexJ][j] );
			}

			// erase class component J
			classesAux.erase( classesAux.begin() + indexJ );

			// update global cluster indexes
			globalIndexes[ indexI ] = aggregateClasses;
			aggregateClasses++;
			globalIndexes.erase( globalIndexes.begin() + indexJ );

			// if threshold reached, build top-portion of tree
			if ( doneWithClasses == true ){
				this->cumulativeClasses.push_back( classesAux[indexI] );
				// class number, by combining indexI and indexJ
				TinyVector< int, 4 > t( globalIndexes[ indexI ], this->currentTree( i, 0 ) - 1, this->currentTree( i, 1 ) - 1, this->currentTree( i, 2 ) );
				this->cumulativeIndexes.push_back( t );
			}

			// Actualize the distance matrix
			Daux.resize( D.shape() );
			Daux = D;

			D.resize( N - ( i + 1 ), N - ( i + 1 ), 2 );

			// eliminate row & column indexJ from distance matrix
			if ( indexJ == 0 ){
				D = Daux( Range( 1, Daux.ubound(firstDim) ), Range( 1, Daux.ubound(secondDim) ) , Range::all() );
			} else if ( indexJ == Daux.ubound(firstDim) ) {
				D = Daux( Range( fromStart, Daux.ubound(firstDim) - 1 ), Range( fromStart, Daux.ubound(secondDim) - 1 ) , Range::all() );
			} else {
				Range I( 0, indexJ - 1 );
				D( I, I, Range::all() ) = Daux( I, I, Range::all() ); 
				Range J( indexJ, D.ubound(firstDim) );
				Range Jaux( indexJ + 1, Daux.ubound(firstDim) );
				D( I, J, Range::all() ) = Daux( I, Jaux, Range::all() );
				D( J, I, Range::all() ) = Daux( Jaux, I, Range::all() );
				D( J, J, Range::all() ) = Daux( Jaux, Jaux, Range::all() );
			}

			// Calculate the new distances and overlaps
			for ( int j = 0; j < D.rows(); j++ ){
				int jaux = j;
				if ( j >= indexJ ){
					jaux += 1;
				}
				if ( j != indexI ){
					// Distance - complete linkage
					// D( indexI, j, 0 ) = this->classDistance( classesAux[indexI], classesAux[j] );
					// Ward linkage
					D( indexI, j, 0 ) = ( ( saveclasses[jaux].size() + saveclasses[indexI].size() ) * Daux( indexI, jaux, 0 ) + 
						                  ( saveclasses[jaux].size() + saveclasses[indexJ].size() ) * Daux( indexJ, jaux, 0 ) - 
										    saveclasses[jaux].size() * Daux( indexI, indexJ, 0 ) ) / 
						                  ( saveclasses[indexI].size() + saveclasses[indexJ].size() + saveclasses[jaux].size() );

					// Overlap - single linkage (only if using overlap constrain)
					if ( this->minOverlap > 0 ){
						D( indexI, j, 1 ) = this->classOverlap( classesAux[indexI], classesAux[j] );
					}
					// Impose symetry
					D( j, indexI, Range::all() ) = D( indexI, j, Range::all() );
				} else {
					D( j, j, Range::all() ) = numeric_limits< Pixel > :: max();
				}
			}
		}
	}

	this->currentTree( Range::all(), 2 ) = sqrt( this->currentTree( Range::all(), 2 ) );

	if ( finalClasses.size() == 0 ){
		finalClasses = classesAux;
	} else {
		classesAux = finalClasses;
	}

	int j = 0;

	cout << "Total number of classes (before size elimination) = " << finalClasses.size() << endl;

	// check if classes are of sufficient size
	for ( int i = 0; i < classesAux.size(); i++ ){
		if ( classesAux[i].size() < this->minElementNumber ){
			finalClasses.erase( finalClasses.begin() + j );
			this->classIndexes.erase( classIndexes.begin() + j );
		} else {
			j++;
		}
	}

	if ( finalClasses.size() == 0 )
		cout << "No classes were found with the selected threshold and minimum number of classes." << endl;
	
	// Display classes and store in matrix format
	classification.resize( finalClasses.size(), N, 1 );
	classification = 0;

	if ( this->classIndexes.size() == 0 ){
		this->classIndexes.push_back(N);
	}

	// transfer classes to Blitz structure
	for ( int i = 0; i < finalClasses.size(); i++ ){
		for ( int j = 0; j < finalClasses[i].size(); j++ ){
			classification( i, finalClasses[i][j], 0 ) = 1;
		}
	}

	// Post-Processor
	if ( false && ( finalClasses.size() > 0 ) ){
		// evaluate current score
		Array< Pixel, 1 > scores( finalClasses.size() );
		for ( int i = 0; i < finalClasses.size(); i++ ){
			scores(i) = 0;
			for ( int j = 0; j < finalClasses[i].size(); j++ ){
				for ( int k = 0; k < finalClasses[i].size(); k++ ){
					scores(i) += pow2( this->alignments( finalClasses[i][j], finalClasses[i][k], 0 ) );
				}
				//Pixel csum = 0;
				//for ( int k = 0; k < finalClasses[i].size(); k++ ){
				//	csum += this->alignments( finalClasses[i][j], finalClasses[i][k], 0 );
				//}
				//scores(i) += csum * csum / finalClasses[i].size() / finalClasses[i].size();
				//scores(i) += csum / finalClasses[i].size();
			}
			//scores(i) = scores(i) / 2.0 / finalClasses[i].size();
			scores(i) = scores(i) / 2.0;
		}
		cout << "Starting classification score " << sum(scores) << endl;

		// do iterations until convergence (or until maxIters reached)
		bool notFinished = true;
		int maxIters = 1000;
		int currentIters = 0;
		while ( ( notFinished == true ) && ( currentIters < maxIters ) ){

			currentIters++;
			notFinished = false;

			// for all volumes
			for ( int i = 0; i < this->alignments.rows(); i++ ){
				Pixel currentScore = sum( scores );
				TinyVector< int, 1 > cpc = maxIndex( classification( Range::all(), i, 0 ) );
				int currentPointClass = cpc[0];
				if ( classification( currentPointClass, i, 0 ) == 0 ){
					continue;
				}
				Pixel bestScore = currentScore;
				int bestClass = currentPointClass;
				Pixel bestScoreOffset = 0;
				Pixel bestNewCurrentScoreOffset = 0;
				// compute change in score if we remove point i from current assigned class
				Pixel newCurrentClassScoreOffset = 0;
				for ( int c = 0; c < finalClasses[currentPointClass].size(); c++ ){
					newCurrentClassScoreOffset += pow2( this->alignments( i, finalClasses[currentPointClass][c], 0 ) );
				}
				//// compute ESS if we remove point i from current assigned class
				//for ( int j = 0; j < finalClasses[currentPointClass].size(); j++ ){
				//	if ( finalClasses[currentPointClass][j] != i ){
				//		Pixel csum = 0;
				//		for ( int k = 0; k < finalClasses[currentPointClass].size(); k++ ){
				//			if ( finalClasses[currentPointClass][k] != i ){
				//				csum += this->alignments( finalClasses[currentPointClass][j], finalClasses[currentPointClass][k], 0 );
				//			}
				//		}
				//		newCurrentClassScoreOffset += csum * csum / (finalClasses[currentPointClass].size()-1) / (finalClasses[currentPointClass].size()-1);
				//		//newCurrentClassScoreOffset += csum / (finalClasses[currentPointClass].size()-1);
				//	}
				//}

				for ( int j = 0; j < finalClasses.size(); j++ ){
					// try shifting volume i to class j (if not already a member)
					if ( j != currentPointClass ){
						Pixel newOtherClassScoreOffset = 0;
						for ( int c = 0; c < finalClasses[j].size(); c++ ){
							newOtherClassScoreOffset += pow2( this->alignments( i, finalClasses[j][c], 0 ) );
						}
						//// compute ESS if we remove point i from current assigned class
						//Pixel csum;
						//for ( int c = 0; c < finalClasses[j].size(); c++ ){
						//	csum = 0;
						//	for ( int k = 0; k < finalClasses[j].size(); k++ ){
						//		csum += this->alignments( finalClasses[j][c], finalClasses[j][k], 0 );
						//	}
						//	csum += this->alignments( finalClasses[j][c], i, 0 );
						//	newOtherClassScoreOffset += csum * csum / (finalClasses[j].size()+1) / (finalClasses[j].size()+1);
						//	//newOtherClassScoreOffset += csum / (finalClasses[j].size()+1);
						//}
						//csum = 0;
						//for ( int k = 0; k < finalClasses[j].size(); k++ ){
						//	csum += this->alignments( i, finalClasses[j][k], 0 );
						//}
						//newOtherClassScoreOffset += csum * csum / (finalClasses[j].size()+1) / (finalClasses[j].size()+1);
						////newOtherClassScoreOffset += csum / (finalClasses[j].size()+1);

						Pixel tentativeScore = currentScore - newCurrentClassScoreOffset + newOtherClassScoreOffset;
						// Pixel tentativeScore = currentScore - scores(currentPointClass) - scores(j) + newCurrentClassScoreOffset + newOtherClassScoreOffset;
						//if ( finalClasses[currentPointClass].size() > 1 ){
						//	Pixel tentativeScore = currentScore - scores(currentPointClass) - scores(j) + 
						//		( scores(currentPointClass) * finalClasses[currentPointClass].size() - newCurrentClassScoreOffset ) / ( finalClasses[currentPointClass].size() - 1 ) + 
						//		( scores(j) * finalClasses[j].size() + newOtherClassScoreOffset ) / ( finalClasses[j].size() + 1 );
							if ( tentativeScore < bestScore ){
								bestScore = tentativeScore;
								bestClass = j;
								bestScoreOffset = newOtherClassScoreOffset;
								bestNewCurrentScoreOffset = newCurrentClassScoreOffset;
							}
						//}
					}
				}
				// make the change permanent
				if ( bestClass != currentPointClass ){

					//scores(currentPointClass) = ( scores(currentPointClass) * finalClasses[currentPointClass].size() - bestNewCurrentScoreOffset ) / ( finalClasses[currentPointClass].size() - 1 );
					//scores(bestClass ) = ( scores(bestClass) * finalClasses[bestClass].size() + bestScoreOffset ) / ( finalClasses[bestClass].size() + 1 );
					// remove from old class
					vector< int > :: iterator result = find( finalClasses[currentPointClass].begin(), finalClasses[currentPointClass].end(), i );
					finalClasses[currentPointClass].erase( result );
					classification( currentPointClass, i, 0 ) = 0;
					// add to new class
					finalClasses[bestClass].push_back( i );
					classification( bestClass, i, 0 ) = 1;
					scores( currentPointClass ) -= bestNewCurrentScoreOffset;
					scores( bestClass ) += bestScoreOffset;
					notFinished = true;
				}
			}
#ifdef WIN32
			// recompute scores to see the improvement
			cout << "Current score = " << sum(scores) << endl;
#endif
		}
		// recompute scores to see the improvement
		cout << "Score at iteration " << currentIters << " = " << sum(scores) << endl;
	}

	for ( int i = 0; i < finalClasses.size(); i++ ){
		cout << "Class " << this->classIndexes[i] << "(" << sum( classification( i, Range::all(), 0 ) ) << "): " ;
		for ( int j = 0; j < classification.cols(); j++ ){
			if ( classification( i, j ) > 0 ){
				cout << j << " " ;
			}
		}
		cout << endl;
	}

}

template< class Pixel >
Pixel nbfHierarchicalClustering< Pixel > :: classDistance( vector< int > & class1, vector< int > & class2 )
{
	Pixel dist = 0;

	// complete linkage
	for ( int i = 0; i < class1.size(); i++ ){		
		for ( int j = 0; j < class2.size(); j++ ){
			dist = max( dist, this->alignments( class1[i], class2[j], 0, 0 ) );
		}
	}
	
	//// mean linkage
	//Pixel count = 0;
	//for ( int i = 0; i < class1.size(); i++ ){		
	//	for ( int j = 0; j < class2.size(); j++ ){
	//		dist += this->alignments( class1[i], class2[j], 0, 0 );
	//		count++;
	//	}
	//}
	//return dist / count;

	return dist;
}

template< class Pixel >
Pixel nbfHierarchicalClustering< Pixel > :: classOverlap( vector< int > & class1, vector< int > & class2 )
{	
	// exact same code as above... may change

	Pixel dist = 0;

	for ( int i = 0; i < class1.size(); i++ ){
		for( int j = 0; j < class2.size(); j++ ){
			dist = max( dist, this->alignments( class1[i], class2[j], 2, 0 ) );
		}
	}
	
	//cout << class1[0] << " " << class2[0] << endl;

	return dist;
}


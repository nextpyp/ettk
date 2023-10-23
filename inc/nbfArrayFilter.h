#ifndef FILE_nbfArrayFilter
#define FILE_nbfArrayFilter

#include <bs/nbfBordStrategyConst.h>

template< class Pixel, int const Dim >
class nbfArrayFilter
{
public:

	// constructor takes weight array as input
	nbfArrayFilter( Array< Pixel, Dim > & );

	virtual ~nbfArrayFilter(){};

	void execute( Array< Pixel, Dim > & ){};

	// Matlab's graythresh (assumes 0 < input < 1 )
	static Pixel graythresh( Array< Pixel, Dim > & );

protected:

	// filter I/O
	Array< Pixel, Dim > * input;
	Array< Pixel, Dim > * output;

};

template< class Pixel, int const Dim >
nbfArrayFilter< Pixel, Dim > :: nbfArrayFilter( Array< Pixel, Dim > & input )
{
	this->input = &input;
}

template< class Pixel, int const Dim >
Pixel nbfArrayFilter< Pixel, Dim > :: graythresh( Array< Pixel, Dim > & input )
{
	// Compute optimal threshold using matlab's graythresh(). 
	// The computation can be restricted to an optional mask.
	// Mask values are taken from the Red component of mask.

    int i;
    float n = 0; // number of points
	
	int num_bins = 256;
	double counts[256];
	for ( i = 0; i < num_bins; i++ ){
		counts[i] = 0;
	}

	n = input.rows() * input.cols() * input.depth();

	Pixel imax = max(input);
	Pixel imin = min(input);

	// counts = imhist(I,num_bins);
	// p = counts / sum(counts);
	Array< Pixel, Dim > :: const_iterator iter = input.begin();
	while( iter != input.end() ){
		int t = vtkMath::Round( (*iter) * ( num_bins - 1 ) );
		counts[ t ] += ( 1.0 / n );
		++iter;
	}

	// omega = cumsum(p);
	// mu = cumsum(p .* (1:num_bins)');
	double omega[256];
	double mu[256];
	omega[0] = counts[0];
	mu[0] = counts[0];
	for( i = 1; i < num_bins; i ++ ) {
		omega[i] = omega[i-1] + counts[i];
		mu[i] = mu[i-1] + counts[i] * ( i + 1 );
	}

	// mu_t = mu(end);
	double mu_t = mu[num_bins-1];

	// sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));
	// maxval = max(sigma_b_squared);
	double sigma_b_squared[256];
	double maxval = 0;
	for( i = 0; i < num_bins; i ++ ) {
		if ( omega[i] > 0 ){
			sigma_b_squared[i] = pow( mu_t * omega[i] - mu[i], 2 ) / ( omega[i] * ( 1 - omega[i] ) );
			if ( sigma_b_squared[i] > maxval ){
				maxval = sigma_b_squared[i];
			}
		}
		else{
			sigma_b_squared[i] = 0;
		}
	}

	// idx = mean(find(sigma_b_squared == maxval));
	double idx = 0;
	int idx_count = 0;
	for( i = 0; i < num_bins; i ++ ) {
		if ( sigma_b_squared[i] == maxval ){
			idx += ( i + 1 );
			idx_count += 1;
		}
	}
	idx /= idx_count;

	// level = (idx - 1) / (num_bins - 1);
    return ( idx - 1 ) / ( num_bins - 1 );
}

#endif /* FILE_nbfArrayFilter */

#pragma once

/** Data structure to handle individual pairs y_j = ( a_ij, x_j ).
*/
template< class Pixel >
class nbfRadonElement
{
public:

	// Constructor takes reference to input/output pixels
	nbfRadonElement( Pixel *, Pixel * );

	~nbfRadonElement(){};

	// Calculate pixel contribution: x * w
	inline Pixel project(){
		return ( (*this->pInput) * this->w );
	}

	// Reverse pixel contribution: v * w
	inline void backProject( Pixel v ){
		(*this->pOutput) += ( v * this->w );
	}

	// Cumulative pixel contributions (used when constructing matrix A)
	void incrementWeight( Pixel w ){
		this->w += w;
	}

	// store pointer to input data
	Pixel * pInput;

	// store pointer to output data
	Pixel * pOutput;

	// store component contribution
	Pixel w;

};

template< class Pixel >
nbfRadonElement< Pixel > :: nbfRadonElement( Pixel * in, Pixel * out )
: pInput(in), pOutput(out), w(0.0)
{
}
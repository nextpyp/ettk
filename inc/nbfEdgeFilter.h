#ifndef FILE_nbfRadon
#define FILE_nbfRadon

#include <vtkMath.h>
#include <nbfLinearInterpolator.h>

template< class Pixel, int const Dim >
class nbfRadon : public nbfArrayFilter< Pixel, Dim >
{
public:

	// constructor takes weight array as input
	nbfRadon( Array< Pixel, Dim > & );

	~nbfRadon(){};

	void setInput( Array< Pixel, 2 > & a ){ this->image.reference(a); };
	void setAngles( Array< Pixel, 1 > & a ){ this->angles.reference(a); };

	// store result in argument
	void execute( Array< Pixel, 2 > & );

protected:

	Array< Pixel, 2 > radon;
	Array< Pixel, 2 > image;
	Array< Pixel, 1 > angles;

};

template< class Pixel, int const Dim >
nbfRadon< Pixel, Dim > :: nbfRadon( Array< Pixel, Dim > & input )
: nbfArrayFilter< Pixel, Dim >( input )
{
}

template< class Pixel, int const Dim >
void nbfRadon< Pixel, Dim > :: forwardProjection( Array< Pixel, 2 > & radon )
{
	Pixel size = this->image.rows();
	this->radon.resize( size, angles.rows() );

	firstIndex i;
	secondIndex j;
	Array< bool, 2 > mask( this->image.shape() );
	Pixel maxR = size * size / 4.0;
	mask = pow2( i - size /2.0 ) + pow( j - size / 2.0) < maxR; 

	Array< bool, 2 > :: const_iterator iMask = mask.begin();
	Array< Pixel, 2 > :: iterator iImage = this->image.begin();
	while( iMask != mask.end() ){
		if (*iMask){
			for( int index = 0; index < this->angles.size(); i++ ){
				Pixel distance = ( TinyVector< Pixel, 2 >( tan(*iAngles), 1.0 ) * iImage->position() ) * cos(*iAngles);
				int fDistance = floor(distance);
				radon( i, fDistance )     += (  fDistance + 1.0 - distance ) * (*iImage);
				radon( i, fDistance + 1 ) += (        distance - fDistance ) * (*iImage);
			}
		}
		++iMask;
		++iImage;
	}
}


template< class Pixel, int const Dim >
void nbfRadon< Pixel, Dim > :: backwardProjection( Array< Pixel, 2 > & angle )
{
	Array< bool, 2 > :: const_iterator iMask = this->mask.begin();
	Array< Pixel, 2 > :: iterator iImage = this->image.begin();
	while( iMask != this->mask.end() ){
		if (*iMask){
			Pixel distance = TinyVector< Pixel, 2 >( tan(angle), 1.0 ) * iImage->position();
			int fDistance = floor(distance);
			(*iImage) = ( fDistance + 1.0 - distance ) * projection( fDistance ) + ( distance - fDistance ) * projection( fDistance + 1 );
		}
		++iMask;
		++iImage;
	}
}

#endif /* FILE_nbfRadon */
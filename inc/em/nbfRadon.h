#ifndef FILE_nbfRadon
#define FILE_nbfRadon

////////////////////////!!!!!!!!!
// REPLACE EDGEFILTER //!!!!!!!!!
////////////////////////!!!!!!!!!

template< class Pixel >
class nbfRadon
{
public:

	// constructor takes weight array as input
	nbfRadon(){};

	~nbfRadon(){};

	void setImage( Array< Pixel, 2 > & );
	void setAngles( Array< Pixel, 1 > & a ){ this->angles.reference(a); };

	// store result in argument
	void forwardProjection(int);
	void fastForwardProjection(int);
	void forwardProjection();
	void backwardProjection(Array< Pixel, 1 > &,int,Array<Pixel,2>&);

	void incrementRadon( Pixel, Pixel );

	void art( Array< Pixel, 2 > &, Array< Pixel, 2 > &, int );

	void getRadon( Array< Pixel, 2 > & a ){ a.reference( this->radon ); };
	void getImage( Array< Pixel, 2 > & a ){ a.reference( this->image ); };

protected:

	Array< Pixel, 2 > radon;
	Array< Pixel, 2 > image;
	Array< bool, 2 > mask;
	Array< Pixel, 1 > angles;
	int size;
	Pixel center;

	Array< Pixel, 1 > singleProjection;
	Pixel * singleProjectionPointer;

	Array< Pixel, 1 > sinPos;
	Array< Pixel, 1 > sinNeg;
	Array< Pixel, 1 > cosPos;
	Array< Pixel, 1 > cosNeg;


};

template< class Pixel >
void nbfRadon< Pixel > :: setImage( Array< Pixel, 2 > & image )
{
	// assign to attribute
	this->image.reference( image );

	this->size = image.rows();
	this->center = ( this->size - 1 ) / 2.0;

	// build circular mask
	firstIndex i;
	secondIndex j;
	this->mask.resize( this->image.shape() );

	// assume square image!
	Pixel maxR = pow2( image.rows() - 2 ) / 4.0;

	// true inside, false outside uniit circle
	Pixel size = ( this->image.rows() - 1 ) / 2.0;
	this->mask = pow2( i - size ) + pow2( j - size ) < maxR; 

	this->sinPos.resize( this->size );
	this->sinNeg.resize( this->size );
	this->cosPos.resize( this->size );
	this->cosNeg.resize( this->size );

	this->singleProjection.resize( this->size );
	this->singleProjectionPointer = this->singleProjection.dataZero();

	nbfMatlabWriter w;
	w.setFileName("test.blitz");
	w.write(this->mask);
}


template< class Pixel >
void nbfRadon< Pixel > :: forwardProjection()
{
	this->radon.resize( this->size, this->angles.rows() );
	this->radon = 0;
	return;
	for ( int i = 0; i < this->angles.numElements(); i++ ){
		this->forwardProjection(i);
		this->radon( Range::all(), i ) = this->singleProjection;
	}
}


template< class Pixel >
void nbfRadon< Pixel > :: art( Array< Pixel, 2 > & targetRadon, Array< Pixel, 2 > & Phi, int numIterations )
{
	// build normalizing vector
	Array< Pixel, 1 > normalize( this->size );
	secondIndex jIndex;
	normalize = sum( this->mask, jIndex );

	this->image = Phi * this->mask;

 	nbfMatlabWriter w;
	w.setFileName("test.blitz");
	//w.write(this->image);

	Array< Pixel, 1 > correction( this->size );
 	for ( int iters = 0; iters < numIterations; iters++ ){
		for ( int thetaIndex = 0; thetaIndex < this->angles.numElements(); thetaIndex++ ){
			this->forwardProjection(thetaIndex);
 			correction = targetRadon( Range::all(), thetaIndex ) - this->singleProjection / normalize / 10.0;
			correction( 0 ) = 0;
			correction( correction.ubound(0) ) = 0;
			this->backwardProjection( correction, thetaIndex, Phi );
			this->image += Phi;
		}
	}
	Phi.reference( this->image );
}

template< class Pixel >
inline void nbfRadon< Pixel > :: incrementRadon( Pixel distance, Pixel imageValue )
{
	int fDistance = distance;	
	Pixel delta = distance - fDistance;
	//this->singleProjection(fDistance) += ( 1.0 - delta ) * imageValue;
	//this->singleProjection(fDistance+1) += delta * imageValue;
	this->singleProjectionPointer[ fDistance++ ] += ( 1.0 - delta ) * imageValue;
	this->singleProjectionPointer[ fDistance ] += delta * imageValue;
}

template< class Pixel >
void nbfRadon< Pixel > :: forwardProjection( int angleIndex )
{
	//cout << angleIndex << endl;
	Pixel * imagePointer = this->image.dataZero();
	bool * maskPointer = this->mask.dataZero();

	int numberOfAngles = this->angles.numElements();

	Pixel distance, angle, sine, cosine, imageValue;

	angle = this->angles(angleIndex);

	Pixel sinposval, sinnegval, cosposval, cosnegval;

	sine = sin( angle );
	cosine = cos( angle );

	//firstIndex fi;
	//sinPos = ( sine * ( fi - this->center + .25 ) ) + this->center;
	//sinNeg = sinPos - .5 * sine;
	//cosPos = ( cosine * ( fi - this->center + .25 ) );
	//cosNeg = cosPos - .5 * cosine;

	Pixel * sinPosPointer = sinPos.dataZero();
	Pixel * sinNegPointer = sinNeg.dataZero();
	Pixel * cosPosPointer = cosPos.dataZero();
	Pixel * cosNegPointer = cosNeg.dataZero();
	for ( int fi = 0; fi < this->size; fi++ ){
		sinPosPointer[fi] = sine * ( fi - this->center + .25 ) + this->center;
		sinNegPointer[fi] = sinPosPointer[fi] - .5 * sine;
		cosPosPointer[fi] = cosine * ( fi - this->center + .25 );
		cosNegPointer[fi] = cosPosPointer[fi] - .5 * cosine;
	}


	this->singleProjection = 0;

	for ( int i = 0; i < this->size; i++ ){
		sinposval = sinPos(i);
		sinnegval = sinNeg(i);
		for ( int j = 0; j < this->size; j++ ){
			int current = i * this->size + j;
			if ( maskPointer[current] ){
				cosposval = cosPos(j);
				cosnegval = cosNeg(j);
				distance = sinposval + cosposval;

				// retrieve image value and weight contribution by 1/4
				imageValue = imagePointer[ current ] * .25;
				
				//imageValue = this->image(i,j);

				this->incrementRadon(distance,imageValue);

				distance = sinposval + cosnegval;
				this->incrementRadon(distance,imageValue);

				distance = sinnegval + cosposval;
				this->incrementRadon(distance,imageValue);

				distance = sinnegval + cosnegval;
				this->incrementRadon(distance,imageValue);
			}
		}
	}
}

template< class Pixel >
void nbfRadon< Pixel > :: fastForwardProjection( int angleIndex )
{
	//cout << angleIndex << endl;
	Pixel * imagePointer = this->image.dataZero();
	bool * maskPointer = this->mask.dataZero();

	int numberOfAngles = this->angles.numElements();

	Pixel distance, angle, sine, cosine, imageValue;

	angle = this->angles(angleIndex);

	Pixel sinposval;

	sine = sin( angle );
	cosine = cos( angle );

	//firstIndex fi;
	//sinPos = ( sine * ( fi - this->center ) ) + this->center;
	//cosPos = ( cosine * ( fi - this->center ) );

	Pixel * sinPosPointer = sinPos.dataZero();
	Pixel * cosPosPointer = cosPos.dataZero();
	for ( int fi = 0; fi < this->size; fi++ ){
		sinPosPointer[fi] = sine * ( fi - this->center ) + this->center;
		cosPosPointer[fi] = cosine * ( fi - this->center );
	}

	this->singleProjection = 0;

	for ( int i = 0; i < this->size; i++ ){
		sinposval = sinPos(i);
		for ( int j = 0; j < this->size; j++ ){
			int current = i * this->size + j;
			if ( maskPointer[current] ){
				distance = sinposval + cosPos(j);

				//imageValue = this->image(i,j);
				imageValue = imagePointer[ current ];

				this->incrementRadon(distance,imageValue);
			}
		}
	}
}

template< class Pixel >
void nbfRadon< Pixel > :: backwardProjection( Array< Pixel, 1 > & p, int angle,Array< Pixel,2>&BP )
{
	Array< Pixel, 2 > :: iterator iImage = BP.begin();
	Array< bool, 2 > :: iterator iMask = this->mask.begin();
	Pixel * pPointer = p.dataZero();
	Pixel angular = this->angles(angle);
	TinyVector< Pixel, 2 > ang( sin(angular), cos(angular) );
	while( iImage != BP.end() ){
		if ( *iMask ){
			Pixel distance = dot( ang, iImage.position() - this->center ) + this->center;
			//Pixel distance = sin(angular) * ( iImage.position()[0] - this->center ) + cos(angular) * ( iImage.position()[1] - this->center ) + this->center;
			int fDistance = distance;
			//if ( fDistance > distance )
			//	fDistance--;
			Pixel delta = distance - fDistance;
			(*iImage) = ( 1.0 - delta ) * pPointer[ fDistance ] + delta * pPointer[ fDistance + 1 ];
		}
		else{
			(*iImage) = 0;
		}
		++iImage;
		++iMask;
	}
}

#endif /* FILE_nbfRadon */
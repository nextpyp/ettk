#ifndef FILE_BordStartegyMirror
#define FILE_BordStartegyMirror

#include <bs/nbfBordStrategy.h>
#include <vector>

/*******************************
 * BordStrategyMirrorSimpleImp *
 *******************************/

template< class Pixel >
class BordStrategyMirrorSimpleImp
{
protected:

  template< int const Dim >
  void refresh( Array< Pixel, Dim > &, // superImage
		Array< Pixel, Dim > & ){}; // imageSource

  void refresh( Array< Pixel, 2 > &, // superImage
		Array< Pixel, 2 > & ); // imageSource

  void refresh( Array< Pixel, 3 > &, // superImage
		Array< Pixel, 3 > & ); // imageSource

public:
  BordStrategyMirrorSimpleImp();
  virtual ~BordStrategyMirrorSimpleImp();

};


/****************************
 * BordStrategyMirrorSimple *
 ****************************/

/** Estrategia Espejo Simple.
    Esta estrategia rellena los bordes con un criterio de espejo simple, es decir, se hace una simetría axial en torno a la fila del borde correspondiente.
 */
template< class Pixel, int const Dim >
class BordStrategyMirrorSimple : public BordStrategyMirrorSimpleImp< Pixel >, public BordStrategy< Pixel, Dim >
{
 public:

  void refresh();

  unsigned getId(){ return BS_ID_MIRROR_SIMPLE;}

  /** Constructor. Constructor que toma como parametro el valor constante de la estrategia de borde. */
  BordStrategyMirrorSimple( Array< Pixel, Dim > &, unsigned );
  ~BordStrategyMirrorSimple();
};


/*******************************
 * BordStrategyMirrorDoubleImp *
 *******************************/

template< class Pixel >
class BordStrategyMirrorDoubleImp
{
protected:

  template< int const Dim >
  void refresh( Array< Pixel, Dim > &, // superImage
		Array< Pixel, Dim > & ){}; // imageSource

  void refresh( Array< Pixel, 2 > &, // superImage
		Array< Pixel, 2 > & ); // imageSource

  void refresh( Array< Pixel, 3 > &, // superImage
		Array< Pixel, 3 > & ); // imageSource

public:
  BordStrategyMirrorDoubleImp();
  virtual ~BordStrategyMirrorDoubleImp();

};


/****************************
 * BordStrategyMirrorDouble *
 ****************************/

/** Estrategia Espejo Doble.
    Esta estrategia rellena los bordes con un criterio de espejo doble, es decir, se periodiza la imagen en todas las direcciones.
 */
template< class Pixel, int const Dim >
class BordStrategyMirrorDouble : public BordStrategyMirrorDoubleImp< Pixel >, public BordStrategy< Pixel, Dim >
{
 public:

  void refresh();

  unsigned getId(){ return BS_ID_MIRROR_DOUBLE;}

  /** Constructor. Constructor que toma como parametro el valor constante de la estrategia de borde. */
  BordStrategyMirrorDouble( Array< Pixel, Dim > &, unsigned );
  ~BordStrategyMirrorDouble();
};


/******************************************************
 * IMPLEMENTACION
 ******************************************************/

/****************************
 * BordStrategyMirrorSimple *
 ****************************/


template< class Pixel, int const Dim >
BordStrategyMirrorSimple< Pixel, Dim > :: BordStrategyMirrorSimple( Array< Pixel, Dim > & A, 
																   unsigned offset = 2 )
: BordStrategy< Pixel, Dim >( A, offset )
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("BordStrategyMirrorSimple");
#endif
}


template< class Pixel, int const Dim >
BordStrategyMirrorSimple< Pixel, Dim > :: ~BordStrategyMirrorSimple()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  DELETING("BordStrategyMirrorSimple");
#endif
}


template< class Pixel, int const Dim >
void BordStrategyMirrorSimple< Pixel, Dim > :: refresh(){
  //cout << this->superImage << endl;
  BordStrategyMirrorSimpleImp< Pixel > :: refresh( this->superImage, *this->imageSource );
  //cout << this->superImage << endl;
}


/*******************************
 * BordStrategyMirrorSimpleImp *
 *******************************/


template< class Pixel >
BordStrategyMirrorSimpleImp< Pixel > :: BordStrategyMirrorSimpleImp()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("BordStrategyMirrorSimpleImp");
#endif
}


template< class Pixel >
BordStrategyMirrorSimpleImp< Pixel > :: ~BordStrategyMirrorSimpleImp()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  DELETING("BordStrategyMirrorSimpleImp");
#endif
}


template< class Pixel >
void BordStrategyMirrorSimpleImp< Pixel > :: refresh( Array< Pixel, 2 > & superImage,
						      Array< Pixel, 2 > & imageSource )
{
  int rows = imageSource.rows();
  int cols = imageSource.cols();
  int superRows = superImage.rows();
  int superCols = superImage.cols();

  int r = ( superRows - rows ) / 2;

  // left
  Array< Pixel, 2 > left( superRows, r );
  left = superImage( Range::all(), Range(r+1,2*r) );
  left.reverseSelf(secondDim);
  superImage( Range::all(), Range(fromStart,r-1) ) = left;

  // right
  Array< Pixel, 2 > right( superRows, r );
  right = superImage( Range::all(), Range(superCols-1-r-1-(r-1),superCols-1-r-1) );
  right.reverseSelf(secondDim);
  superImage( Range::all(), Range(superCols-r,toEnd) ) = right;
  
  // top
  Array< Pixel, 2 > top( r, superCols );
  top = superImage( Range(r+1,2*r), Range::all() );
  top.reverseSelf(firstDim);
  superImage( Range(fromStart,r-1), Range::all() ) = top;

  // bottom
  Array< Pixel, 2 > bottom( r, superCols );
  bottom = superImage( Range(superRows-1-r-1-(r-1), superRows-1-r-1 ), Range::all() );
  bottom.reverseSelf(firstDim);
  superImage( Range(superRows-r,toEnd), Range::all() ) = bottom;

}


template< class Pixel >
void BordStrategyMirrorSimpleImp< Pixel > :: refresh( Array< Pixel, 3 > & superImage,
						      Array< Pixel, 3 > & imageSource )
{
  int rows = imageSource.rows();
  int cols = imageSource.cols();
  int superRows = superImage.rows();
  int superCols = superImage.cols();
  int depth = imageSource.depth();
  int superDepth = superImage.depth();

  int r = ( superRows - rows ) / 2;

  // left
  Array< Pixel, 3 > left( superRows, superCols, r );
  left = superImage( Range::all(), Range::all(), Range(r+1,2*r) );
  left.reverseSelf(thirdDim);
  superImage( Range::all(), Range::all(), Range(fromStart,r-1) ) = left;

  // right
  Array< Pixel, 3 > right( superRows, superCols, r );
  right = superImage( Range::all(), Range::all(), Range(superDepth-1-r-1-(r-1),superDepth-1-r-1) );
  right.reverseSelf(thirdDim);
  superImage( Range::all(), Range::all(), Range(superDepth-r,toEnd) ) = right;
  
  // top
  Array< Pixel, 3 > top( superRows, r, superDepth );
  top = superImage( Range::all(), Range(r+1,2*r), Range::all() );
  top.reverseSelf(secondDim);
  superImage( Range::all(), Range(fromStart,r-1), Range::all() ) = top;

  // bottom
  Array< Pixel, 3 > bottom( superRows, r, superDepth );
  bottom = superImage( Range::all(), Range(superCols-1-r-1-(r-1), superCols-1-r-1 ), Range::all() );
  bottom.reverseSelf(secondDim);
  superImage( Range::all(), Range(superCols-r,toEnd), Range::all() ) = bottom;

  // front
  Array< Pixel, 3 > front( r, superCols, superDepth );
  front = superImage( Range(r+1,2*r), Range::all(), Range::all() );
  front.reverseSelf(firstDim);
  superImage( Range(fromStart,r-1), Range::all(), Range::all() ) = front;

  // back
  Array< Pixel, 3 > back( r, superCols, superDepth );
  back = superImage( Range(superRows-1-r-1-(r-1), superRows-1-r-1 ), Range::all(), Range::all() );
  back.reverseSelf(firstDim);
  superImage( Range(superRows-r,toEnd), Range::all(), Range::all() ) = back;
}


/****************************
 * BordStrategyMirrorDouble *
 ****************************/

template< class Pixel, int const Dim >
BordStrategyMirrorDouble< Pixel, Dim > :: BordStrategyMirrorDouble( Array< Pixel, Dim > & A,
																   unsigned offset = 2 )
: BordStrategy< Pixel, Dim >( A, offset )
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("BordStrategyMirrorDouble");
#endif
}


template< class Pixel, int const Dim >
BordStrategyMirrorDouble< Pixel, Dim > :: ~BordStrategyMirrorDouble()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  DELETING("BordStrategyMirrorDouble");
#endif
}


template< class Pixel, int const Dim >
void BordStrategyMirrorDouble< Pixel, Dim > :: refresh()
{
  BordStrategyMirrorDoubleImp< Pixel > :: refresh( this->superImage, *this->imageSource );
  // cout << this->superImage << endl;
}


/*******************************
 * BordStrategyMirrorDoubleImp *
 *******************************/


template< class Pixel >
BordStrategyMirrorDoubleImp< Pixel > :: BordStrategyMirrorDoubleImp()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("BordStrategyMirrorDoubleImp");
#endif
}


template< class Pixel >
BordStrategyMirrorDoubleImp< Pixel > :: ~BordStrategyMirrorDoubleImp()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  DELETING("BordStrategyMirrorDoubleImp");
#endif
}


template< class Pixel >
void BordStrategyMirrorDoubleImp< Pixel > :: refresh( Array< Pixel, 2 > & superImage,
						      Array< Pixel, 2 > & imageSource )
{
  int rows = imageSource.rows();
  int cols = imageSource.cols();
  int superRows = superImage.rows();
  int superCols = superImage.cols();

  vector< int > I, J;
  int r = ( superRows - rows ) / 2;

  for ( int i = 0; i < rows; i++){

    // COORDENADA I
    // valor hacia atras
    if ( r - i - 1 >= 0){
      I.push_back( r - i - 1);
    }

    // valor actual
    I.push_back( i + r );

    // valor hacia adelante
    if ( 2 * ( superRows - r) - i - r - 1 < superRows ){
      I.push_back( 2 * ( superRows - r ) - i - r - 1 );
    }

    for ( int j = 0; j < cols; j++){

      // if en la zona de donde saco valores
      if ( !( ( r <= i ) &&
	      ( i < rows - r ) &&
	      ( r <= j ) &&
	      ( j < cols - r )
	      ) ) {


	// COORDENADA J
	// valor hacia atras
	if ( r - j - 1 >= 0 ){
	  J.push_back( r - j - 1 );
	}

	// valor actual
	J.push_back( j + r );

	// valor hacia adelante
	if ( 2 * ( superCols - r ) - j - r - 1 < superCols ){
	  J.push_back( 2 * ( superCols - r) - j - r - 1 );
	}

	superImage[ indexSet( I, J) ] = imageSource(i,j);

	// vacio los indices
	J.clear();
      }

      // sino se esta en la condicion, simplemente copio la chica en la grande.
      else superImage( i + r, j + r) = imageSource(i,j);
    }

    I.clear();
  }
}


template< class Pixel >
void BordStrategyMirrorDoubleImp< Pixel > :: refresh( Array< Pixel, 3 > & superImage,
						      Array< Pixel, 3 > & imageSource )
{
  int rows = imageSource.rows();
  int cols = imageSource.cols();
  int superRows = superImage.rows();
  int superCols = superImage.cols();
  int depth = imageSource.depth();
  int superDepth = superImage.depth();

  vector< int > I, J, K;
  int r = ( superRows - rows ) / 2;

  for ( int i = 0; i < rows; i++){

    // COORDENADA I
    // valor hacia atras
    if ( r - i - 1 >= 0){
      I.push_back( r - i - 1 );
    }

    // valor actual
    I.push_back( i + r );

    // valor hacia adelante
    if ( 2 * ( superRows - r) - i - r - 1 < superRows ){
      I.push_back( 2 * ( superRows - r ) - i - r - 1 );
    }

    for ( int j = 0; j < cols; j++){

      // COORDENADA J
      // valor hacia atras
      if ( r - j - 1 >= 0 ){
	J.push_back( r - j - 1 );
      }

      // valor actual
      J.push_back( j + r );

      // valor hacia adelante
      if ( 2 * ( superCols - r ) - j - r - 1 < superCols ){
	J.push_back( 2 * ( superCols - r) - j - r - 1 );
      }

      for ( int k = 0; k < depth; k++){

	// if en la zona de donde saco valores
	if ( !( ( r <= i ) &&
		( i < rows - r ) &&
		( r <= j ) &&
		( j < cols - r ) &&
		( r <= k ) &&
		( k < depth - r )
		) ) {

	  // COORDENADA K
	  // valor hacia atras
	  if ( r - k - 1 >= 0 ){
	    K.push_back( r - k - 1 );
	  }

	  // valor actual
	  K.push_back( k + r );

	  // valor hacia adelante
	  if ( 2 * ( superDepth - r ) - k - r - 1 < superDepth ){
	    K.push_back( 2 * ( superDepth - r) - k - r - 1 );
	  }

	  superImage[ indexSet( I, J, K) ] = imageSource(i,j,k);

	  // vacio los indices
	  K.clear();
	}

	// sino se esta en la condicion, simplemente copio la chica en la grande.
	else superImage( i + r, j + r, k + r) = imageSource(i,j,k);
      }

      J.clear();
    }

    I.clear();
  }
}


#endif /* FILE_BordStrategyMirror */

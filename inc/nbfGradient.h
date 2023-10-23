#ifndef FILE_Gradient
#define FILE_Gradient

#include <bs/nbfBordStrategyConst.h>
#include <bs/nbfBordStrategyMirror.h>
#include <blitz/array/stencil-et.h>

// Se usan MACROS para pasar las funciones a la nomeclatura Blitz++.

// Funcion Switch

template< class Pixel >
Pixel switchFunction( Pixel x, Pixel y )
{
  if ( x * y < 0 ) return (0);
  if ( abs(x) > abs(y) ) return (y);
  return (x);
}

BZ_DECLARE_FUNCTION2(switchFunction);


/************
 * Gradient *
 ************/

  /** Cálculo de gradientes.
      En esta clase se concentra todo el cálculo de diferencias sobre {\tt Array}s de cualquier dimensión. La idea es asociar un gradiente a una determinada imagen, de forma que todas las diferencias que queramos calcular las hagamos sobre una misma imagen. \\
      El cálculo de las diferencias se hace asumiendo alguna estrategia de borde. La misma puede especificarse en el constructor, o utilizar una por defecto. Para el manejo de la estrategia de borde lo que hacemos es generar una imagen mas grande, en cuyo centro copiamos la imagen original y generamos los bordes de acuerdo a la estrategia. El orden de derivadas que queramos calcular va a determinar el tamaño del stencil (vecindario), y esto determinara el tamaño de la imagen grande que vamos a generar. De acuerdo a esto, las diferencias se calcularán siempre sobre la imagen grande. \\
      En caso de que cambiemos los datos de la imagen original, debemos notificar al {\tt RefreshManager}, de modo que se reflejen los cambios en los bordes de la imagen grande, sino, estaremos calculando las diferencias sin haber actualizado los valores del borde de acuerdo a la estrategia adecuada. \\

      {\bf Métodos disponibles de {\tt Gradient}}. Todos los métodos de esta clase usan los esquemas de diferencias utilizados en los stencils de {\bf Blitz++}. Por esta razón utilizamos exactamente la misma nomeclatura que ellos, en favor de la extensibilidad y la claridad de la programación. \\
      Los nombres de los métodos de esta clase se forman de un sufijo, más el nombre del esquema en diferencias de {\bf Blitz++}. Todos los posibles esquemas pueden consultarse directamente en la documentación de {\bf Blitz++}, o en el código, más precisamente en el archivo {\tt blitz/array/stencilops.h}. Para cada uno de los esquemas, hay uno equivalente en esta clase. \\

      \begin{itemize}

      \item {\bf Diferencias Finitas} - cálculo de diferencias finitas de todos los órdenes y modalidades. Los métodos son de la forma: \\
      \begin{center} {\tt void forward11( Array< Pixel, Dim > \& ); } \end{center}
      reemplazando por el esquema en diferencias que se desee utilizar, en este caso se utilizó un esquema hacia adelante : {\tt forward}. También existen esquemas hacia atrás : {\tt backward} y centrados: {\tt central}. El primer número indica el orden de derivación (derivada primera (1), derivada segunda (2), etc.). El segundo número indica el orden de aproximación del esquema: primer orden (1), segundo orden (2), etcétera.

      \item {\bf Gradientes } - estos métodos calculan el gradiente de la imagen, escribiendo en un Array (unidimensional) de Arrays (multidimensionales), cada una de las componentes del gradiente. Los nombres de los métodos comienzan por {\tt grad}, y se les agrega el sufijo del tipo de diferencia (la primer letra con mayúscula): \\
      \begin{center} {\tt void gradForward11( Array< Array< Pixel, Dim >, 1 > \& ); } \end{center}

      \item {\bf Gradientes Normalizados} - cálculo de gradientes normalizados, es lo mismo que arriba pero dividiendo por el correspondiente módulo, es decir $\frac{\nabla I}{|\nabla I|}$. La nomeclatura es {\tt gradN} seguida del esquema de diferencias: \\
      \begin{center} {\tt void gradNForward11( Array< Array< Pixel, Dim >, 1 > \& ); } \end{center}

      \item {\bf Módulos L1 de Gradiente} - se calculan las {\tt Dim} componentes del gradiente y se halla su módulo L1, es decir: $\sum |\frac{\partial I}{\partial x_i}|$. La nomeclatura es {\tt modGradL1} seguido del esquema de diferencias:
      \begin{center} {\tt void modGradL1Forward11( Array< Array< Pixel, Dim >, 1 > \& ); } \end{center}

      \item {\bf Módulos L2 de Gradiente} - se calculan las {\tt Dim} componentes del gradiente y se halla su módulo L2, es decir: $\sqrt{\sum (\frac{\partial I}{\partial x_i})^2 }$. La nomeclatura es {\tt modGradL2} seguido del esquema de diferencias:
      \begin{center} {\tt void modGradL2Forward11( Array< Array< Pixel, Dim >, 1 > \& ); } \end{center}

      \end{itemize}
  */
template< class Pixel, int const Dim >
class Gradient : public GradientImp< Pixel >
{
  friend class GradientManager< Pixel, Dim >;

private:


  // Estrategia de borde.
  BordStrategy< Pixel, Dim > * ptrBordStrat;

  // Si es true quiere decir que gradiente es el encargado de crear y destruir la estrategia de borde. Si es false simplemente me pasan el puntero y no me preocupo de crearlo y destruirlo.
  bool createBordStrat;

  // Imagen agrandada incluyendo la estrategia de borde.
  Array< Pixel, Dim > superImage;

  // Vista central de la imagen agrandada.
  Array< Pixel, Dim > * superImageView;

  // imagen auxiliar.
  Array< Pixel, Dim > imAux;

  // Tamanio del vecindario sobre el cual se calcula la estrategia de borde.
  // Depende del orden de derivadas que estemos utilizando (del tamanio del stencil).
  unsigned vecindad;

  // Creadores y Destructores manejados exclusivamente por el GradientManager.

  Gradient();

  Gradient( Array< Pixel, Dim > & Image,
	    BordStrategy< Pixel, Dim > * BS,
	    unsigned vecindad );
  ~Gradient();

public:

  /** Construcción de Gradiente.  Crea un gradiente, asociado a una
      determinada imagen. A cada objeto gradiente asociamos una {\em
      Estrategia de Borde} y una {\em Imagen}. En este sentido, al
      crear un gradiente para una imagen, se crea además una imagen
      más grande que incluye la imagen original (en su centro), más
      una banda entorno al borde (para la estrategia de borde). Es por
      esto, que no tendrá sentido tener más de un objeto gradiente por
      cada imagen, pues, en caso de tenerlos, habría una super imagen
      por cada uno de ellos, con los perjuicios que esto implica en
      cuanto a requerimientos de memoria.  En este sentido, este
      método se encarga de que exista un único gradiente para una
      misma imagen. Éste utilizara una determinada estrategia de
      borde. Cuando se pretenda crear otro gradiente sobre la misma
      imagen, utilizando una estrategia de borde distinta, ésta última
      es la que tendrá sentido, y la anterior dejará de valer.  @param
      1 Especifica la imagen sobre la cual se calcularán las
      diferencias.  @param 2 Especifica la estrategia de borde a
      utilizar. Si no le pasamos nada, se toma una estrategia {\tt
      BordStrategyMirrorSimple}.  @param 3 Especifica el tamaño de la
      vecindad que queremos utilizar, estará determinado por el orden
      de derivadas que necesitemos. Evidentemente, si no escojemos
      ninguna estrategia, no podemos especificar este parámetro. Por
      defecto toma el valor 2.  */
  static Gradient< Pixel, Dim > * create( Array< Pixel, Dim > &,
					  BordStrategy< Pixel, Dim > * = 0,
					  unsigned = 2 );

  /** Destrucción de objetos gradiente.  Este proceso es controlado
      internamente por la clase {\tt GradientManager}. Como existe una
      sola instancia de la clase {\tt Gradiente}, ésta se destruirá,
      únicamente, cuando no queden referencias a la misma.  Si
      destruyeramos la instancia (usando {\tt delete}) cuando aún
      quedan referencias al objeto, evidentemente habría problemas.  */
  void destroy();

 /****************************************************************************
  * Finite differences.
  ****************************************************************************/

  Array<Pixel,Dim> & averaver( Array<Pixel,Dim> & A, const int dim ){
	  imAux = sqrt( pow2( blitz::forward11( *superImageView, dim ) ) + 
		        pow2( blitz::forward11( *superImageView, dim ) ) );
	  // A = sqrt( A );
	  return imAux;
  }

 /****************************************************************************
  * Gradient operators
  ****************************************************************************/

  // calculo de curvaturas

  /** Curvatura genérica. */
  void genericK( Array< Pixel, Dim > & ); // no pintó

  /** Curvatura Media */
  void meanK( Array< Pixel, Dim > & );

  /** Curvatura Gaussiana. */
  void gaussK( Array< Pixel, Dim > & );

  // OPERADORES NABLA

  /** Operador Nabla Mas.
      Este método implementa el operador $\nabla^+$ definido en el capítulo de implementación numérica de la documentación del proyecto {\bf Flujos}. Este operador da un esquema para el cálculo de $|\nabla I|$, para más detalles, ver la referencia anterior.

      @param 1 Imagen donde se guarda el resultado de aplicar el operador.
      @param 2 Orden con que se hace el cálculo del operador. Es posible hacer el cálculo para primer y para segundo orden, tomando 1 y 2, respectivamente.
  */
  void nablaForward( Array< Pixel, Dim > &, int );

  /** Operador Nabla Menos.
      Este método implementa el operador $\nabla^-$ definido en el capítulo de implementación numérica de la documentación del proyecto {\bf Flujos}. Este operador da un esquema para el cálculo de $|\nabla I|$, para más detalles, ver la referencia anterior.

      @param 1 Imagen donde se guarda el resultado de aplicar el operador.
      @param 2 Orden con que se hace el cálculo del operador. Es posible hacer el cálculo para primer y para segundo orden, tomando 1 y 2, respectivamente.
  */
  void nablaBackward( Array< Pixel, Dim > &, int );

  void nablaGodunov( Array< Pixel, Dim > & );

  /** Calculo de Curvaturas Principales.
      Este método implementa el operador $\nabla^+$ definido en el capítulo de implementación numérica de la documentación del proyecto {\bf Flujos}. Este operador da un esquema para el cálculo de $|\nabla I|$, para más detalles, ver la referencia anterior.

      @param 1 Imagen donde se guarda el resultado de aplicar el operador.
      @param 2 Orden con que se hace el cálculo del operador. Es posible hacer el cálculo para primer y para segundo orden, tomando 1 y 2, respectivamente.
  */
  void principalDirections( Array< Array< Pixel, 3 >, 1 > &,
			    Array< Array< Pixel, 3 >, 2 > & );

  void frobeniusDot( Array< Array< Pixel, Dim >, 2 > & H,
		     Array< Array< Pixel, Dim >, 1 > & V,
		     Array< Array< Pixel, Dim >, 1 > & W,
		     Array< Pixel, Dim > & F );

};


/******************************************************
 * IMPLEMENTACION
 ******************************************************/

/************
 * Gradient *
 ************/


template< class Pixel, int const Dim >
Gradient< Pixel, Dim > :: Gradient()
  : GradientImp< Pixel >()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("Gradient");
#endif

  FLUJOS_EXCEPTION("Datos incompletos : Este constructor no deberia invocarse.",
		   "Gradient<Pixel,Dim>::Gradient()");
}


template< class Pixel, int const Dim >
Gradient< Pixel, Dim > :: Gradient( Array< Pixel, Dim > & imageSource,
				    BordStrategy< Pixel, Dim > * ptrStrat,
				    unsigned vec )
  : GradientImp< Pixel >(),
    vecindad( vec ),
    imAux( imageSource.shape() )
{

#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("Gradient");
#endif

  // Estrategia de bordes
  if ( ptrStrat == 0 ){
    ptrBordStrat = new BordStrategyMirrorSimple< Pixel, Dim >();
    createBordStrat = true;
  }
  else{
    ptrBordStrat = ptrStrat;
    createBordStrat = false;
  }

  // aca se CREA la SuperImage y se inicializa correctamente.
  ptrBordStrat->generateSuper( imageSource, superImage, vecindad );
  superImageView = &imageSource;
}


template< class Pixel, int const Dim >
Gradient< Pixel, Dim > * Gradient< Pixel, Dim > :: create( Array< Pixel, Dim > & Image,
							   BordStrategy< Pixel, Dim > * BS = 0,
							   unsigned vecindad = 2 )
{
  // busco si ya hay un gradiente para esta imagen
  Gradient< Pixel, Dim > * result = GradientManager< Pixel, Dim > :: askForGradient( &Image );

  if ( result == 0 ){
    // si no hay, creo uno nuevo.
    result = new Gradient< Pixel, Dim >( Image, BS, vecindad );
    GradientManager< Pixel, Dim > :: subscribe( result );
    return result;
  }
  else{ // si ya hay un gradiente creado para la imagen:

    // chequeo la estrategia de borde. Si BS = 0, dejo la estrategia existente.
    // sino, cambio la vieja por la nueva
    if ( BS != 0 ){
      if ( BS->getId() != result->ptrBordStrat->getId() ){

	// borro la vieja estrategia (si es que la cree yo)
	if ( result->createBordStrat )
	  delete result->ptrBordStrat; // se llama solo el unsubscribe
	else{
	  // sino, debo llamar a prepo el unsubscribe
	  RefreshManager< Pixel, Dim > :: unsubscribe( result->ptrBordStrat );
	}

	// asigno los nuevos atributos
	result->ptrBordStrat = BS;
	result->createBordStrat = false;
	result->ptrBordStrat->generateSuper( Image, result->superImage, result->vecindad );
	result->superImageView = &Image; // no deberia cambiar, pero por las dudas ...
      }
    }

    // chequeo el tamanio de la vecindad
    if ( vecindad > result->vecindad ) { // si es mas grande, tengo que recrear la superImage
      // asigno la nueva vecindad
      result->vecindad = vecindad;

      // genero la nueva superImage y actualizo la vista.
      result->ptrBordStrat->generateSuper( Image, result->superImage, result->vecindad );
      result->superImageView = &Image;
    }
  }
  return result;
}


template< class Pixel, int const Dim >
void Gradient< Pixel, Dim > :: destroy()
{
  GradientManager< Pixel, Dim > :: unsubscribe( this );
}


template< class Pixel, int const Dim >
Gradient< Pixel, Dim > :: ~Gradient()
{
  if ( superImage.size() != 0 ) superImage.free();
  if ( imAux.size() != 0 ) imAux.free();

  // Si cree la estrategia la destruyo.
  if ( createBordStrat && ( ptrBordStrat != 0 ) ) delete ptrBordStrat;

#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  DELETING("Gradient");
#endif
}


// Curvaturas


template<class Pixel, int const Dim>
void Gradient< Pixel, Dim > :: genericK( Array< Pixel, Dim > & K )
{
/*
  Array< Array< Pixel, Dim >, 1 > grad( Dim );
  for ( int i = firstDim; i < Dim; i++ )
    grad( i ).resize( K.shape() );
  gradNCentral12n( grad );
  K = 0;
  for ( int i = firstDim; i < Dim; i++ )
    K += blitz :: central12n( grad( i ), i );
*/
}


template<class Pixel, int const Dim>
void Gradient< Pixel, Dim > :: meanK( Array< Pixel, Dim > & K )
{
  GradientImp< Pixel > :: meanKImp( *superImageView, K );
}


template<class Pixel, int const Dim>
void Gradient< Pixel, Dim > :: gaussK( Array< Pixel, Dim > & K )
{
  GradientImp< Pixel > :: gaussKImp( *superImageView, K );
}


BZ_DECLARE_STENCIL_OPERATOR1(nablaForward2D,A)
	return sqrt( pow2( max( backward11n(A,firstDim), 0 ) ) + 
	             pow2( min( forward11n(A,firstDim), 0 ) ) +
	             pow2( max( backward11n(A,secondDim), 0 ) ) + 
				 pow2( min( forward11n(A,secondDim), 0 ) ) );
BZ_END_STENCIL_OPERATOR

BZ_ET_STENCIL(nablaForward2D,P_numtype)

// Nablas

template< class Pixel, int const Dim >
void Gradient< Pixel, Dim > :: nablaForward( Array<Pixel, Dim> & B, int orden = 1 )
{
  B = 0;

  if ( orden == 1 ){

    for ( int i = firstDim; i < Dim; i++ ){

      imAux = blitz :: backward11n( *superImageView, i );
      B += where( imAux > 0, pow2(imAux), 0 );

      imAux = blitz :: forward11n( *superImageView, i );
      B += where( imAux < 0, pow2(imAux), 0 );
    }
  }

  else{

    for ( int i = firstDim; i < Dim; i++ ){

      // A, C, E, ... ( del sethian )
      imAux = blitz :: backward11n( *superImageView, i )
	+ recip_2 * switchFunction( blitz :: backward22n( *superImageView, i ),
				    blitz :: central22n( *superImageView, i ) );

      B += where( imAux > 0, pow2(imAux), 0 ); // max(x2,0)

      // B, D, F, ... ( del sethian )
      imAux = blitz :: forward11n( *superImageView, i )
	- recip_2 * switchFunction( blitz :: forward22n( *superImageView, i ),
				    blitz :: central22n( *superImageView, i ) );

      B += where( imAux < 0, pow2(imAux), 0 ); // min(x2,0)
    }
  }
  B = sqrt( B );
}


template< class Pixel, int const Dim >
void Gradient< Pixel, Dim > :: nablaBackward( Array< Pixel, Dim > & B, int orden = 1 )
{
  B = 0;

  if ( orden == 1 ){

    for ( int i = firstDim; i < Dim; i++ ){

      imAux = blitz :: forward11n( *superImageView, i );
      B += where( imAux > 0, pow2(imAux), 0 ); // max(x2,0)

      imAux = blitz :: backward11n( *superImageView, i );
      B += where( imAux < 0, pow2(imAux), 0 ); // min(x2,0)
    }
  }

  else{				// orden 2

    for ( int i = firstDim; i < Dim; i++ ){

      // A, C, E, ... ( del sethian )
      imAux = blitz :: backward11n( *superImageView, i );
      + recip_2 * switchFunction( blitz :: backward22n( *superImageView, i ),
				  blitz :: central22n( *superImageView, i ) );

      B += where( imAux < 0, pow2(imAux), 0 ); // min(x2,0)

      // B, D, F, ... ( del sethian )
      imAux = blitz :: forward11n( *superImageView, i );
      - recip_2 * switchFunction( blitz :: forward22n( *superImageView, i ),
				  blitz :: central22n( *superImageView, i ) );

      B += where( imAux > 0, pow2(imAux), 0 ); // max(x2,0)
    }
  }
  B = sqrt( B );
}


template< class Pixel1, class Pixel2 >
Pixel1 maxSqr( Pixel1 x, Pixel2 y )
{
  // z = max^2( x, -y, 0 )
  if ( ( x > 0 ) || ( y < 0 ) )
    if ( x > -y ) return ( pow2(x) );
    else return ( pow2(y) );
  else return (0);
}

BZ_DECLARE_FUNCTION2(maxSqr)


template< class Pixel, int const Dim >
void Gradient< Pixel, Dim > :: nablaGodunov( Array< Pixel, Dim > & B )
{
  // Este operador implementa:
  //
  // sqrt( max( max( DerMenosX, 0)^2, min( DerMasX, 0)^2 )
  //     + max( max( DerMenosY, 0)^2, min( DerMasY, 0)^2 ) )
  //
  // Que es lo mismo que:
  //
  // sqrt( max( DerMenosX, -DerMasX, 0 )^2 + max( DerMenosY, -DerMasY, 0 )^2 )
  //

  B = 0;

  for ( int i = firstDim; i < Dim; i++ ){
    // imAux = max^2( backward11n, -forward11n, 0 )
    imAux = maxSqr( blitz :: backward11n( *superImageView, i ),
		    blitz :: forward11n( *superImageView, i ) );
    B += imAux;
  }
  B = sqrt( B );
}


template< class Pixel, int const Dim >
void Gradient< Pixel, Dim > :: principalDirections( Array< Array< Pixel, 3 >, 1 > & curvatures,
						    Array< Array< Pixel, 3 >, 2 > & directions )
{
  Array< Array< Pixel, 3 >, 1 > h(3);
  Array< Array< Pixel, 3 >, 1 > f(3);
  Array< Array< Pixel, 3 >, 1 > gradient(3);
  
  for ( int i = 0; i < 3; i ++ ){
    gradient(i).resize( superImageView->shape() );
    h(i).resize( superImageView->shape() );
    f(i).resize( superImageView->shape() );
  }
  
  Array< Pixel, 3 > gamma( superImageView->shape() );
  Array< Pixel, 3 > delta( superImageView->shape() );
  
  this->gradCentral14n( gradient );
  gamma = sqrt( pow2( gradient(0) ) + pow2( gradient(1) ) );
  delta = sqrt( pow2( gamma ) + pow2( gradient(2) ) );
  
  //  cout << "sqrt(I_x^2+I_y^2+I_z^2) = 0 : " << count( where( delta == 0, 1, 0 ) ) << endl;

  h(0) = gradient(1) / gamma;
  h(1) = (-1) * gradient(0) / gamma;
  h(2) = 0;
  
  f(0) = gradient(2) * gradient(0) / gamma / delta;
  f(1) = gradient(2) * gradient(1) / gamma / delta;
  f(2) = (-1) * gamma / delta;
  
  Array< Array< Pixel, 3 >, 2 > H(3,3);
  for ( int i = 0; i < 3; i ++ ){
    for ( int j = 0; j < 3; j ++ ){
      H(i,j).resize( superImageView->shape() );
    }
  }
  
  OperadoresDiferenciales< Pixel, 3 > op( gradient );
  op.jacobianCentral14n( H );
  
  this->frobeniusDot( H, h, h, gamma );	// gamma = h' H h
  this->frobeniusDot( H, f, f, delta );	// delta = f' F f
  
  Array< Pixel, 3 > cruz( superImageView->shape() );
  this->frobeniusDot( H, h, f, cruz ); // cruz = h' H f
  
  //  cout << "f'Hh = 0 : " << count( where( cruz == 0, 1, 0 ) ) << endl;

  // modGrad = |\/g|
  Array< Pixel, 3 > modGrad( superImageView->shape() );
  //  this->modGradL2Central14n( modGrad );
  modGrad = 1;

  curvatures(0) = ( gamma + delta + sqrt( pow2( gamma - delta ) + 4 * pow2( cruz ) ) ) / 2;
  // / ( modGrad + EPSILON );
  curvatures(1) = ( gamma + delta - sqrt( pow2( gamma - delta ) + 4 * pow2( cruz ) ) ) / 2;
  // / ( modGrad + EPSILON );
  
  for ( int i = 0; i < 2; i++ ){
    for ( int j = 0; j < 3; j++ ){
      directions(i,j) = h(j) + f(j) * ( modGrad * curvatures(i) - gamma ) / cruz;
    }
  }
  
  // normalizacion de las direcciones
  Array< Pixel, 3 > norm( superImageView->shape() );
  for ( int i = 0; i < 2; i++ ){
    norm = 0;
    for ( int j = 0; j < 3; j++ ){
      norm += pow2( directions(i,j) );
    }
    for ( int j = 0; j < 3; j++ ){
      directions(i,j) /= sqrt( norm );
    }
  }

  // elimino los puntos mal condicionados (los pongo en 0)
  for ( int i = 0; i < 2; i++ ){
    for ( int j = 0; j < 3; j++ ){
      directions(i,j) = where( ( abs( gamma ) > EPSILON ) &&
			       ( abs( delta ) > EPSILON ) &&
			       ( abs( cruz ) > EPSILON ),
			       directions(i,j),
			       0 );
    }
  }
  int counter = count( where( directions(0,0) + directions(0,1) + directions(0,2) == 0 , 1, 0 ) );
  cout << "Ill posed points (1st field): " << counter << endl;
  counter = count( where( directions(1,0) + directions(1,1) + directions(1,2) == 0 , 1, 0 ) );  cout << "Ill posed points (2nd field): " << counter << endl;

}

template< class Pixel, int const Dim >
void Gradient< Pixel, Dim > :: frobeniusDot( Array< Array< Pixel, Dim >, 2 > & H,
					     Array< Array< Pixel, Dim >, 1 > & V,
					     Array< Array< Pixel, Dim >, 1 > & W,
					     Array< Pixel, Dim > & F )
{
  int lado = H.shape()(0);
  F = 0;
  for ( int i = 0; i < lado; i++ ){
    for ( int j = 0; j < lado; j++ ){
      F += H(i,j) * V(i) * W(j);
    }
  }
}


#endif /* FILE Gradient */

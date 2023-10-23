#ifndef FILE_BordStartegy //
#define FILE_BordStartegy

#include <MetaBounds.hh>

// BordStrategy Id
#define BS_ID_CONST 0
#define BS_ID_MIRROR_SIMPLE 1
#define BS_ID_MIRROR_DOUBLE 2
	    
/****************
 * BordStrategy *
 ****************/

/** Estrategia de Borde.  Con esta clase nos referimos a las distintas
    estrategias de borde que podemos utilizar. Cada uno de los casos
    particulares, serán clases que heredan de ésta, redefiniendo el
    método {\tt refresh()}.  En general podemos utilizar una misma
    estartegia de borde, es decir una misma instancia de un objeto de
    esta clase, para varias imágenes. La idea es simplemente utilizar
    el método que se encarga de actualizar los datos sobre el borde
    cuando se ha producido algún cambio en la imagen original.  La
    clase tiene métodos que permiten manipular la imagen agrandada
    (datos + borde), de acuerdo a las necesidades de cada caso, pero
    principalmente una función de actualización, que permite
    recalcular los datos del borde, una vez que han cambiado los datos
    sobre la imagen original.  A estos efectos hay un objeto de
    control, de la clase {\tt RefreshManager}, que lleva un registro
    de que imágenes están asociadas a que estartegias, de forma de
    poder realizar la actualización correctamente. Por más detalles
    sobre este procedimiento, consultar la documentación de dicha
    clase.  */
template< class Pixel, int const Dim >
class BordStrategy
{
protected:

  Array< Pixel, Dim > * imageSource; // pointer to source array
  Array< Pixel, Dim > superImage;	 // internal enlarged array
  
  unsigned padSize;                  // size of side padding

  TinyVector< int, Dim > lboundDim;
  TinyVector< int, Dim > uboundDim;

public:

  /** Genera una imagen extendida.
      Esta función crea una imagen que incluye la imagen original en su centro, y los bordes se rellenan con la estrategia de borde deseada. En general, esta función se llama una única vez, y luego simplemente se utiliza la función de actualización. Sin embargo, el hecho de no conocer sobre que imagen se va a trabajar, impide llamar este método directamente desde del constructor, pero como contrapartida podemos utilizar una misma instancia para varias imágenes.
  */
  void generateSuper();
  
  /** Actualiza la estrategia de borde.  Esta función se encarga de
      actualizar la imagen agrandada, a partir de la imagen original y
      la estrategia de borde correspondiente. Cada una de las
      estartegias de borde, clases que heredan de {\tt BordStrategy}
      implementarán adecuadamente este método, que aquí se define como
      virtual puro.  */
  virtual void refresh() = 0;

  /** Identificación.  Identifica cada una de las estrategias de borde
      que heredan de esta clase con una constante de la forma
      BS\_ID\_xxx. Es utilizada por el {\tt GradientManager} para
      comparar las distintas estrategias y saber si son del mismo
      tipo. */
  virtual unsigned getId() = 0;

  /** Nueva Estrategia de Borde.  El constructor no tiene parámetros,
      pues es independiente de la imagen que esté considerando. Puedo
      utilizar la misma estrategia para varias imágenes.  */
  BordStrategy( Array< Pixel, Dim > &, unsigned );
  
  virtual ~BordStrategy();
};


/******************************************************
 * IMPLEMENTACION
 ******************************************************/

/****************
 * BordStrategy *
 ****************/

template< class Pixel, int const Dim >
void BordStrategy< Pixel, Dim > :: generateSuper()
{
  // creo la imagen grande con el valor dado de padSize
  TinyVector< int, Dim > rangeDim;
  rangeDim = this->imageSource->shape();
  rangeDim += ( 2 * this->padSize );
  
  // creo la imagen y asigno el artributo.
  this->superImage.resize( rangeDim );
  
  // lleno el vector de bounds
  TinyVector< int, Dim > pad;
  pad = this->padSize;
  this->lboundDim = this->superImage.lbound() + pad;
  this->uboundDim = this->superImage.ubound() - pad;

  // Inicializo el objeto de Indexacion
  RectDomain< Dim > rDomain( this->lboundDim, this->uboundDim );

  // pongo en el centro de la super los datos de la source
  superImage( rDomain ) = ( *imageSource );

  // libero la imagen original y le re-asigno el centro de la superImage
  imageSource->free();
  
  // re-assign the original image so it refers to the center of the enlarged one
  Array< Pixel, Dim > Ai = superImage( rDomain );
  imageSource->reference( Ai );
}


template< class Pixel, int const Dim>
BordStrategy< Pixel, Dim > :: BordStrategy( Array< Pixel, Dim > & image_o,
										    unsigned padSize = 2 )
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("BordStrategy");
#endif

  // assign attributes
  this->imageSource = &image_o;
  this->padSize = padSize;

  this->generateSuper();
}


template< class Pixel, int const Dim>
BordStrategy< Pixel, Dim > :: ~BordStrategy()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  DELETING("BordStrategy");
#endif  
}


#endif /* FILE_BordStrategy */

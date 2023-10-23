#ifndef FILE_BordStartegyConst
#define FILE_BordStartegyConst

#include <bs/nbfBordStrategy.h>

/*********************
 * BordStrategyConst *
 *********************/

/** Estrategia de Bordes Constante.  Esta estrategia consiste en
    extender la imagen, poniendo un valor constante en todos los
    pixeles del borde. El valor de la constante puede especificarse en
    el constructor de la clase.  */
template< class Pixel, int const Dim >
class BordStrategyConst :  public BordStrategy< Pixel, Dim >
{
public:
  
/** Implementación.  Este método implementa la actualización de la
    estrategia de bordes. Simplemente copia el valor de la constante
    en todos los pixeles del borde.  */
  void refresh();

  unsigned getId(){ return BS_ID_CONST;}

  /** Constructor.  Constructor que toma como parametro el valor
      constante de la estrategia de borde. No existe otra forma de
      especificar la constante que no sea a través de este método.  */
  BordStrategyConst( Array< Pixel, Dim > &, unsigned, Pixel );
  ~BordStrategyConst();
  
private:  
  Pixel valor;			// atributo para el valor constante de los bordes.

};


/******************************************************
 * IMPLEMENTACION
 ******************************************************/


/*********************
 * BordStrategyConst *
 *********************/

template< class Pixel, int const Dim >
BordStrategyConst< Pixel, Dim > :: BordStrategyConst( Array< Pixel, Dim > & A, unsigned offset = 2, Pixel val = 0 )
: BordStrategy< Pixel, Dim >( A, offset ), valor(val)
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  CREATING("BordStrategyConst");
#endif
}


template< class Pixel, int const Dim >
BordStrategyConst< Pixel, Dim > :: ~BordStrategyConst()
{
#ifdef FLUJOS_WITH_REFERENCE_COUNTING
  DELETING("BordStrategyConst");
#endif
}


template< class Pixel, int const Dim >
void BordStrategyConst< Pixel, Dim > :: refresh()
{
  typename Array< Pixel, Dim > :: iterator begin = this->superImage.begin();
  typename Array< Pixel, Dim > :: iterator end = this->superImage.end();
  TinyVector< int, Dim > pos;
  
//  TinyVector< int, Dim > upperBound( this->superImage.shape() - this->padSize - 1 );
//  TinyVector< int, Dim > lowerBound( this->padSize );
  
  while ( begin != end ){
    pos = begin.position();
    if ( !checkBounds( pos, this->lboundDim, this->uboundDim ) )
      // si esta en la banda exterior, pongo el valor de la constante
      (*begin) = valor;
    ++begin;
  }
}


#endif /* FILE_BordStrategyConst */

#ifndef FILE_MetaBounds
#define FILE_MetaBounds

/* Meta Funcion que permite saber si un determinado TV esta dentro de
   una Bounding Region dada por una par de TVs. Retorna TRUE si esta
   dentro de la region o en el borde de la misma. */

template< int J >
struct meta_bounds
{
  template< class Position,  int   Dim >
  static  bool f(   const TinyVector< Position, Dim >& x, 
		    const TinyVector< Position, Dim >& lowerBound,
		    const TinyVector< Position, Dim >& upperBound )
  {
    return( ( x(J) >= lowerBound(J) ) && ( x(J) <= upperBound(J) ) 
	    && meta_bounds< J-1 > :: f( x, lowerBound, upperBound ) );    
  }
};


template< >
struct meta_bounds< 0 >
{
  template< class Position, int   Dim >
  static  bool f(   const TinyVector< Position, Dim >& x, 
		    const TinyVector< Position, Dim >& lowerBound, 
		    const TinyVector< Position, Dim >& upperBound )
  {
    return( ( x(0) >= lowerBound(0) ) && ( x(0) <= upperBound(0) ) );
  }
};

//---------------------------------------------------------------------------------

template< class Position, int   Dim >
inline bool checkBounds(   const TinyVector< Position, Dim >& x, 
			   const TinyVector< Position, Dim >& lowerBound,
			   const TinyVector< Position, Dim >& upperBound )
  
{
 return  meta_bounds< Dim - 1 > :: f( x, lowerBound, upperBound );
}

#endif /* FILE_MetaBounds */

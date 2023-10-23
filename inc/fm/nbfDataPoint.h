#ifndef FILE_nbfDataPoint
#define FILE_nbfDataPoint

template< class Pixel, int const Dim >
class nbfDataPoint
{
public:

	Pixel * distance;
	TinyVector< int, Dim > position;

	nbfDataPoint(){}

	nbfDataPoint( const TinyVector< int, Dim > & tv )
	{
        this->position = tv;  
		this->distance == 0;
	}

	nbfDataPoint( const nbfDataPoint< Pixel, Dim > & PuntoDatoRef )
	{
        this->position = PuntoDatoRef.position;
		this->distance = PuntoDatoRef.distance;
	}
  
	nbfDataPoint< Pixel, Dim > operator = ( const nbfDataPoint< Pixel, Dim > & input )
    {
	  this->position = input.position;
	  this->distance = input.distance;
      return (*this);
    }

~nbfDataPoint(){};
};

#endif // FILE_nbfDataPoint

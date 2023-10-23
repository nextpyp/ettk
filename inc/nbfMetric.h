#ifndef FILE_nbfMetric
#define FILE_nbfMetric

class nbfMetric
{
public:

	// constructor takes weight array as input
	nbfMetric(){};

	virtual ~nbfMetric(){};

	// scalar types
	template< class Scalar >
	static inline Scalar mod( Scalar a ){ return sqrt( nbfMetric::modSqr(a) ); }

	template< class Scalar >
	static inline Scalar mod( Scalar a, Scalar b ){ return sqrt( nbfMetric::modSqr(a,b) ); }

	template< class Scalar >
	static inline Scalar modSqr( Scalar a ){ return a*a; }

	template< class Scalar >
	static inline Scalar modSqr( Scalar a, Scalar b ){ return pow2(a-b); }

	template< class Scalar >
	static inline Scalar sum( Scalar a, Scalar b ){ return (a+b); }

	template< class Scalar >
	static inline Scalar sub( Scalar a, Scalar b ){ return (a-b); }

	// vectorial types
	template< class Scalar, int const Dim >
		static inline Scalar mod( TinyVector< Scalar, Dim > & a ){ return sqrt( nbfMetric::modSqr(a) ); }

	template< class Scalar, int const Dim >
	static inline Scalar mod( TinyVector< Scalar, Dim > & a, TinyVector< Scalar, Dim > & b ){ 
		return sqrt( nbfMetric::modSqr(a,b) );
	}

	template< class Scalar, int const Dim >
	static inline Scalar modSqr( TinyVector< Scalar, Dim > & a ){ return dot(a,a); }

	template< class Scalar, int const Dim >
	static inline Scalar modSqr( TinyVector< Scalar, Dim > & a, TinyVector< Scalar, Dim > & b ){

		TinyVector< Scalar, Dim > bn;
		bn = b / sqrt(dot(b,b));
		return pow2( 1.0 - fabs(dot(a,bn)) );

		TinyVector< Scalar, Dim > dir;
		dir = a - b;
		dir = dir / sqrt(dot(dir,dir));
		return ( 1.0 - sqrt(pow2(dot(a,bn))) ) + .1 * ( 1 - fabs( dot(dir,bn) ) );

		// define anisotropy weight in the FA term
		Scalar w = 2.0;
        Scalar fa, fb, theta, tmp, mod;
		fa = sqrt(dot(a,a));
		fb = sqrt(dot(b,b));
   		if ( dot(a,b) > 0 ){
			tmp = dot(a/fa,b/fb);
			if ( fa > fb ){
				mod = log(fa/fb);
			}
			else{
				mod = log(fb/fa);
			}
		}
		else{
			tmp = dot(a/fa,-b/fb);
			if ( fa > fb ){
				mod = log(fa/fb);
			}
			else{
				mod = log(fb/fa);
			}
		}
		if ( tmp < 1 ){
			theta = fabs( acos( tmp ) );
		}
		else{
			theta = 0;
		}
        return (w * mod + theta);

		// isotropic metric (NOT LONGER USED)
		if ( dot(a,b) > 0 ){
			return dot(a-b,a-b); 
		}
		else{
			return dot(a+b,a+b);
		}
	}

	template< class Scalar, int const Dim >
	static inline TinyVector< Scalar, Dim > sum( TinyVector< Scalar, Dim > & a, TinyVector< Scalar, Dim > & b ){ 
		if ( dot(a,b) > 0 ){
			return (a+b);
		}
		else{
			return (a-b);
		}
	}

	template< class Scalar, int const Dim >
	static inline TinyVector< Scalar, Dim > sub( TinyVector< Scalar, Dim > & a, TinyVector< Scalar, Dim > & b ){ 
		if ( dot(a,b) > 0 ){
			return (a-b);
		}
		else{
			return (a+b);
		}
	}

	// Tensor types
	template< class Scalar >
	static inline Scalar mod( Array< Scalar, 2 > & a ){ return sum(sqrt(a*a)); }

	template< class Scalar >
	static inline Scalar mod( Array< Scalar, 2 > & a, Array< Scalar, 2 > & b ){ 
		return sum(sqrt(a*b));
	}

	template< class Scalar >
	static inline Scalar modSqr( Array< Scalar, 2 > & a ){ return sum(a*a); }

	template< class Scalar >
	static inline Scalar modSqr( Array< Scalar, 2 > & a, Array< Scalar, 2 > & b ){ 
		return sum((a-b)*(a-b)); 
	}

	template< class Scalar >
	static inline Array< Scalar, 2 > sum( Array< Scalar, 2 > & a, Array< Scalar, 2 > & b ){ 
		return (a+b);
	}

	template< class Scalar >
	static inline Array< Scalar, 2 > sub( Array< Scalar, 2 > & a, Array< Scalar, 2 > & b ){ 
		return (a-b);
	}
};

#endif /* FILE_nbfMetric */

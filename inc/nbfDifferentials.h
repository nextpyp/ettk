#ifndef FILE_nbfDifferentials
#define FILE_nbfDifferentials

#include "nbf-stencil-et.h"
#include <vtkMath.h>

#define EPSILON 1.19209290E-07F

#define BZ_DECLARE_STENCIL_OPERATOR4(name,A,B,C,D) \
  template<class T>                              \
  inline _bz_typename T::T_numtype name(T& A, T& B, T& C, T& D)    \
  {

BZ_DECLARE_STENCIL_OPERATOR1(nablaPlus2D,A)
return sqrt( pow2( blitz::minmax::max( backward11n(A,firstDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,firstDim), 0 ) ) +
			 pow2( blitz::minmax::max( backward11n(A,secondDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,secondDim), 0 ) ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(nablaMinus2D,A)
return sqrt( pow2( blitz::minmax::max( forward11n(A,firstDim), 0 ) ) + 
	         pow2( blitz::minmax::min( backward11n(A,firstDim), 0 ) ) +
			 pow2( blitz::minmax::max( backward11n(A,secondDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,secondDim), 0 ) ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(nablaPlus3D,A)
return sqrt( pow2( blitz::minmax::max( backward11n(A,firstDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,firstDim), 0 ) ) +
			 pow2( blitz::minmax::max( backward11n(A,secondDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,secondDim), 0 ) ) +
			 pow2( blitz::minmax::max( backward11n(A,thirdDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,thirdDim), 0 ) ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(nablaMinus3D,A)
return sqrt( pow2( blitz::minmax::max( forward11n(A,firstDim), 0 ) ) + 
	         pow2( blitz::minmax::min( backward11n(A,firstDim), 0 ) ) +
			 pow2( blitz::minmax::max( backward11n(A,secondDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,secondDim), 0 ) ) +
			 pow2( blitz::minmax::max( backward11n(A,thirdDim), 0 ) ) + 
	         pow2( blitz::minmax::min( forward11n(A,thirdDim), 0 ) ) );
BZ_END_STENCIL_OPERATOR

// nablaGodunov
//
// sqrt( max( max( DerMenosX, 0)^2, min( DerMasX, 0)^2 )
//     + max( max( DerMenosY, 0)^2, min( DerMasY, 0)^2 ) )
//
// Que es lo mismo que:
//
// sqrt( max( DerMenosX, -DerMasX, 0 )^2 + max( DerMenosY, -DerMasY, 0 )^2 )
 
BZ_DECLARE_STENCIL_OPERATOR1(nablaGodunov2D,A)
return sqrt( pow2( blitz::minmax::max( 
			       blitz::minmax::max( backward11n(A,firstDim), -forward11n(A,firstDim) )
				   , 0 ) ) 
			+ 
	         pow2( blitz::minmax::max(
			       blitz::minmax::max( backward11n(A,secondDim), -forward11n(A,secondDim) )
				   , 0 ) ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(nablaGodunov3D,A)
return sqrt( pow2( blitz::minmax::max( 
			       blitz::minmax::max( backward11n(A,firstDim), -forward11n(A,firstDim) )
				   , 0 ) ) 
			+ 
	         pow2( blitz::minmax::max(
			       blitz::minmax::max( backward11n(A,secondDim), -forward11n(A,secondDim) )
				   , 0 ) )
			+ 
	         pow2( blitz::minmax::max(
			       blitz::minmax::max( backward11n(A,thirdDim), -forward11n(A,thirdDim) )
				   , 0 ) ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(meanK2D,A)
return ( ( central22n(A,firstDim) * pow2( central12n(A,secondDim) )
	     - 2 * central12n(A,firstDim ) * central12n(A,secondDim)
	         * mixed22n(A,firstDim,secondDim )
		 + central22n(A,secondDim) * pow2( central12n(A,firstDim) ) )
         / sqrt( pow3( ( pow2( central12n(A,firstDim) )
		               + pow2( central12n(A,secondDim) ) ) ) + EPSILON ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(meanK3D,A)
return ( (
		( central22n(A,secondDim) + central22n(A,thirdDim) ) * pow2( central12n(A,firstDim) )
      + ( central22n(A,firstDim) + central22n(A,thirdDim) ) * pow2( central12n(A,secondDim) )
      + ( central22n(A,firstDim) + central22n(A,secondDim) ) * pow2( central12n(A,thirdDim) )
       - 2 * central12n(A,firstDim) * central12n(A,secondDim) * mixed22n(A,firstDim,secondDim)
       - 2 * central12n(A,firstDim) * central12n(A,thirdDim) * mixed22n(A,firstDim,thirdDim)
       - 2 * central12n(A,secondDim) * central12n(A,thirdDim) * mixed22n(A,secondDim,thirdDim)
      ) / sqrt( pow3( ( pow2( central12n(A,firstDim) )
	                  + pow2( central12n(A,secondDim) )
	                  + pow2( central12n(A,thirdDim) ) ) ) + EPSILON ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR2(expansion2D,A,F)
return ( blitz::minmax::max(1*F,0) * nablaPlus2D(A)
	   + blitz::minmax::min(1*F,0) * nablaMinus2D(A) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR2(expansion3D,F,A)
return ( blitz::minmax::max(F,0) * nablaPlus3D(A)
	   + blitz::minmax::min(F,0) * nablaMinus3D(A) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(curvature2D,A)
return ( ( central22n(A,firstDim) * pow2( central12n(A,secondDim) )
	     - 2 * central12n(A,firstDim ) * central12n(A,secondDim)
	         * mixed22n(A,firstDim,secondDim )
		 + central22n(A,secondDim) * pow2( central12n(A,firstDim) ) )
         / ( pow2( central12n(A,firstDim) ) + pow2( central12n(A,secondDim) ) + EPSILON ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR1(curvature3D,A)
return ( (
		( central22n(A,secondDim) + central22n(A,thirdDim) ) * pow2( central12n(A,firstDim) )
      + ( central22n(A,firstDim) + central22n(A,thirdDim) ) * pow2( central12n(A,secondDim) )
      + ( central22n(A,firstDim) + central22n(A,secondDim) ) * pow2( central12n(A,thirdDim) )
       - 2 * central12n(A,firstDim) * central12n(A,secondDim) * mixed22n(A,firstDim,secondDim)
       - 2 * central12n(A,firstDim) * central12n(A,thirdDim) * mixed22n(A,firstDim,thirdDim)
       - 2 * central12n(A,secondDim) * central12n(A,thirdDim) * mixed22n(A,secondDim,thirdDim)
      ) / ( pow2( central12n(A,firstDim) )
	      + pow2( central12n(A,secondDim) )
	      + pow2( central12n(A,thirdDim) ) + EPSILON ) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR3(advection2D,A,Ux,Uy)
return ( blitz::minmax::max(1*Ux,0) * backward11n(A,firstDim)
	   + blitz::minmax::min(1*Ux,0) * forward11n(A,firstDim)
	   + blitz::minmax::max(1*Uy,0) * backward11n(A,secondDim)
	   + blitz::minmax::min(1*Uy,0) * forward11n(A,secondDim) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR4(advection3D,A,Ux,Uy,Uz)
return ( blitz::minmax::max(1*Ux,0) * backward11n(A,firstDim)
	   + blitz::minmax::min(1*Ux,0) * forward11n(A,firstDim)
	   + blitz::minmax::max(1*Uy,0) * backward11n(A,secondDim)
	   + blitz::minmax::min(1*Uy,0) * forward11n(A,secondDim)
	   + blitz::minmax::max(1*Uz,0) * backward11n(A,thirdDim)
	   + blitz::minmax::min(1*Uz,0) * forward11n(A,thirdDim) );
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR2(diffusion2D,G,A)
return central12n(G,firstDim) * central12n(A,firstDim) + G * central22n(A,firstDim) +
	   central12n(G,secondDim) * central12n(A,secondDim) + G * central22n(A,secondDim);
BZ_END_STENCIL_OPERATOR

BZ_DECLARE_STENCIL_OPERATOR2(diffusion3D,G,A)
return central12n(G,firstDim) * central12n(A,firstDim) + G * central22n(A,firstDim) +
	   central12n(G,secondDim) * central12n(A,secondDim) + G * central22n(A,secondDim) +
       central12n(G,thirdDim) * central12n(A,thirdDim) + G * central22n(A,thirdDim);
BZ_END_STENCIL_OPERATOR

// Redistance Function Helpers

// Sign function from Sussman & Fatemi
float sussmanSign( float x )
{
	float Dx = 1;
	if ( x > Dx ) return (1);
	if ( x < -Dx ) return (-1);
	return ( x / Dx + 1 / vtkMath::Pi() * sin( vtkMath::Pi() * x / Dx ) );
}

BZ_DECLARE_FUNCTION(sussmanSign)

// Gradient function from Sussman & Fatemi
template<class T1>
inline _bz_typename T1::T_numtype
sussmanGrad(T1& A, T1& S, int dim) {
  T1::T_numtype forward = forward11n(A,dim);
  T1::T_numtype backward = backward11n(A,dim);
  return ( ( forward * S <= 0 ) * ( ( backward + forward ) * S <= 0 ) * forward +
           ( backward * S > 0 ) * ( ( forward + backward ) * S > 0 ) * backward );
}

// 2D Morel-Sussman redistance PDE
BZ_DECLARE_STENCIL_OPERATOR2(morelSussman2D,A,S)
return (  S * 
	         ( 1 - sqrt( pow2( sussmanGrad(A,S,firstDim) ) + 
			             pow2( sussmanGrad(A,S,secondDim) ) ) ) );
BZ_END_STENCIL_OPERATOR

// 3D Morel-Sussman redistance PDE
BZ_DECLARE_STENCIL_OPERATOR2(morelSussman3D,A,S)
	return (  S * 
	         ( 1 - sqrt( pow2( sussmanGrad(A,S,firstDim) ) + 
			             pow2( sussmanGrad(A,S,secondDim) ) +
						 pow2( sussmanGrad(A,S,thirdDim) ) ) ) );
BZ_END_STENCIL_OPERATOR

// 2D Morel-Osher redistance PDE
BZ_DECLARE_STENCIL_OPERATOR2(morelOsher2D,A,S)
return S * ( S >= 0 ) * 
( 1 - sqrt( pow2( blitz::minmax::max( backward11n(A,firstDim), 0 ) ) +
            pow2( blitz::minmax::min(  forward11n(A,firstDim), 0 ) ) + 
            pow2( blitz::minmax::max( backward11n(A,secondDim), 0 ) ) +
            pow2( blitz::minmax::min(  forward11n(A,secondDim), 0 ) ) ) )
+ S * ( S < 0 ) * 
( 1 - sqrt( pow2( blitz::minmax::min( backward11n(A,firstDim), 0 ) ) +
            pow2( blitz::minmax::max(  forward11n(A,firstDim), 0 ) ) + 
            pow2( blitz::minmax::min( backward11n(A,secondDim), 0 ) ) +
            pow2( blitz::minmax::max(  forward11n(A,secondDim), 0 ) ) ) );
BZ_END_STENCIL_OPERATOR

// 3D Morel-Osher redistance PDE
BZ_DECLARE_STENCIL_OPERATOR2(morelOsher3D,A,S)
return S * ( S >= 0 ) * 
( 1 - sqrt( pow2( blitz::minmax::max( backward11n(A,firstDim), 0 ) ) +
            pow2( blitz::minmax::min(  forward11n(A,firstDim), 0 ) ) + 
            pow2( blitz::minmax::max( backward11n(A,secondDim), 0 ) ) +
            pow2( blitz::minmax::min(  forward11n(A,secondDim), 0 ) ) + 
            pow2( blitz::minmax::max( backward11n(A,thirdDim), 0 ) ) +
            pow2( blitz::minmax::min(  forward11n(A,thirdDim), 0 ) ) ) )
+ S * ( S < 0 ) * 
( 1 - sqrt( pow2( blitz::minmax::min( backward11n(A,firstDim), 0 ) ) +
            pow2( blitz::minmax::max(  forward11n(A,firstDim), 0 ) ) + 
            pow2( blitz::minmax::min( backward11n(A,secondDim), 0 ) ) +
            pow2( blitz::minmax::max(  forward11n(A,secondDim), 0 ) ) + 
            pow2( blitz::minmax::min( backward11n(A,thirdDim), 0 ) ) +
            pow2( blitz::minmax::max(  forward11n(A,thirdDim), 0 ) ) ) );
BZ_END_STENCIL_OPERATOR

BZ_ET_STENCIL(nablaPlus2D,P_numtype)
BZ_ET_STENCIL(nablaPlus3D,P_numtype)
BZ_ET_STENCIL(nablaMinus2D,P_numtype)
BZ_ET_STENCIL(nablaMinus3D,P_numtype)
BZ_ET_STENCIL(nablaGodunov2D,P_numtype)
BZ_ET_STENCIL(nablaGodunov3D,P_numtype)

BZ_ET_STENCIL(meanK2D,P_numtype)
BZ_ET_STENCIL(meanK3D,P_numtype)

BZ_ET_STENCIL2(expansion2D,P_numtype)
BZ_ET_STENCIL2(expansion3D,P_numtype)
BZ_ET_STENCIL(curvature2D,P_numtype)
BZ_ET_STENCIL(curvature3D,P_numtype)
BZ_ET_STENCIL3(advection2D,P_numtype)
BZ_ET_STENCIL4(advection3D,P_numtype)

BZ_ET_STENCIL2(diffusion2D,P_numtype)
BZ_ET_STENCIL2(diffusion3D,P_numtype)

BZ_ET_STENCIL2(morelSussman2D,P_numtype)
BZ_ET_STENCIL2(morelSussman3D,P_numtype)
BZ_ET_STENCIL2(morelOsher2D,P_numtype)
BZ_ET_STENCIL2(morelOsher3D,P_numtype)

BZ_ET_STENCIL3(divn,P_numtype)

#endif /* FILE_nbfDifferentials */

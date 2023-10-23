#ifndef FILE_nbfMorphologyFilters
#define FILE_nbfMorphologyFilters

// 3D erosion

template<class T>													
inline _bz_typename T::T_numtype									
erode3D(T& A )													
{					
	int size = 1;
	int res = 0;											
	for ( int i = -size; i <= size; i++ ){								
		for ( int j = -size; j <= size; j++ ){							
			for ( int k = -size; k <= size; k++ ){						
				//if ( A(i,j,k) == true ){
				//	res += 1;
				//}
				res += (int)(A(i,j,k));	
			}														
		}															
	}																
	return res > 9;														
}																	

// allow direct call of stencil operator
BZ_ET_STENCIL(erode3D,P_numtype)

template<class T>													
inline _bz_typename T::T_numtype									
dilate3D(T& A )													
{					
	int size = 1;
	int res = 0;											
	for ( int i = -size; i <= size; i++ ){								
		for ( int j = -size; j <= size; j++ ){							
			for ( int k = -size; k <= size; k++ ){
				res += (int)A(i,j,k);	
			}														
		}															
	}																
	return res > 0;														
}																	

// allow direct call of stencil operator
BZ_ET_STENCIL(dilate3D,P_numtype)


template<class T>													
inline _bz_typename T::T_numtype									
dilate2D(T& A )													
{					
	int size = 1;
	int res = 0;											
	for ( int i = -size; i <= size; i++ ){								
		for ( int j = -size; j <= size; j++ ){							
			res += (int)A(i,j);	
		}														
	}																															
	return res > 0;														
}																	

// allow direct call of stencil operator
BZ_ET_STENCIL(dilate2D,P_numtype)

#endif /* FILE_nbfVeselnessFilter */

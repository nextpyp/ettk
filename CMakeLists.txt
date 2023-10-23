#*=========================================================================
#
#  Program:   ETTK
#  Language:  C++
#  Date:      $Date: $WCDATE$ $
#  Version:   $Revision: $WCREV$
#  Author:    Alberto Bartesaghi
#
#=========================================================================*

CMAKE_MINIMUM_REQUIRED(VERSION 1.8)

PROJECT(ETTK)

# Configure External Packages

FIND_PACKAGE(VTK REQUIRED PATHS external/VTK5.10.1 )
FIND_PACKAGE(ITK REQUIRED PATHS external/InsightToolkit-4.2.1/build )

# Eliminate cmake warning message
IF(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
	cmake_policy(SET CMP0015 NEW)
ENDIF(COMMAND cmake_policy)

SET(CMAKE_BUILD_TYPE Release)

SET(CMAKE_VERBOSE_MAKEFILE FALSE)

# set to true to display the full compile and link commands for debugging
OPTION(CMAKE_VERBOSE_MAKEFILE "set this to true if you are using makefiles and want to see the full compile and link commands." OFF)

# Use mpi c++ compiler
SET(CMAKE_CXX_COMPILER g++)
# Set maximum optimization, eliminate warnings
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -w")
# Set blitz debug option: BZDEBUG
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBZDEBUG -w")

#SET(LIBRARY_OUTPUT_PATH bin)
#SET(RUNTIME_OUTPUT_DIRECTORY bin)
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR_BIN})

# Configure Blitz
FIND_PATH(BLITZ_DIR blitz/array.h ${CMAKE_CURRENT_SOURCE_DIR}/external/blitz-0.9)
SET(BLITZ_INCLUDE_DIR ${BLITZ_DIR})
SET(BLITZ_LIBRARY_DIR ${BLITZ_DIR}/lib)
# SET(BLITZ_LIBRARIES blitz)

# Configure K-means
FIND_PATH(KMEANS_DIRECTORY KCtree.h ${CMAKE_CURRENT_SOURCE_DIR}/external/kmlocal-1.7.1/src)

# Configure Soft
FIND_PATH(SOFT_DIRECTORY soft_fftw.h ${CMAKE_CURRENT_SOURCE_DIR}/external/soft-1.0)

# Configure svdlibc
FIND_PATH(SVDLIBC_DIRECTORY svdlib.h ${CMAKE_CURRENT_SOURCE_DIR}/external/svdlibc-1.34)

# Configure MXML
FIND_PATH(MXML_DIR mxml.h ${CMAKE_CURRENT_SOURCE_DIR}/external/mxml-2.2.2)
SET(MXML_LIBRARY_DIR ${MXML_DIR})
SET(MXML_LIBRARIES mxml)

# Configure MPI
SET(MPI_LIBRARY_DIR /usr/local/mpich2/lib)
SET(MPI_INCLUDE_DIR /usr/local/mpich2/include)
#SET(MPI_LIBRARIES mpichcxx mpich pthread rt)

# Configure FFTW
SET(FFTW_LIBRARY_DIR /usr/lib)
SET(FFTW_INCLUDE_DIR /usr/include)
SET(FFTW_LIBRARIES fftw3f fftw3)

SET(ITK_INCLUDES
${ITK_INCLUDE_DIRS_BUILD_TREE}
${ITK_INCLUDE_DIRS_BUILD_TREE_CXX}
${ITK_INCLUDE_DIRS_SYSTEM}
)

SET(ITK_LIBRARY_PATH
${ITK_DIR}/bin
)

# External include directories
INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}/inc
${BLITZ_INCLUDE_DIR}
${SOFT_DIRECTORY} 
${KMEANS_DIRECTORY} 
${MPI_INCLUDE_DIR}
${FFTW_INCLUDE_DIR}
${VTK_INCLUDE_DIRS}
${ITK_INCLUDE_DIRS}
${MXML_DIR}
${SVDLIBC_DIRECTORY}
)

SET(KMEANS_SOURCES
${KMEANS_DIRECTORY}/KCtree.cpp
${KMEANS_DIRECTORY}/KCutil.cpp
${KMEANS_DIRECTORY}/KM_ANN.cpp
${KMEANS_DIRECTORY}/KMcenters.cpp
${KMEANS_DIRECTORY}/KMdata.cpp
${KMEANS_DIRECTORY}/KMeans.cpp
${KMEANS_DIRECTORY}/KMfilterCenters.cpp
${KMEANS_DIRECTORY}/KMlocal.cpp
${KMEANS_DIRECTORY}/KMrand.cpp
${KMEANS_DIRECTORY}/KMterm.cpp
)

SET(SVDLIBC_SOURCES
${SVDLIBC_DIRECTORY}/las2.c
${SVDLIBC_DIRECTORY}/svdlib.c
${SVDLIBC_DIRECTORY}/svdutil.c
)

SET(SOFT_SOURCES
${SOFT_DIRECTORY}/cospmls.c
${SOFT_DIRECTORY}/csecond.c
${SOFT_DIRECTORY}/fft_grids.c
${SOFT_DIRECTORY}/fft_grids_so3.c
${SOFT_DIRECTORY}/FFTcode.c
${SOFT_DIRECTORY}/FST_semi_memo.c
${SOFT_DIRECTORY}/indextables.c
${SOFT_DIRECTORY}/makeWigner.c
${SOFT_DIRECTORY}/naive_synthesis.c
${SOFT_DIRECTORY}/newFCT.c
${SOFT_DIRECTORY}/oddweights.c
${SOFT_DIRECTORY}/OURmods.c
${SOFT_DIRECTORY}/OURperms.c
${SOFT_DIRECTORY}/permroots.c
${SOFT_DIRECTORY}/primitive.c
${SOFT_DIRECTORY}/primitive_FST.c
${SOFT_DIRECTORY}/rotate_so3.c
${SOFT_DIRECTORY}/rotate_so3_mem.c
${SOFT_DIRECTORY}/seminaive.c
${SOFT_DIRECTORY}/so3_correlate_fftw.c
${SOFT_DIRECTORY}/so3_correlate_sym.c
${SOFT_DIRECTORY}/soft.c
${SOFT_DIRECTORY}/soft_fftw.c
${SOFT_DIRECTORY}/soft_fftw_pc.c
${SOFT_DIRECTORY}/soft_fftw_wo.c
${SOFT_DIRECTORY}/soft_sym.c
${SOFT_DIRECTORY}/utils_so3.c
${SOFT_DIRECTORY}/utils_vec_cx.c
${SOFT_DIRECTORY}/weights.c
${SOFT_DIRECTORY}/wignerTransforms.c
${SOFT_DIRECTORY}/wignerTransforms_fftw.c
${SOFT_DIRECTORY}/wignerTransforms_sym.c
)

LINK_DIRECTORIES(
${MPI_LIBRARY_DIR}
${VTK_LIBRARY_DIRS}
${ITK_LIBRARY_PATH}
${BLITZ_LIBRARY_DIR}
${FFTW_LIBRARY_DIR}
${MXML_LIBRARY_DIR}
${CMAKE_CURRENT_SOURCE_DIR_BIN}
)

SET(ITK_LIBRARIES
ITKCommon
)
# ITKIO

SET(VTK_LIBRARIES
vtkCommon vtkImaging vtkRendering vtkFiltering
)

LINK_LIBRARIES(
${MPI_LIBRARIES}
${VTK_LIBRARIES}
${ITK_LIBRARIES}
${BLITZ_LIBRARIES}
${FFTW_LIBRARIES}
${MXML_LIBRARIES}
NBFObjects
)

# Compile all common objects in single library
ADD_LIBRARY(NBFObjects inc/io/nbfVTKInterface.cpp ${SOFT_SOURCES} ${SVDLIBC_SOURCES} ${KMEANS_SOURCES})

SET(EXECUTABLE_OUTPUT_PATH
${CMAKE_CURRENT_SOURCE_DIR}/bin
)

# Generate rules for all executables
FOREACH(entry Correlation3DNew
			  CutVolumes3DFromPositions
			  LoopCreateVolumeList
			  #MPI_Classification
			  Test_Metric_Filter
			  virus_segment_membrane
			  virus_segment_membrane_select_threshold)
    ADD_EXECUTABLE(${entry} src/${entry}.cpp)
ENDFOREACH(entry)

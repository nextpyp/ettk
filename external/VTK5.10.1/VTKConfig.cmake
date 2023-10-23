#-----------------------------------------------------------------------------
#
# VTKConfig.cmake - VTK CMake configuration file for external projects.
#
# This file is configured by VTK and used by the UseVTK.cmake module
# to load VTK's settings for an external project.


# The set of VTK libraries
SET(VTK_LIBRARIES "vtkCommon;vtkFiltering;vtkImaging;vtkGraphics;vtkGenericFiltering;vtkIO;vtkRendering;vtkVolumeRendering;vtkHybrid;vtkWidgets;vtkInfovis;vtkGeovis;vtkViews;vtkCharts")

# The VTK include file directories.
SET(VTK_INCLUDE_DIRS "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/VolumeRendering;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Rendering;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Charts;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Chemistry;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkalglib;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Infovis;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Geovis;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Views;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/VolumeRendering;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Hybrid;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Widgets;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Rendering;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Charts;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Chemistry;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Rendering/Testing/Cxx;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/IO;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Imaging;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Graphics;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/GenericFiltering;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Filtering;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/Cxx;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtknetcdf;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtknetcdf;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtknetcdf/include;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtklibproj4;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtklibproj4;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/DICOMParser;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/DICOMParser;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkfreetype/include;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkfreetype/include;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/LSDyna;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/LSDyna;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/MaterialLibrary;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkmetaio;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkmetaio;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/verdict;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/verdict;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5/src;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5/src;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5/hl/src;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5/hl/src;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/utf8/source;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Infovis;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkalglib;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Geovis;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Views")

# The VTK library directories.
SET(VTK_LIBRARY_DIRS "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin")

# The VTK binary executable directories.  Note that if
# VTK_CONFIGURATION_TYPES is set (see below) then these directories
# will be the parent directories under which there will be a directory
# of runtime binaries for each configuration type.
SET(VTK_EXECUTABLE_DIRS "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin")

# The VTK runtime library directories.  Note that if
# VTK_CONFIGURATION_TYPES is set (see below) then these directories
# will be the parent directories under which there will be a directory
# of runtime libraries for each configuration type.
SET(VTK_RUNTIME_LIBRARY_DIRS "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin")

# The runtime library path variable name e.g. LD_LIBRARY_PATH,
# this environment variable should be set to VTK_RUNTIME_LIBRARY_DIRS
SET(VTK_RUNTIME_PATH_VAR_NAME "LD_LIBRARY_PATH")

# The C and C++ flags added by VTK to the cmake-configured flags.
SET(VTK_REQUIRED_C_FLAGS "")
SET(VTK_REQUIRED_CXX_FLAGS " -Wno-deprecated")
SET(VTK_REQUIRED_EXE_LINKER_FLAGS "")
SET(VTK_REQUIRED_SHARED_LINKER_FLAGS "")
SET(VTK_REQUIRED_MODULE_LINKER_FLAGS "")

# The VTK version number
SET(VTK_MAJOR_VERSION "5")
SET(VTK_MINOR_VERSION "10")
SET(VTK_BUILD_VERSION "1")

# The location of the UseVTK.cmake file.
SET(VTK_USE_FILE "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/UseVTK.cmake")

# The build settings file.
SET(VTK_BUILD_SETTINGS_FILE "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/VTKBuildSettings.cmake")

# The directory containing class list files for each kit.
SET(VTK_KITS_DIR "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities")

# The wrapping hints file.
SET(VTK_WRAP_HINTS "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Wrapping/hints")

# CMake extension module directory and macro file.
SET(VTK_LOAD_CMAKE_EXTENSIONS_MACRO "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/CMake/vtkLoadCMakeExtensions.cmake")
SET(VTK_CMAKE_DIR "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/CMake")
SET(VTK_CMAKE_EXTENSIONS_DIR "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/CMake")

# The list of available kits.
SET(VTK_KITS "COMMON;FILTERING;IO;GRAPHICS;GENERIC_FILTERING;IMAGING;RENDERING;VOLUMERENDERING;HYBRID;WIDGETS;INFOVIS;GEOVIS;VIEWS;CHARTS")

# The list of available languages.
SET(VTK_LANGUAGES "")

# VTK Configuration options.
SET(VTK_BUILD_SHARED_LIBS "OFF")
SET(VTK_DEBUG_LEAKS "OFF")
SET(VTK_HAVE_VP1000 "")
SET(VTK_LEGACY_REMOVE "OFF")
SET(VTK_LEGACY_SILENT "OFF")
SET(VTK_OPENGL_HAS_OSMESA "OFF")
SET(VTK_USE_64BIT_IDS "ON")
SET(VTK_USE_ANSI_STDLIB "ON")
SET(VTK_USE_BOOST "OFF")
SET(VTK_USE_QT "OFF")
SET(VTK_USE_CARBON "OFF")
SET(VTK_USE_CG_SHADERS "OFF")
SET(VTK_USE_CHARTS "ON")
SET(VTK_USE_CHEMISTRY "OFF")
SET(VTK_USE_COCOA "OFF")
SET(VTK_USE_GEOVIS "ON")
SET(VTK_USE_GL2PS "OFF")
SET(VTK_USE_GLSL_SHADERS "ON")
SET(VTK_USE_GNU_R "OFF")
SET(VTK_USE_GUISUPPORT "OFF")
SET(VTK_USE_INFOVIS "ON")
SET(VTK_USE_MANGLED_MESA "OFF")
SET(VTK_USE_MATLAB_MEX "OFF")
SET(VTK_USE_MATROX_IMAGING "OFF")
SET(VTK_USE_METAIO "ON")
SET(VTK_USE_MFC "")
SET(VTK_USE_MPI "OFF")
SET(VTK_USE_MYSQL "OFF")
SET(VTK_USE_N_WAY_ARRAYS "ON")
SET(VTK_USE_ODBC "OFF")
SET(VTK_USE_PARALLEL "OFF")
SET(VTK_USE_PARALLEL_BGL "OFF")
SET(VTK_USE_POSTGRES "OFF") 
SET(VTK_USE_QVTK "OFF")
SET(VTK_USE_QVTK_OPENGL "")
SET(VTK_USE_RENDERING "ON")
SET(VTK_USE_TDX "OFF")
SET(VTK_USE_TEXT_ANALYSIS "OFF")
SET(VTK_USE_TK "OFF")
SET(VTK_USE_TK "OFF")
SET(VTK_USE_VIDEO_FOR_WINDOWS "")
SET(VTK_USE_VIEWS "ON")
SET(VTK_USE_VOLUMEPRO_1000 "OFF")
SET(VTK_USE_X "ON")
SET(VTK_WRAP_JAVA "OFF")
SET(VTK_WRAP_PYTHON "OFF")
SET(VTK_WRAP_TCL "OFF")
SET(VTK_WRAP_PYTHON_SIP "OFF")

# The Hybrid and VolumeRendering kits are now switched with Rendering.
SET(VTK_USE_HYBRID "ON")
SET(VTK_USE_VOLUMERENDERING "ON")

# The MPI configuration
SET(VTK_MPIRUN_EXE "")
SET(VTK_MPI_CLIENT_POSTFLAGS "")
SET(VTK_MPI_CLIENT_PREFLAGS "")
SET(VTK_MPI_MAX_NUMPROCS "")
SET(VTK_MPI_PRENUMPROC_FLAGS "")
SET(VTK_MPI_NUMPROC_FLAG "")
SET(VTK_MPI_POSTFLAGS "")
SET(VTK_MPI_PREFLAGS "")
SET(VTK_MPI_SERVER_POSTFLAGS "")
SET(VTK_MPI_SERVER_PREFLAGS "")
SET(VTK_MPI_INCLUDE_DIR "")
SET(VTK_MPI_LIBRARIES ";")

# The Tcl/Tk configuration.
SET(VTK_TCL_TK_STATIC "")
SET(VTK_TCL_TK_COPY_SUPPORT_LIBRARY "")
SET(VTK_TCL_SUPPORT_LIBRARY_PATH "")
SET(VTK_TK_SUPPORT_LIBRARY_PATH "")
SET(VTK_TCL_TK_MACROS_MODULE "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/CMake/vtkTclTkMacros.cmake")
SET(VTK_TCL_HOME "")
IF (NOT CMAKE_CROSSCOMPILING)
  SET(VTK_WRAP_TCL_EXE "")
  SET(VTK_WRAP_TCL_INIT_EXE "")
ENDIF (NOT CMAKE_CROSSCOMPILING)
SET(VTK_TK_INTERNAL_DIR "")
SET(VTK_TK_RESOURCES_DIR "")
SET(VTK_TCL_INCLUDE_DIR "")
SET(VTK_TCL_LIBRARY "")
SET(VTK_TK_INCLUDE_DIR "")
SET(VTK_TK_LIBRARY "")

# The Java configuration.
SET(VTK_JAVA_JAR "")
IF (NOT CMAKE_CROSSCOMPILING)
  SET(VTK_PARSE_JAVA_EXE "")
  SET(VTK_WRAP_JAVA_EXE "")
ENDIF (NOT CMAKE_CROSSCOMPILING)
SET(VTK_JAVA_INCLUDE_DIR ";")
SET(VTK_JAVA_AWT_LIBRARY "")
SET(VTK_JVM_LIBRARY "")

# The Matlab configuration.
SET(VTK_MATLAB_ROOT_DIR "")
SET(VTK_MATLAB_INCLUDE_DIR "")
SET(VTK_MATLAB_LIB_DIR "")

# The Python configuration.
# If VTK_CONFIGURATION_TYPES is set (see below) then the VTK_PYTHONPATH_DIRS
# will have subdirectories for each configuration type.
SET(VTK_PYTHONPATH_DIRS "")
IF (NOT CMAKE_CROSSCOMPILING)
  SET(VTK_WRAP_PYTHON_EXE "")
  SET(VTK_WRAP_PYTHON_INIT_EXE "")
ENDIF (NOT CMAKE_CROSSCOMPILING)
SET(VTK_PYTHON_INCLUDE_DIR "")
SET(VTK_PYTHON_LIBRARY "")

# Other executables
IF (NOT CMAKE_CROSSCOMPILING)
  SET(VTK_ENCODESTRING_EXE "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/vtkEncodeString")
ENDIF (NOT CMAKE_CROSSCOMPILING)

# The Doxygen configuration.
SET(VTK_DOXYGEN_HOME "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/Doxygen")

# The VTK test script locations.
SET(VTK_HEADER_TESTING_PY "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py")
SET(VTK_FIND_STRING_TCL "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/Tcl/FindString.tcl")
SET(VTK_PRINT_SELF_CHECK_TCL "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/Tcl/PrintSelfCheck.tcl")
SET(VTK_RT_IMAGE_TEST_TCL "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/Tcl/rtImageTest.tcl")
SET(VTK_PRT_IMAGE_TEST_TCL "")

# The names of utility libraries used by VTK.
SET(VTK_PNG_LIBRARIES        "vtkpng")
SET(VTK_PNG__INCLUDE_DIR     "")
SET(VTK_ZLIB_LIBRARIES       "vtkzlib")
SET(VTK_ZLIB_INCLUDE_DIR     "")
SET(VTK_JPEG_LIBRARIES       "vtkjpeg")
SET(VTK_JPEG_INCLUDE_DIR     "")
SET(VTK_TIFF_LIBRARIES       "vtktiff")
SET(VTK_TIFF_INCLUDE_DIR     "")
SET(VTK_EXPAT_LIBRARIES      "vtkexpat")
SET(VTK_EXPAT_INCLUDE_DIR    "")
SET(VTK_FREETYPE_LIBRARIES   "vtkfreetype")
SET(VTK_FREETYPE_INCLUDE_DIR "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkfreetype/include;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkfreetype/include")
SET(VTK_LIBXML2_LIBRARIES    "vtklibxml2")
SET(VTK_LIBXML2_INCLUDE_DIR  "")
SET(VTK_LIBPROJ4_LIBRARIES   "vtkproj4")
SET(VTK_LIBPROJ4_INCLUDE_DIR "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtklibproj4;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtklibproj4")
SET(VTK_HDF5_LIBRARIES       "vtkhdf5;vtkhdf5_hl")
SET(VTK_HDF5_INCLUDE_DIR     "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5/hl/src;/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Utilities/vtkhdf5/hl/src")

# The VTK Qt configuration.
IF(VTK_USE_QVTK)
  INCLUDE(${VTK_DIR}/VTKConfigQt.cmake)
ENDIF(VTK_USE_QVTK)

# Relative install paths in the VTK install tree
SET(VTK_INSTALL_BIN_DIR "/bin")
SET(VTK_INSTALL_INCLUDE_DIR "/include/vtk-5.10")
SET(VTK_INSTALL_LIB_DIR "/lib/vtk-5.10")
SET(VTK_INSTALL_PACKAGE_DIR "/lib/vtk-5.10")

# A VTK install tree always provides one build configuration.  A VTK
# build tree may provide either one or multiple build configurations
# depending on the CMake generator used.  Since VTK can be used either
# from a build tree or an install tree it is useful for outside
# projects to know the configurations available.  If this
# VTKConfig.cmake is in a VTK install tree VTK_CONFIGURATION_TYPES
# will be empty and VTK_BUILD_TYPE will be set to the value of
# CMAKE_BUILD_TYPE used to build VTK.  If VTKConfig.cmake is in a VTK
# build tree then VTK_CONFIGURATION_TYPES and VTK_BUILD_TYPE will have
# values matching CMAKE_CONFIGURATION_TYPES and CMAKE_BUILD_TYPE for
# that build tree (only one will ever be set).
SET(VTK_CONFIGURATION_TYPES )
SET(VTK_BUILD_TYPE Debug)

# The VTK targets file.
IF(NOT VTK_INSTALL_EXPORT_NAME AND NOT TARGET vtkCommon
   AND EXISTS "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/VTKTargets.cmake")
  INCLUDE("/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/VTKTargets.cmake")
ENDIF()

# The old, less clear name for VTK_RUNTIME_LIBRARY_DIRS.  Kept here
# for compatibility.
SET(VTK_RUNTIME_DIRS ${VTK_RUNTIME_LIBRARY_DIRS})

# The name of the encoders used by VTK
SET( VTK_USE_FFMPEG_ENCODER "OFF" )
SET( VTK_USE_OGGTHEORA_ENCODER "OFF")


# For backward compatability.  DO NOT USE.
SET(VTK_SOURCE_DIR "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1")
IF(NOT TCL_LIBRARY)
  SET(TCL_LIBRARY "" CACHE FILEPATH "Location of Tcl library imported from VTK.  This may mean your project is depending on VTK to get this setting.  Consider using FindTCL.cmake.")
ENDIF(NOT TCL_LIBRARY)
IF(NOT TK_LIBRARY)
  SET(TK_LIBRARY "" CACHE FILEPATH "Location of Tk library imported from VTK.  This may mean your project is depending on VTK to get this setting.  Consider using FindTCL.cmake.")
ENDIF(NOT TK_LIBRARY)
MARK_AS_ADVANCED(TCL_LIBRARY TK_LIBRARY)


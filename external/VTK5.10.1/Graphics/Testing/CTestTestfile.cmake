# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Graphics/Testing
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Graphics/Testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(HeaderTesting-Graphics "/opt/apps/rhel7/anaconda2/bin/python" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Graphics" "VTK_GRAPHICS_EXPORT" "vtkButterflySubdivisionFilter.h" "vtkDijkstraGraphInternals.h" "vtkLinearSubdivisionFilter.h" "vtkLoopSubdivisionFilter.h" "vtkOutlineFilter.h" "vtkProgrammableFilter.h" "vtkProgrammableSource.h" "vtkStructuredGridOutlineFilter.h" "vtkStructuredPointsGeometryFilter.h" "vtkTableBasedClipCases.h")
subdirs("Cxx")

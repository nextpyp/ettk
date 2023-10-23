# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Geovis/Testing
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Geovis/Testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(HeaderTesting-Geovis "/opt/apps/rhel7/anaconda2/bin/python" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Geovis" "VTK_GEOVIS_EXPORT" "vtkGeoPatch.h")
subdirs("Cxx")

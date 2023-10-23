# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/GenericFiltering/Testing
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/GenericFiltering/Testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(HeaderTesting-GenericFiltering "/opt/apps/rhel7/anaconda2/bin/python" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/GenericFiltering" "VTK_GENERIC_FILTERING_EXPORT")
subdirs("Cxx")

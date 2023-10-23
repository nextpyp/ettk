# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/IO/Testing
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/IO/Testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(HeaderTesting-IO "/opt/apps/rhel7/anaconda2/bin/python" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/IO" "VTK_IO_EXPORT" "vtkBase64Utilities.h" "vtkMINC.h" "vtkMySQLDatabasePrivate.h" "vtkODBCInternals.h" "vtkOffsetsManagerArray.h" "vtkPLY.h" "vtkPostgreSQLDatabasePrivate.h" "vtkXMLUtilities.h" "vtkXMLWriterC.h" "vtkXMLWriterF.h")
subdirs("Cxx")

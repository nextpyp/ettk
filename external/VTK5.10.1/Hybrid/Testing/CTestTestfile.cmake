# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Hybrid/Testing
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Hybrid/Testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(HeaderTesting-Hybrid "/opt/apps/rhel7/anaconda2/bin/python" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Hybrid" "VTK_HYBRID_EXPORT" "vtk3DS.h" "vtkExodusIIReaderParser.h" "vtkExodusIIReaderPrivate.h" "vtkExodusIIReaderVariableCheck.h" "vtkVRML.h" "vtkX3D.h" "vtkX3DExporterFIWriterHelper.h")
subdirs("Cxx")

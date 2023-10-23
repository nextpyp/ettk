# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Infovis/Testing
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Infovis/Testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(HeaderTesting-Infovis "/opt/apps/rhel7/anaconda2/bin/python" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Infovis" "VTK_INFOVIS_EXPORT" "vtkArrayNorm.h" "vtkMultiCorrelativeStatisticsAssessFunctor.h" "vtkKMeansAssessFunctor.h" "vtkBoostGraphAdapter.h" "vtkBoostRandomSparseArraySource.h" "vtkSQLDatabaseGraphSource.h" "vtkSQLDatabaseTableSource.h" "vtkStatisticsAlgorithmPrivate.h")
subdirs("Cxx")

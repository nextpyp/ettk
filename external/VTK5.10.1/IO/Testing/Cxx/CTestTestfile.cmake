# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/IO/Testing/Cxx
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/IO/Testing/Cxx
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TestSimplePointsReaderWriter "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/IOCxxTests" "TestSimplePointsReaderWriter")
add_test(TestSQLDatabaseSchema "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/IOCxxTests" "TestSQLDatabaseSchema")
add_test(TestSQLiteDatabase "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/IOCxxTests" "TestSQLiteDatabase")
set_tests_properties(TestSQLiteDatabase PROPERTIES  RUN_SERIAL "ON")
add_test(TestDataObjectIO "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/IOCxxTests" "TestDataObjectIO")
add_test(TestDataObjectXMLIO "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/IOCxxTests" "TestDataObjectXMLIO")
set_tests_properties(TestDataObjectXMLIO PROPERTIES  RUN_SERIAL "ON")
add_test(Array-TestArraySerialization "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/ArrayIOCxxTests" "TestArraySerialization")
add_test(Array-TestArrayDataWriter "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/ArrayIOCxxTests" "TestArrayDataWriter")
add_test(Array-TestArrayDenormalized "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/bin/ArrayIOCxxTests" "TestArrayDenormalized")

itk_module_test()
set(ITKIOMetaTests
itkMetaImageIOMetaDataTest.cxx
itkMetaImageIOGzTest.cxx
itkMetaImageIOTest.cxx
itkMetaImageIOTest2.cxx
itkLargeMetaImageWriteReadTest.cxx
testMetaArray.cxx
testMetaBlob.cxx
testMetaCommand.cxx
testMetaGroup.cxx
testMetaImage.cxx
testMetaLandmark.cxx
testMetaLine.cxx
testMetaMesh.cxx
testMetaObject.cxx
testMetaScene.cxx
testMetaSurface.cxx
testMetaTube.cxx
testMetaUtils.cxx
itkMetaImageStreamingIOTest.cxx
itkMetaImageStreamingWriterIOTest.cxx
)

CreateTestDriver(ITKIOMeta  "${ITKIOMeta-Test_LIBRARIES}" "${ITKIOMetaTests}")

itk_add_test(NAME itkMetaImageIOMetaDataTest
      COMMAND ITKIOMetaTestDriver itkMetaImageIOMetaDataTest
              ${ITK_TEST_OUTPUT_DIR}/MetaImageIOMetaDataTest.mhd)
itk_add_test(NAME itkMetaImageIOGzTest
      COMMAND ITKIOMetaTestDriver itkMetaImageIOGzTest
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME itkMetaImageIOTest
      COMMAND ITKIOMetaTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/HeadMRVolume.mhd,HeadMRVolume.raw}
              ${ITK_TEST_OUTPUT_DIR}/HeadMRVolume.mhd
    itkMetaImageIOTest DATA{${ITK_DATA_ROOT}/Input/HeadMRVolume.mhd,HeadMRVolume.raw} ${ITK_TEST_OUTPUT_DIR}/HeadMRVolume.mhd)
itk_add_test(NAME itkMetaImageIOTest2
      COMMAND ITKIOMetaTestDriver itkMetaImageIOTest2
      ${ITK_TEST_OUTPUT_DIR}/itkMetaImageIOTest2.mha)
itk_add_test(NAME itkMetaImageIOShouldFailTest
      COMMAND ITKIOMetaTestDriver itkMetaImageIOTest
              DATA{${ITK_DATA_ROOT}/Input/MetaImageError.mhd} 1)
itk_add_test(NAME testMetaGroup
      COMMAND ITKIOMetaTestDriver testMetaGroup
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaImage
      COMMAND ITKIOMetaTestDriver testMetaImage)
itk_add_test(NAME testMetaLandmark
      COMMAND ITKIOMetaTestDriver testMetaLandmark
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaLine
      COMMAND ITKIOMetaTestDriver testMetaLine
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaMesh
      COMMAND ITKIOMetaTestDriver testMetaMesh
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaObject
      COMMAND ITKIOMetaTestDriver testMetaObject
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaScene
      COMMAND ITKIOMetaTestDriver testMetaScene
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaSurface
      COMMAND ITKIOMetaTestDriver testMetaSurface
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaTube
      COMMAND ITKIOMetaTestDriver testMetaTube
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaUtils
      COMMAND ITKIOMetaTestDriver testMetaUtils
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME MetaCommandRequiredOptionsWithMultipleStrings
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 7 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7)
itk_add_test(NAME MetaCommandRequiredOptionsWithMultipleStringsAtEnd
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 7 -r2_ri 7 -r1c_rs StringValue1 StringValue2 StringValue3)
itk_add_test(NAME MetaCommandOptionalWithRequired
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 107 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7 -o1_ri 1 -o2_ri 99)
itk_add_test(NAME MetaCommandOptionalOptoinsWithOutArguments
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 207 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7 -o1_oi -o2_oi)
itk_add_test(NAME MetaCommandOptionalOptoinsWithOutArguments1
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 107 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7 -o1_oi 1 -o2_oi 99)
itk_add_test(NAME MetaCommandRequiredAndOptionalOptoinsWithArguments
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 207 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7 -o1_ri 1 -o2_ri 99 -o3_ri_oi 50 50)
itk_add_test(NAME MetaCommandRequiredAndOptionalOptoinsWithOutArguments
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 257 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7 -o1_ri 1 -o2_ri 99 -o3_ri_oi 50)
itk_add_test(NAME MetaCommandRequiredAndOptionalOptoinsWithOutArguments-1
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 0 -sum 457 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7 -o1_ri 1 -o2_ri 99 -o3_ri_oi 50 -o1_oi -o2_oi)
itk_add_test(NAME MetaCommandRequiredOptionsWithMissingMultipleStrings
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 1 -sum 7 -r1c_rs -r2_ri 7)
itk_add_test(NAME MetaCommandRequiredOptionsWithMissingMultipleStringsAtEnd
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 1 -sum 7 -r2_ri 7 -r1c_rs)
itk_add_test(NAME MetaCommandOptionalWithRequiredMissing
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 1 -sum 107 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri -o1_ri 1 -o2_ri 99)
itk_add_test(NAME MetaCommandRequiredAndOptionalOptoinsWithArgumentsMissing
      COMMAND ITKIOMetaTestDriver testMetaCommand
              -ExpectedFailStatus 1 -sum 207 -r1c_rs StringValue1 StringValue2 StringValue3 -r2_ri 7 -o3_ri_oi)
itk_add_test(NAME testMetaArray
      COMMAND ITKIOMetaTestDriver testMetaArray
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME testMetaBlob
      COMMAND ITKIOMetaTestDriver testMetaBlob
              ${ITK_TEST_OUTPUT_DIR})
itk_add_test(NAME itkMetaImageStreamingIOTest
      COMMAND ITKIOMetaTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/HeadMRVolume.mhd,HeadMRVolume.raw}
              ${ITK_TEST_OUTPUT_DIR}/HeadMRVolumeStreamed.mhd
    itkMetaImageStreamingIOTest DATA{${ITK_DATA_ROOT}/Input/HeadMRVolume.mhd,HeadMRVolume.raw} ${ITK_TEST_OUTPUT_DIR}/HeadMRVolumeStreamed.mhd)
itk_add_test(NAME itkMetaImageCompressedStreamingIOTest
      COMMAND ITKIOMetaTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/HeadMRVolume.mhd,HeadMRVolume.raw}
              ${ITK_TEST_OUTPUT_DIR}/HeadMRVolumeCompressedStreamed.mha
    itkMetaImageStreamingIOTest DATA{${ITK_DATA_ROOT}/Input/HeadMRVolumeCompressed.mha} ${ITK_TEST_OUTPUT_DIR}/HeadMRVolumeCompressedStreamed.mha)
itk_add_test(NAME itkMetaImageStreamingWriterIOTest
      COMMAND ITKIOMetaTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/HeadMRVolume.mhd,HeadMRVolume.raw}
              ${ITK_TEST_OUTPUT_DIR}/MetaImageStreamingWriterIOTest.mha
    itkMetaImageStreamingWriterIOTest DATA{${ITK_DATA_ROOT}/Input/HeadMRVolume.mhd,HeadMRVolume.raw} ${ITK_TEST_OUTPUT_DIR}/MetaImageStreamingWriterIOTest.mha)
itk_add_test(NAME itkMetaImageStreamingWriterIOTest2
      COMMAND ITKIOMetaTestDriver
    --compare ${ITK_DATA_ROOT}/Input/mri3D.mhd
              ${ITK_TEST_OUTPUT_DIR}/mri3DWriteStreamed.mha
    itkMetaImageStreamingWriterIOTest ${ITK_DATA_ROOT}/Input/mri3D.mhd ${ITK_TEST_OUTPUT_DIR}/mri3DWriteStreamed.mha)

if( "${ITK_COMPUTER_MEMORY_SIZE}" GREATER 5 )

  # Image of 1.67 Gigabytes (pixels size is 16bits)
  itk_add_test(NAME itkLargeMetaImageWriteReadTest1
      COMMAND ITKIOMetaTestDriver
    itkLargeMetaImageWriteReadTest ${ITK_TEST_OUTPUT_DIR}/LargeImage01.mhd  30000L )

  # Image of 2.98 Gigabytes (pixels size is 16bits)
  itk_add_test(NAME itkLargeMetaImageWriteReadTest2
      COMMAND ITKIOMetaTestDriver
    itkLargeMetaImageWriteReadTest ${ITK_TEST_OUTPUT_DIR}/LargeImage02.mhd  40000L )

  # Image of 4.9 Gigabytes (pixels size is 8-bits out and 16-bit in)
  itk_add_test(NAME itkLargeMetaImageWriteReadTest3
      COMMAND ITKIOMetaTestDriver
    itkLargeMetaImageWriteReadTest ${ITK_TEST_OUTPUT_DIR}/LargeImage03.mhd  50000L )

  # Due to the large memory requirements this tests must be run one by one
  set_tests_properties(itkLargeMetaImageWriteReadTest1
                       itkLargeMetaImageWriteReadTest2
                       itkLargeMetaImageWriteReadTest3
                       PROPERTIES
                       RUN_SERIAL 1
                      )
endif()


if( "${ITK_COMPUTER_MEMORY_SIZE}" GREATER 12 )

  # Image of 9.12 Gigabytes (pixels size is 16bits)
  itk_add_test(NAME itkLargeMetaImageWriteReadTest4
      COMMAND ITKIOMetaTestDriver
    itkLargeMetaImageWriteReadTest ${ITK_TEST_OUTPUT_DIR}/LargeImage04.mhd  70000L )

  # Due to the large memory requirements this tests must be run one by one
  set_tests_properties(itkLargeMetaImageWriteReadTest4
    PROPERTIES
    RUN_SERIAL 1
    )

endif()

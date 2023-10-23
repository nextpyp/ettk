# CMake generated Testfile for 
# Source directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Imaging/Testing
# Build directory: /dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Imaging/Testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(HeaderTesting-Imaging "/opt/apps/rhel7/anaconda2/bin/python" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Common/Testing/HeaderTesting.py" "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Imaging" "VTK_IMAGING_EXPORT" "vtkImageAppendComponents.h" "vtkImageCityBlockDistance.h" "vtkImageDivergence.h" "vtkImageDotProduct.h" "vtkImageFFT.h" "vtkImageFourierCenter.h" "vtkImageFourierFilter.h" "vtkImageHybridMedian2D.h" "vtkImageLuminance.h" "vtkImageMagnitude.h" "vtkImageMapToRGBA.h" "vtkImageMirrorPad.h" "vtkImageNormalize.h" "vtkImageRFFT.h" "vtkImageStencilIterator.h" "vtkImageWrapPad.h" "vtkSimpleImageFilterExample.h")
subdirs("Cxx")

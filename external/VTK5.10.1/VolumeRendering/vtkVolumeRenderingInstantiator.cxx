#include "vtkVolumeRenderingInstantiator.h"
  
extern vtkObject* vtkInstantiatorvtkEncodedGradientShaderNew();
extern vtkObject* vtkInstantiatorvtkFiniteDifferenceGradientEstimatorNew();
extern vtkObject* vtkInstantiatorvtkFixedPointRayCastImageNew();
extern vtkObject* vtkInstantiatorvtkFixedPointVolumeRayCastCompositeGOHelperNew();
extern vtkObject* vtkInstantiatorvtkFixedPointVolumeRayCastCompositeGOShadeHelperNew();
extern vtkObject* vtkInstantiatorvtkFixedPointVolumeRayCastCompositeHelperNew();
extern vtkObject* vtkInstantiatorvtkFixedPointVolumeRayCastCompositeShadeHelperNew();
extern vtkObject* vtkInstantiatorvtkFixedPointVolumeRayCastMIPHelperNew();
extern vtkObject* vtkInstantiatorvtkFixedPointVolumeRayCastMapperNew();
extern vtkObject* vtkInstantiatorvtkVolumeRayCastSpaceLeapingImageFilterNew();
extern vtkObject* vtkInstantiatorvtkGPUVolumeRayCastMapperNew();
extern vtkObject* vtkInstantiatorvtkHAVSVolumeMapperNew();
extern vtkObject* vtkInstantiatorvtkProjectedAAHexahedraMapperNew();
extern vtkObject* vtkInstantiatorvtkProjectedTetrahedraMapperNew();
extern vtkObject* vtkInstantiatorvtkRecursiveSphereDirectionEncoderNew();
extern vtkObject* vtkInstantiatorvtkSmartVolumeMapperNew();
extern vtkObject* vtkInstantiatorvtkSphericalDirectionEncoderNew();
extern vtkObject* vtkInstantiatorvtkVolumeOutlineSourceNew();
extern vtkObject* vtkInstantiatorvtkVolumePickerNew();
extern vtkObject* vtkInstantiatorvtkVolumeProMapperNew();
extern vtkObject* vtkInstantiatorvtkVolumeRayCastCompositeFunctionNew();
extern vtkObject* vtkInstantiatorvtkVolumeRayCastIsosurfaceFunctionNew();
extern vtkObject* vtkInstantiatorvtkVolumeRayCastMIPFunctionNew();
extern vtkObject* vtkInstantiatorvtkVolumeRayCastMapperNew();
extern vtkObject* vtkInstantiatorvtkVolumeRenderingFactoryNew();
extern vtkObject* vtkInstantiatorvtkVolumeTextureMapper2DNew();
extern vtkObject* vtkInstantiatorvtkVolumeTextureMapper3DNew();
extern vtkObject* vtkInstantiatorvtkUnstructuredGridBunykRayCastFunctionNew();
extern vtkObject* vtkInstantiatorvtkUnstructuredGridHomogeneousRayIntegratorNew();
extern vtkObject* vtkInstantiatorvtkUnstructuredGridLinearRayIntegratorNew();
extern vtkObject* vtkInstantiatorvtkUnstructuredGridPartialPreIntegrationNew();
extern vtkObject* vtkInstantiatorvtkUnstructuredGridPreIntegrationNew();
extern vtkObject* vtkInstantiatorvtkUnstructuredGridVolumeRayCastMapperNew();
extern vtkObject* vtkInstantiatorvtkUnstructuredGridVolumeZSweepMapperNew();
extern vtkObject* vtkInstantiatorvtkOpenGLGPUVolumeRayCastMapperNew();
extern vtkObject* vtkInstantiatorvtkOpenGLHAVSVolumeMapperNew();
extern vtkObject* vtkInstantiatorvtkOpenGLProjectedAAHexahedraMapperNew();
extern vtkObject* vtkInstantiatorvtkOpenGLProjectedTetrahedraMapperNew();
extern vtkObject* vtkInstantiatorvtkOpenGLRayCastImageDisplayHelperNew();
extern vtkObject* vtkInstantiatorvtkOpenGLVolumeTextureMapper2DNew();
extern vtkObject* vtkInstantiatorvtkOpenGLVolumeTextureMapper3DNew();


  
void vtkVolumeRenderingInstantiator::ClassInitialize()
{
  
  vtkInstantiator::RegisterInstantiator("vtkEncodedGradientShader", vtkInstantiatorvtkEncodedGradientShaderNew);
  vtkInstantiator::RegisterInstantiator("vtkFiniteDifferenceGradientEstimator", vtkInstantiatorvtkFiniteDifferenceGradientEstimatorNew);
  vtkInstantiator::RegisterInstantiator("vtkFixedPointRayCastImage", vtkInstantiatorvtkFixedPointRayCastImageNew);
  vtkInstantiator::RegisterInstantiator("vtkFixedPointVolumeRayCastCompositeGOHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeGOHelperNew);
  vtkInstantiator::RegisterInstantiator("vtkFixedPointVolumeRayCastCompositeGOShadeHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeGOShadeHelperNew);
  vtkInstantiator::RegisterInstantiator("vtkFixedPointVolumeRayCastCompositeHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeHelperNew);
  vtkInstantiator::RegisterInstantiator("vtkFixedPointVolumeRayCastCompositeShadeHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeShadeHelperNew);
  vtkInstantiator::RegisterInstantiator("vtkFixedPointVolumeRayCastMIPHelper", vtkInstantiatorvtkFixedPointVolumeRayCastMIPHelperNew);
  vtkInstantiator::RegisterInstantiator("vtkFixedPointVolumeRayCastMapper", vtkInstantiatorvtkFixedPointVolumeRayCastMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeRayCastSpaceLeapingImageFilter", vtkInstantiatorvtkVolumeRayCastSpaceLeapingImageFilterNew);
  vtkInstantiator::RegisterInstantiator("vtkGPUVolumeRayCastMapper", vtkInstantiatorvtkGPUVolumeRayCastMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkHAVSVolumeMapper", vtkInstantiatorvtkHAVSVolumeMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkProjectedAAHexahedraMapper", vtkInstantiatorvtkProjectedAAHexahedraMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkProjectedTetrahedraMapper", vtkInstantiatorvtkProjectedTetrahedraMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkRecursiveSphereDirectionEncoder", vtkInstantiatorvtkRecursiveSphereDirectionEncoderNew);
  vtkInstantiator::RegisterInstantiator("vtkSmartVolumeMapper", vtkInstantiatorvtkSmartVolumeMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkSphericalDirectionEncoder", vtkInstantiatorvtkSphericalDirectionEncoderNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeOutlineSource", vtkInstantiatorvtkVolumeOutlineSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumePicker", vtkInstantiatorvtkVolumePickerNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeProMapper", vtkInstantiatorvtkVolumeProMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeRayCastCompositeFunction", vtkInstantiatorvtkVolumeRayCastCompositeFunctionNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeRayCastIsosurfaceFunction", vtkInstantiatorvtkVolumeRayCastIsosurfaceFunctionNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeRayCastMIPFunction", vtkInstantiatorvtkVolumeRayCastMIPFunctionNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeRayCastMapper", vtkInstantiatorvtkVolumeRayCastMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeRenderingFactory", vtkInstantiatorvtkVolumeRenderingFactoryNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeTextureMapper2D", vtkInstantiatorvtkVolumeTextureMapper2DNew);
  vtkInstantiator::RegisterInstantiator("vtkVolumeTextureMapper3D", vtkInstantiatorvtkVolumeTextureMapper3DNew);
  vtkInstantiator::RegisterInstantiator("vtkUnstructuredGridBunykRayCastFunction", vtkInstantiatorvtkUnstructuredGridBunykRayCastFunctionNew);
  vtkInstantiator::RegisterInstantiator("vtkUnstructuredGridHomogeneousRayIntegrator", vtkInstantiatorvtkUnstructuredGridHomogeneousRayIntegratorNew);
  vtkInstantiator::RegisterInstantiator("vtkUnstructuredGridLinearRayIntegrator", vtkInstantiatorvtkUnstructuredGridLinearRayIntegratorNew);
  vtkInstantiator::RegisterInstantiator("vtkUnstructuredGridPartialPreIntegration", vtkInstantiatorvtkUnstructuredGridPartialPreIntegrationNew);
  vtkInstantiator::RegisterInstantiator("vtkUnstructuredGridPreIntegration", vtkInstantiatorvtkUnstructuredGridPreIntegrationNew);
  vtkInstantiator::RegisterInstantiator("vtkUnstructuredGridVolumeRayCastMapper", vtkInstantiatorvtkUnstructuredGridVolumeRayCastMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkUnstructuredGridVolumeZSweepMapper", vtkInstantiatorvtkUnstructuredGridVolumeZSweepMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkOpenGLGPUVolumeRayCastMapper", vtkInstantiatorvtkOpenGLGPUVolumeRayCastMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkOpenGLHAVSVolumeMapper", vtkInstantiatorvtkOpenGLHAVSVolumeMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkOpenGLProjectedAAHexahedraMapper", vtkInstantiatorvtkOpenGLProjectedAAHexahedraMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkOpenGLProjectedTetrahedraMapper", vtkInstantiatorvtkOpenGLProjectedTetrahedraMapperNew);
  vtkInstantiator::RegisterInstantiator("vtkOpenGLRayCastImageDisplayHelper", vtkInstantiatorvtkOpenGLRayCastImageDisplayHelperNew);
  vtkInstantiator::RegisterInstantiator("vtkOpenGLVolumeTextureMapper2D", vtkInstantiatorvtkOpenGLVolumeTextureMapper2DNew);
  vtkInstantiator::RegisterInstantiator("vtkOpenGLVolumeTextureMapper3D", vtkInstantiatorvtkOpenGLVolumeTextureMapper3DNew);

  
}
          
void vtkVolumeRenderingInstantiator::ClassFinalize()
{ 

  vtkInstantiator::UnRegisterInstantiator("vtkEncodedGradientShader", vtkInstantiatorvtkEncodedGradientShaderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFiniteDifferenceGradientEstimator", vtkInstantiatorvtkFiniteDifferenceGradientEstimatorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFixedPointRayCastImage", vtkInstantiatorvtkFixedPointRayCastImageNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFixedPointVolumeRayCastCompositeGOHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeGOHelperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFixedPointVolumeRayCastCompositeGOShadeHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeGOShadeHelperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFixedPointVolumeRayCastCompositeHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeHelperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFixedPointVolumeRayCastCompositeShadeHelper", vtkInstantiatorvtkFixedPointVolumeRayCastCompositeShadeHelperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFixedPointVolumeRayCastMIPHelper", vtkInstantiatorvtkFixedPointVolumeRayCastMIPHelperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFixedPointVolumeRayCastMapper", vtkInstantiatorvtkFixedPointVolumeRayCastMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeRayCastSpaceLeapingImageFilter", vtkInstantiatorvtkVolumeRayCastSpaceLeapingImageFilterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGPUVolumeRayCastMapper", vtkInstantiatorvtkGPUVolumeRayCastMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkHAVSVolumeMapper", vtkInstantiatorvtkHAVSVolumeMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkProjectedAAHexahedraMapper", vtkInstantiatorvtkProjectedAAHexahedraMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkProjectedTetrahedraMapper", vtkInstantiatorvtkProjectedTetrahedraMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRecursiveSphereDirectionEncoder", vtkInstantiatorvtkRecursiveSphereDirectionEncoderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkSmartVolumeMapper", vtkInstantiatorvtkSmartVolumeMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkSphericalDirectionEncoder", vtkInstantiatorvtkSphericalDirectionEncoderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeOutlineSource", vtkInstantiatorvtkVolumeOutlineSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumePicker", vtkInstantiatorvtkVolumePickerNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeProMapper", vtkInstantiatorvtkVolumeProMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeRayCastCompositeFunction", vtkInstantiatorvtkVolumeRayCastCompositeFunctionNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeRayCastIsosurfaceFunction", vtkInstantiatorvtkVolumeRayCastIsosurfaceFunctionNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeRayCastMIPFunction", vtkInstantiatorvtkVolumeRayCastMIPFunctionNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeRayCastMapper", vtkInstantiatorvtkVolumeRayCastMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeRenderingFactory", vtkInstantiatorvtkVolumeRenderingFactoryNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeTextureMapper2D", vtkInstantiatorvtkVolumeTextureMapper2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVolumeTextureMapper3D", vtkInstantiatorvtkVolumeTextureMapper3DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkUnstructuredGridBunykRayCastFunction", vtkInstantiatorvtkUnstructuredGridBunykRayCastFunctionNew);
  vtkInstantiator::UnRegisterInstantiator("vtkUnstructuredGridHomogeneousRayIntegrator", vtkInstantiatorvtkUnstructuredGridHomogeneousRayIntegratorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkUnstructuredGridLinearRayIntegrator", vtkInstantiatorvtkUnstructuredGridLinearRayIntegratorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkUnstructuredGridPartialPreIntegration", vtkInstantiatorvtkUnstructuredGridPartialPreIntegrationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkUnstructuredGridPreIntegration", vtkInstantiatorvtkUnstructuredGridPreIntegrationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkUnstructuredGridVolumeRayCastMapper", vtkInstantiatorvtkUnstructuredGridVolumeRayCastMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkUnstructuredGridVolumeZSweepMapper", vtkInstantiatorvtkUnstructuredGridVolumeZSweepMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkOpenGLGPUVolumeRayCastMapper", vtkInstantiatorvtkOpenGLGPUVolumeRayCastMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkOpenGLHAVSVolumeMapper", vtkInstantiatorvtkOpenGLHAVSVolumeMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkOpenGLProjectedAAHexahedraMapper", vtkInstantiatorvtkOpenGLProjectedAAHexahedraMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkOpenGLProjectedTetrahedraMapper", vtkInstantiatorvtkOpenGLProjectedTetrahedraMapperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkOpenGLRayCastImageDisplayHelper", vtkInstantiatorvtkOpenGLRayCastImageDisplayHelperNew);
  vtkInstantiator::UnRegisterInstantiator("vtkOpenGLVolumeTextureMapper2D", vtkInstantiatorvtkOpenGLVolumeTextureMapper2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkOpenGLVolumeTextureMapper3D", vtkInstantiatorvtkOpenGLVolumeTextureMapper3DNew);

  
}

vtkVolumeRenderingInstantiator::vtkVolumeRenderingInstantiator()
{
  if(++vtkVolumeRenderingInstantiator::Count == 1)
    { 
    vtkVolumeRenderingInstantiator::ClassInitialize(); 
    }
}

vtkVolumeRenderingInstantiator::~vtkVolumeRenderingInstantiator()
{
  if(--vtkVolumeRenderingInstantiator::Count == 0)
    { 
    vtkVolumeRenderingInstantiator::ClassFinalize(); 
    }
}

// Number of translation units that include this class's header.
// Purposely not initialized.  Default is static initialization to 0.
unsigned int vtkVolumeRenderingInstantiator::Count;

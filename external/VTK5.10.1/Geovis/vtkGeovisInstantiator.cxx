#include "vtkGeovisInstantiator.h"
  
extern vtkObject* vtkInstantiatorvtkCompassRepresentationNew();
extern vtkObject* vtkInstantiatorvtkCompassWidgetNew();
extern vtkObject* vtkInstantiatorvtkGeoAdaptiveArcsNew();
extern vtkObject* vtkInstantiatorvtkGeoAlignedImageRepresentationNew();
extern vtkObject* vtkInstantiatorvtkGeoAlignedImageSourceNew();
extern vtkObject* vtkInstantiatorvtkGeoArcsNew();
extern vtkObject* vtkInstantiatorvtkGeoAssignCoordinatesNew();
extern vtkObject* vtkInstantiatorvtkGeoCameraNew();
extern vtkObject* vtkInstantiatorvtkGeoFileImageSourceNew();
extern vtkObject* vtkInstantiatorvtkGeoFileTerrainSourceNew();
extern vtkObject* vtkInstantiatorvtkGeoGlobeSourceNew();
extern vtkObject* vtkInstantiatorvtkGeoGraticuleNew();
extern vtkObject* vtkInstantiatorvtkGeoImageNodeNew();
extern vtkObject* vtkInstantiatorvtkGeoInteractorStyleNew();
extern vtkObject* vtkInstantiatorvtkGeoProjectionSourceNew();
extern vtkObject* vtkInstantiatorvtkGeoProjectionNew();
extern vtkObject* vtkInstantiatorvtkGeoRandomGraphSourceNew();
extern vtkObject* vtkInstantiatorvtkGeoSampleArcsNew();
extern vtkObject* vtkInstantiatorvtkGeoSphereTransformNew();
extern vtkObject* vtkInstantiatorvtkGeoTerrainNew();
extern vtkObject* vtkInstantiatorvtkGeoTerrain2DNew();
extern vtkObject* vtkInstantiatorvtkGeoTerrainNodeNew();
extern vtkObject* vtkInstantiatorvtkGeoTransformNew();
extern vtkObject* vtkInstantiatorvtkGeoTreeNodeNew();
extern vtkObject* vtkInstantiatorvtkGeoTreeNodeCacheNew();
extern vtkObject* vtkInstantiatorvtkGeoViewNew();
extern vtkObject* vtkInstantiatorvtkGeoView2DNew();
extern vtkObject* vtkInstantiatorvtkGlobeSourceNew();


  
void vtkGeovisInstantiator::ClassInitialize()
{
  
  vtkInstantiator::RegisterInstantiator("vtkCompassRepresentation", vtkInstantiatorvtkCompassRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkCompassWidget", vtkInstantiatorvtkCompassWidgetNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoAdaptiveArcs", vtkInstantiatorvtkGeoAdaptiveArcsNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoAlignedImageRepresentation", vtkInstantiatorvtkGeoAlignedImageRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoAlignedImageSource", vtkInstantiatorvtkGeoAlignedImageSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoArcs", vtkInstantiatorvtkGeoArcsNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoAssignCoordinates", vtkInstantiatorvtkGeoAssignCoordinatesNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoCamera", vtkInstantiatorvtkGeoCameraNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoFileImageSource", vtkInstantiatorvtkGeoFileImageSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoFileTerrainSource", vtkInstantiatorvtkGeoFileTerrainSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoGlobeSource", vtkInstantiatorvtkGeoGlobeSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoGraticule", vtkInstantiatorvtkGeoGraticuleNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoImageNode", vtkInstantiatorvtkGeoImageNodeNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoInteractorStyle", vtkInstantiatorvtkGeoInteractorStyleNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoProjectionSource", vtkInstantiatorvtkGeoProjectionSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoProjection", vtkInstantiatorvtkGeoProjectionNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoRandomGraphSource", vtkInstantiatorvtkGeoRandomGraphSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoSampleArcs", vtkInstantiatorvtkGeoSampleArcsNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoSphereTransform", vtkInstantiatorvtkGeoSphereTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoTerrain", vtkInstantiatorvtkGeoTerrainNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoTerrain2D", vtkInstantiatorvtkGeoTerrain2DNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoTerrainNode", vtkInstantiatorvtkGeoTerrainNodeNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoTransform", vtkInstantiatorvtkGeoTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoTreeNode", vtkInstantiatorvtkGeoTreeNodeNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoTreeNodeCache", vtkInstantiatorvtkGeoTreeNodeCacheNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoView", vtkInstantiatorvtkGeoViewNew);
  vtkInstantiator::RegisterInstantiator("vtkGeoView2D", vtkInstantiatorvtkGeoView2DNew);
  vtkInstantiator::RegisterInstantiator("vtkGlobeSource", vtkInstantiatorvtkGlobeSourceNew);

  
}
          
void vtkGeovisInstantiator::ClassFinalize()
{ 

  vtkInstantiator::UnRegisterInstantiator("vtkCompassRepresentation", vtkInstantiatorvtkCompassRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkCompassWidget", vtkInstantiatorvtkCompassWidgetNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoAdaptiveArcs", vtkInstantiatorvtkGeoAdaptiveArcsNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoAlignedImageRepresentation", vtkInstantiatorvtkGeoAlignedImageRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoAlignedImageSource", vtkInstantiatorvtkGeoAlignedImageSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoArcs", vtkInstantiatorvtkGeoArcsNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoAssignCoordinates", vtkInstantiatorvtkGeoAssignCoordinatesNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoCamera", vtkInstantiatorvtkGeoCameraNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoFileImageSource", vtkInstantiatorvtkGeoFileImageSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoFileTerrainSource", vtkInstantiatorvtkGeoFileTerrainSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoGlobeSource", vtkInstantiatorvtkGeoGlobeSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoGraticule", vtkInstantiatorvtkGeoGraticuleNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoImageNode", vtkInstantiatorvtkGeoImageNodeNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoInteractorStyle", vtkInstantiatorvtkGeoInteractorStyleNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoProjectionSource", vtkInstantiatorvtkGeoProjectionSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoProjection", vtkInstantiatorvtkGeoProjectionNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoRandomGraphSource", vtkInstantiatorvtkGeoRandomGraphSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoSampleArcs", vtkInstantiatorvtkGeoSampleArcsNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoSphereTransform", vtkInstantiatorvtkGeoSphereTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoTerrain", vtkInstantiatorvtkGeoTerrainNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoTerrain2D", vtkInstantiatorvtkGeoTerrain2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoTerrainNode", vtkInstantiatorvtkGeoTerrainNodeNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoTransform", vtkInstantiatorvtkGeoTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoTreeNode", vtkInstantiatorvtkGeoTreeNodeNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoTreeNodeCache", vtkInstantiatorvtkGeoTreeNodeCacheNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoView", vtkInstantiatorvtkGeoViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGeoView2D", vtkInstantiatorvtkGeoView2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGlobeSource", vtkInstantiatorvtkGlobeSourceNew);

  
}

vtkGeovisInstantiator::vtkGeovisInstantiator()
{
  if(++vtkGeovisInstantiator::Count == 1)
    { 
    vtkGeovisInstantiator::ClassInitialize(); 
    }
}

vtkGeovisInstantiator::~vtkGeovisInstantiator()
{
  if(--vtkGeovisInstantiator::Count == 0)
    { 
    vtkGeovisInstantiator::ClassFinalize(); 
    }
}

// Number of translation units that include this class's header.
// Purposely not initialized.  Default is static initialization to 0.
unsigned int vtkGeovisInstantiator::Count;

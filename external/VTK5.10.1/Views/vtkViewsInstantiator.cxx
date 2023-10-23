#include "vtkViewsInstantiator.h"
  
extern vtkObject* vtkInstantiatorvtkConvertSelectionDomainNew();
extern vtkObject* vtkInstantiatorvtkDataRepresentationNew();
extern vtkObject* vtkInstantiatorvtkEmptyRepresentationNew();
extern vtkObject* vtkInstantiatorvtkGraphLayoutViewNew();
extern vtkObject* vtkInstantiatorvtkHierarchicalGraphPipelineNew();
extern vtkObject* vtkInstantiatorvtkHierarchicalGraphViewNew();
extern vtkObject* vtkInstantiatorvtkIcicleViewNew();
extern vtkObject* vtkInstantiatorvtkInteractorStyleAreaSelectHoverNew();
extern vtkObject* vtkInstantiatorvtkInteractorStyleTreeMapHoverNew();
extern vtkObject* vtkInstantiatorvtkRenderedSurfaceRepresentationNew();
extern vtkObject* vtkInstantiatorvtkRenderedGraphRepresentationNew();
extern vtkObject* vtkInstantiatorvtkRenderedRepresentationNew();
extern vtkObject* vtkInstantiatorvtkRenderedTreeAreaRepresentationNew();
extern vtkObject* vtkInstantiatorvtkRenderedHierarchyRepresentationNew();
extern vtkObject* vtkInstantiatorvtkRenderViewNew();
extern vtkObject* vtkInstantiatorvtkRenderViewBaseNew();
extern vtkObject* vtkInstantiatorvtkTreeAreaViewNew();
extern vtkObject* vtkInstantiatorvtkTreeMapViewNew();
extern vtkObject* vtkInstantiatorvtkTreeRingViewNew();
extern vtkObject* vtkInstantiatorvtkViewNew();
extern vtkObject* vtkInstantiatorvtkViewUpdaterNew();
extern vtkObject* vtkInstantiatorvtkParallelCoordinatesHistogramRepresentationNew();
extern vtkObject* vtkInstantiatorvtkParallelCoordinatesRepresentationNew();
extern vtkObject* vtkInstantiatorvtkParallelCoordinatesViewNew();


  
void vtkViewsInstantiator::ClassInitialize()
{
  
  vtkInstantiator::RegisterInstantiator("vtkConvertSelectionDomain", vtkInstantiatorvtkConvertSelectionDomainNew);
  vtkInstantiator::RegisterInstantiator("vtkDataRepresentation", vtkInstantiatorvtkDataRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkEmptyRepresentation", vtkInstantiatorvtkEmptyRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkGraphLayoutView", vtkInstantiatorvtkGraphLayoutViewNew);
  vtkInstantiator::RegisterInstantiator("vtkHierarchicalGraphPipeline", vtkInstantiatorvtkHierarchicalGraphPipelineNew);
  vtkInstantiator::RegisterInstantiator("vtkHierarchicalGraphView", vtkInstantiatorvtkHierarchicalGraphViewNew);
  vtkInstantiator::RegisterInstantiator("vtkIcicleView", vtkInstantiatorvtkIcicleViewNew);
  vtkInstantiator::RegisterInstantiator("vtkInteractorStyleAreaSelectHover", vtkInstantiatorvtkInteractorStyleAreaSelectHoverNew);
  vtkInstantiator::RegisterInstantiator("vtkInteractorStyleTreeMapHover", vtkInstantiatorvtkInteractorStyleTreeMapHoverNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderedSurfaceRepresentation", vtkInstantiatorvtkRenderedSurfaceRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderedGraphRepresentation", vtkInstantiatorvtkRenderedGraphRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderedRepresentation", vtkInstantiatorvtkRenderedRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderedTreeAreaRepresentation", vtkInstantiatorvtkRenderedTreeAreaRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderedHierarchyRepresentation", vtkInstantiatorvtkRenderedHierarchyRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderView", vtkInstantiatorvtkRenderViewNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderViewBase", vtkInstantiatorvtkRenderViewBaseNew);
  vtkInstantiator::RegisterInstantiator("vtkTreeAreaView", vtkInstantiatorvtkTreeAreaViewNew);
  vtkInstantiator::RegisterInstantiator("vtkTreeMapView", vtkInstantiatorvtkTreeMapViewNew);
  vtkInstantiator::RegisterInstantiator("vtkTreeRingView", vtkInstantiatorvtkTreeRingViewNew);
  vtkInstantiator::RegisterInstantiator("vtkView", vtkInstantiatorvtkViewNew);
  vtkInstantiator::RegisterInstantiator("vtkViewUpdater", vtkInstantiatorvtkViewUpdaterNew);
  vtkInstantiator::RegisterInstantiator("vtkParallelCoordinatesHistogramRepresentation", vtkInstantiatorvtkParallelCoordinatesHistogramRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkParallelCoordinatesRepresentation", vtkInstantiatorvtkParallelCoordinatesRepresentationNew);
  vtkInstantiator::RegisterInstantiator("vtkParallelCoordinatesView", vtkInstantiatorvtkParallelCoordinatesViewNew);

  
}
          
void vtkViewsInstantiator::ClassFinalize()
{ 

  vtkInstantiator::UnRegisterInstantiator("vtkConvertSelectionDomain", vtkInstantiatorvtkConvertSelectionDomainNew);
  vtkInstantiator::UnRegisterInstantiator("vtkDataRepresentation", vtkInstantiatorvtkDataRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkEmptyRepresentation", vtkInstantiatorvtkEmptyRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGraphLayoutView", vtkInstantiatorvtkGraphLayoutViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkHierarchicalGraphPipeline", vtkInstantiatorvtkHierarchicalGraphPipelineNew);
  vtkInstantiator::UnRegisterInstantiator("vtkHierarchicalGraphView", vtkInstantiatorvtkHierarchicalGraphViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkIcicleView", vtkInstantiatorvtkIcicleViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkInteractorStyleAreaSelectHover", vtkInstantiatorvtkInteractorStyleAreaSelectHoverNew);
  vtkInstantiator::UnRegisterInstantiator("vtkInteractorStyleTreeMapHover", vtkInstantiatorvtkInteractorStyleTreeMapHoverNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderedSurfaceRepresentation", vtkInstantiatorvtkRenderedSurfaceRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderedGraphRepresentation", vtkInstantiatorvtkRenderedGraphRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderedRepresentation", vtkInstantiatorvtkRenderedRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderedTreeAreaRepresentation", vtkInstantiatorvtkRenderedTreeAreaRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderedHierarchyRepresentation", vtkInstantiatorvtkRenderedHierarchyRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderView", vtkInstantiatorvtkRenderViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderViewBase", vtkInstantiatorvtkRenderViewBaseNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTreeAreaView", vtkInstantiatorvtkTreeAreaViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTreeMapView", vtkInstantiatorvtkTreeMapViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTreeRingView", vtkInstantiatorvtkTreeRingViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkView", vtkInstantiatorvtkViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkViewUpdater", vtkInstantiatorvtkViewUpdaterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkParallelCoordinatesHistogramRepresentation", vtkInstantiatorvtkParallelCoordinatesHistogramRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkParallelCoordinatesRepresentation", vtkInstantiatorvtkParallelCoordinatesRepresentationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkParallelCoordinatesView", vtkInstantiatorvtkParallelCoordinatesViewNew);

  
}

vtkViewsInstantiator::vtkViewsInstantiator()
{
  if(++vtkViewsInstantiator::Count == 1)
    { 
    vtkViewsInstantiator::ClassInitialize(); 
    }
}

vtkViewsInstantiator::~vtkViewsInstantiator()
{
  if(--vtkViewsInstantiator::Count == 0)
    { 
    vtkViewsInstantiator::ClassFinalize(); 
    }
}

// Number of translation units that include this class's header.
// Purposely not initialized.  Default is static initialization to 0.
unsigned int vtkViewsInstantiator::Count;

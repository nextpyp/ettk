#include "vtkHybridInstantiator.h"
  
extern vtkObject* vtkInstantiatorvtk3DSImporterNew();
extern vtkObject* vtkInstantiatorvtkArcPlotterNew();
extern vtkObject* vtkInstantiatorvtkAnnotatedCubeActorNew();
extern vtkObject* vtkInstantiatorvtkAxisActorNew();
extern vtkObject* vtkInstantiatorvtkAxesActorNew();
extern vtkObject* vtkInstantiatorvtkAxisFollowerNew();
extern vtkObject* vtkInstantiatorvtkBarChartActorNew();
extern vtkObject* vtkInstantiatorvtkBSplineTransformNew();
extern vtkObject* vtkInstantiatorvtkCaptionActor2DNew();
extern vtkObject* vtkInstantiatorvtkCornerAnnotationNew();
extern vtkObject* vtkInstantiatorvtkCubeAxesActorNew();
extern vtkObject* vtkInstantiatorvtkCubeAxesActor2DNew();
extern vtkObject* vtkInstantiatorvtkDepthSortPolyDataNew();
extern vtkObject* vtkInstantiatorvtkEarthSourceNew();
extern vtkObject* vtkInstantiatorvtkFacetReaderNew();
extern vtkObject* vtkInstantiatorvtkGreedyTerrainDecimationNew();
extern vtkObject* vtkInstantiatorvtkGridTransformNew();
extern vtkObject* vtkInstantiatorvtkImageDataLIC2DNew();
extern vtkObject* vtkInstantiatorvtkImageDataLIC2DExtentTranslatorNew();
extern vtkObject* vtkInstantiatorvtkImageToPolyDataFilterNew();
extern vtkObject* vtkInstantiatorvtkImplicitModellerNew();
extern vtkObject* vtkInstantiatorvtkIterativeClosestPointTransformNew();
extern vtkObject* vtkInstantiatorvtkLandmarkTransformNew();
extern vtkObject* vtkInstantiatorvtkLegendBoxActorNew();
extern vtkObject* vtkInstantiatorvtkLegendScaleActorNew();
extern vtkObject* vtkInstantiatorvtkMNIObjectReaderNew();
extern vtkObject* vtkInstantiatorvtkMNIObjectWriterNew();
extern vtkObject* vtkInstantiatorvtkMNITagPointReaderNew();
extern vtkObject* vtkInstantiatorvtkMNITagPointWriterNew();
extern vtkObject* vtkInstantiatorvtkMNITransformReaderNew();
extern vtkObject* vtkInstantiatorvtkMNITransformWriterNew();
extern vtkObject* vtkInstantiatorvtkPCAAnalysisFilterNew();
extern vtkObject* vtkInstantiatorvtkPieChartActorNew();
extern vtkObject* vtkInstantiatorvtkPolarAxesActorNew();
extern vtkObject* vtkInstantiatorvtkPolyDataSilhouetteNew();
extern vtkObject* vtkInstantiatorvtkPolyDataToImageStencilNew();
extern vtkObject* vtkInstantiatorvtkProcrustesAlignmentFilterNew();
extern vtkObject* vtkInstantiatorvtkProjectedTerrainPathNew();
extern vtkObject* vtkInstantiatorvtkRIBExporterNew();
extern vtkObject* vtkInstantiatorvtkRIBLightNew();
extern vtkObject* vtkInstantiatorvtkRIBPropertyNew();
extern vtkObject* vtkInstantiatorvtkRenderLargeImageNew();
extern vtkObject* vtkInstantiatorvtkSpiderPlotActorNew();
extern vtkObject* vtkInstantiatorvtkTemporalDataSetCacheNew();
extern vtkObject* vtkInstantiatorvtkTemporalInterpolatorNew();
extern vtkObject* vtkInstantiatorvtkTemporalShiftScaleNew();
extern vtkObject* vtkInstantiatorvtkTemporalSnapToTimeStepNew();
extern vtkObject* vtkInstantiatorvtkThinPlateSplineTransformNew();
extern vtkObject* vtkInstantiatorvtkTransformToGridNew();
extern vtkObject* vtkInstantiatorvtkVRMLImporterNew();
extern vtkObject* vtkInstantiatorvtkVectorTextNew();
extern vtkObject* vtkInstantiatorvtkVideoSourceNew();
extern vtkObject* vtkInstantiatorvtkWeightedTransformFilterNew();
extern vtkObject* vtkInstantiatorvtkXYPlotActorNew();
extern vtkObject* vtkInstantiatorvtkX3DExporterNew();
extern vtkObject* vtkInstantiatorvtkExodusIICacheNew();
extern vtkObject* vtkInstantiatorvtkExodusIIReaderNew();
extern vtkObject* vtkInstantiatorvtkDSPFilterDefinitionNew();
extern vtkObject* vtkInstantiatorvtkExodusModelNew();
extern vtkObject* vtkInstantiatorvtkDSPFilterGroupNew();
extern vtkObject* vtkInstantiatorvtkExodusReaderNew();


  
void vtkHybridInstantiator::ClassInitialize()
{
  
  vtkInstantiator::RegisterInstantiator("vtk3DSImporter", vtkInstantiatorvtk3DSImporterNew);
  vtkInstantiator::RegisterInstantiator("vtkArcPlotter", vtkInstantiatorvtkArcPlotterNew);
  vtkInstantiator::RegisterInstantiator("vtkAnnotatedCubeActor", vtkInstantiatorvtkAnnotatedCubeActorNew);
  vtkInstantiator::RegisterInstantiator("vtkAxisActor", vtkInstantiatorvtkAxisActorNew);
  vtkInstantiator::RegisterInstantiator("vtkAxesActor", vtkInstantiatorvtkAxesActorNew);
  vtkInstantiator::RegisterInstantiator("vtkAxisFollower", vtkInstantiatorvtkAxisFollowerNew);
  vtkInstantiator::RegisterInstantiator("vtkBarChartActor", vtkInstantiatorvtkBarChartActorNew);
  vtkInstantiator::RegisterInstantiator("vtkBSplineTransform", vtkInstantiatorvtkBSplineTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkCaptionActor2D", vtkInstantiatorvtkCaptionActor2DNew);
  vtkInstantiator::RegisterInstantiator("vtkCornerAnnotation", vtkInstantiatorvtkCornerAnnotationNew);
  vtkInstantiator::RegisterInstantiator("vtkCubeAxesActor", vtkInstantiatorvtkCubeAxesActorNew);
  vtkInstantiator::RegisterInstantiator("vtkCubeAxesActor2D", vtkInstantiatorvtkCubeAxesActor2DNew);
  vtkInstantiator::RegisterInstantiator("vtkDepthSortPolyData", vtkInstantiatorvtkDepthSortPolyDataNew);
  vtkInstantiator::RegisterInstantiator("vtkEarthSource", vtkInstantiatorvtkEarthSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkFacetReader", vtkInstantiatorvtkFacetReaderNew);
  vtkInstantiator::RegisterInstantiator("vtkGreedyTerrainDecimation", vtkInstantiatorvtkGreedyTerrainDecimationNew);
  vtkInstantiator::RegisterInstantiator("vtkGridTransform", vtkInstantiatorvtkGridTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkImageDataLIC2D", vtkInstantiatorvtkImageDataLIC2DNew);
  vtkInstantiator::RegisterInstantiator("vtkImageDataLIC2DExtentTranslator", vtkInstantiatorvtkImageDataLIC2DExtentTranslatorNew);
  vtkInstantiator::RegisterInstantiator("vtkImageToPolyDataFilter", vtkInstantiatorvtkImageToPolyDataFilterNew);
  vtkInstantiator::RegisterInstantiator("vtkImplicitModeller", vtkInstantiatorvtkImplicitModellerNew);
  vtkInstantiator::RegisterInstantiator("vtkIterativeClosestPointTransform", vtkInstantiatorvtkIterativeClosestPointTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkLandmarkTransform", vtkInstantiatorvtkLandmarkTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkLegendBoxActor", vtkInstantiatorvtkLegendBoxActorNew);
  vtkInstantiator::RegisterInstantiator("vtkLegendScaleActor", vtkInstantiatorvtkLegendScaleActorNew);
  vtkInstantiator::RegisterInstantiator("vtkMNIObjectReader", vtkInstantiatorvtkMNIObjectReaderNew);
  vtkInstantiator::RegisterInstantiator("vtkMNIObjectWriter", vtkInstantiatorvtkMNIObjectWriterNew);
  vtkInstantiator::RegisterInstantiator("vtkMNITagPointReader", vtkInstantiatorvtkMNITagPointReaderNew);
  vtkInstantiator::RegisterInstantiator("vtkMNITagPointWriter", vtkInstantiatorvtkMNITagPointWriterNew);
  vtkInstantiator::RegisterInstantiator("vtkMNITransformReader", vtkInstantiatorvtkMNITransformReaderNew);
  vtkInstantiator::RegisterInstantiator("vtkMNITransformWriter", vtkInstantiatorvtkMNITransformWriterNew);
  vtkInstantiator::RegisterInstantiator("vtkPCAAnalysisFilter", vtkInstantiatorvtkPCAAnalysisFilterNew);
  vtkInstantiator::RegisterInstantiator("vtkPieChartActor", vtkInstantiatorvtkPieChartActorNew);
  vtkInstantiator::RegisterInstantiator("vtkPolarAxesActor", vtkInstantiatorvtkPolarAxesActorNew);
  vtkInstantiator::RegisterInstantiator("vtkPolyDataSilhouette", vtkInstantiatorvtkPolyDataSilhouetteNew);
  vtkInstantiator::RegisterInstantiator("vtkPolyDataToImageStencil", vtkInstantiatorvtkPolyDataToImageStencilNew);
  vtkInstantiator::RegisterInstantiator("vtkProcrustesAlignmentFilter", vtkInstantiatorvtkProcrustesAlignmentFilterNew);
  vtkInstantiator::RegisterInstantiator("vtkProjectedTerrainPath", vtkInstantiatorvtkProjectedTerrainPathNew);
  vtkInstantiator::RegisterInstantiator("vtkRIBExporter", vtkInstantiatorvtkRIBExporterNew);
  vtkInstantiator::RegisterInstantiator("vtkRIBLight", vtkInstantiatorvtkRIBLightNew);
  vtkInstantiator::RegisterInstantiator("vtkRIBProperty", vtkInstantiatorvtkRIBPropertyNew);
  vtkInstantiator::RegisterInstantiator("vtkRenderLargeImage", vtkInstantiatorvtkRenderLargeImageNew);
  vtkInstantiator::RegisterInstantiator("vtkSpiderPlotActor", vtkInstantiatorvtkSpiderPlotActorNew);
  vtkInstantiator::RegisterInstantiator("vtkTemporalDataSetCache", vtkInstantiatorvtkTemporalDataSetCacheNew);
  vtkInstantiator::RegisterInstantiator("vtkTemporalInterpolator", vtkInstantiatorvtkTemporalInterpolatorNew);
  vtkInstantiator::RegisterInstantiator("vtkTemporalShiftScale", vtkInstantiatorvtkTemporalShiftScaleNew);
  vtkInstantiator::RegisterInstantiator("vtkTemporalSnapToTimeStep", vtkInstantiatorvtkTemporalSnapToTimeStepNew);
  vtkInstantiator::RegisterInstantiator("vtkThinPlateSplineTransform", vtkInstantiatorvtkThinPlateSplineTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkTransformToGrid", vtkInstantiatorvtkTransformToGridNew);
  vtkInstantiator::RegisterInstantiator("vtkVRMLImporter", vtkInstantiatorvtkVRMLImporterNew);
  vtkInstantiator::RegisterInstantiator("vtkVectorText", vtkInstantiatorvtkVectorTextNew);
  vtkInstantiator::RegisterInstantiator("vtkVideoSource", vtkInstantiatorvtkVideoSourceNew);
  vtkInstantiator::RegisterInstantiator("vtkWeightedTransformFilter", vtkInstantiatorvtkWeightedTransformFilterNew);
  vtkInstantiator::RegisterInstantiator("vtkXYPlotActor", vtkInstantiatorvtkXYPlotActorNew);
  vtkInstantiator::RegisterInstantiator("vtkX3DExporter", vtkInstantiatorvtkX3DExporterNew);
  vtkInstantiator::RegisterInstantiator("vtkExodusIICache", vtkInstantiatorvtkExodusIICacheNew);
  vtkInstantiator::RegisterInstantiator("vtkExodusIIReader", vtkInstantiatorvtkExodusIIReaderNew);
  vtkInstantiator::RegisterInstantiator("vtkDSPFilterDefinition", vtkInstantiatorvtkDSPFilterDefinitionNew);
  vtkInstantiator::RegisterInstantiator("vtkExodusModel", vtkInstantiatorvtkExodusModelNew);
  vtkInstantiator::RegisterInstantiator("vtkDSPFilterGroup", vtkInstantiatorvtkDSPFilterGroupNew);
  vtkInstantiator::RegisterInstantiator("vtkExodusReader", vtkInstantiatorvtkExodusReaderNew);

  
}
          
void vtkHybridInstantiator::ClassFinalize()
{ 

  vtkInstantiator::UnRegisterInstantiator("vtk3DSImporter", vtkInstantiatorvtk3DSImporterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkArcPlotter", vtkInstantiatorvtkArcPlotterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkAnnotatedCubeActor", vtkInstantiatorvtkAnnotatedCubeActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkAxisActor", vtkInstantiatorvtkAxisActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkAxesActor", vtkInstantiatorvtkAxesActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkAxisFollower", vtkInstantiatorvtkAxisFollowerNew);
  vtkInstantiator::UnRegisterInstantiator("vtkBarChartActor", vtkInstantiatorvtkBarChartActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkBSplineTransform", vtkInstantiatorvtkBSplineTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkCaptionActor2D", vtkInstantiatorvtkCaptionActor2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkCornerAnnotation", vtkInstantiatorvtkCornerAnnotationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkCubeAxesActor", vtkInstantiatorvtkCubeAxesActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkCubeAxesActor2D", vtkInstantiatorvtkCubeAxesActor2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkDepthSortPolyData", vtkInstantiatorvtkDepthSortPolyDataNew);
  vtkInstantiator::UnRegisterInstantiator("vtkEarthSource", vtkInstantiatorvtkEarthSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkFacetReader", vtkInstantiatorvtkFacetReaderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGreedyTerrainDecimation", vtkInstantiatorvtkGreedyTerrainDecimationNew);
  vtkInstantiator::UnRegisterInstantiator("vtkGridTransform", vtkInstantiatorvtkGridTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkImageDataLIC2D", vtkInstantiatorvtkImageDataLIC2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkImageDataLIC2DExtentTranslator", vtkInstantiatorvtkImageDataLIC2DExtentTranslatorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkImageToPolyDataFilter", vtkInstantiatorvtkImageToPolyDataFilterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkImplicitModeller", vtkInstantiatorvtkImplicitModellerNew);
  vtkInstantiator::UnRegisterInstantiator("vtkIterativeClosestPointTransform", vtkInstantiatorvtkIterativeClosestPointTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkLandmarkTransform", vtkInstantiatorvtkLandmarkTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkLegendBoxActor", vtkInstantiatorvtkLegendBoxActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkLegendScaleActor", vtkInstantiatorvtkLegendScaleActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkMNIObjectReader", vtkInstantiatorvtkMNIObjectReaderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkMNIObjectWriter", vtkInstantiatorvtkMNIObjectWriterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkMNITagPointReader", vtkInstantiatorvtkMNITagPointReaderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkMNITagPointWriter", vtkInstantiatorvtkMNITagPointWriterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkMNITransformReader", vtkInstantiatorvtkMNITransformReaderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkMNITransformWriter", vtkInstantiatorvtkMNITransformWriterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPCAAnalysisFilter", vtkInstantiatorvtkPCAAnalysisFilterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPieChartActor", vtkInstantiatorvtkPieChartActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPolarAxesActor", vtkInstantiatorvtkPolarAxesActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPolyDataSilhouette", vtkInstantiatorvtkPolyDataSilhouetteNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPolyDataToImageStencil", vtkInstantiatorvtkPolyDataToImageStencilNew);
  vtkInstantiator::UnRegisterInstantiator("vtkProcrustesAlignmentFilter", vtkInstantiatorvtkProcrustesAlignmentFilterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkProjectedTerrainPath", vtkInstantiatorvtkProjectedTerrainPathNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRIBExporter", vtkInstantiatorvtkRIBExporterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRIBLight", vtkInstantiatorvtkRIBLightNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRIBProperty", vtkInstantiatorvtkRIBPropertyNew);
  vtkInstantiator::UnRegisterInstantiator("vtkRenderLargeImage", vtkInstantiatorvtkRenderLargeImageNew);
  vtkInstantiator::UnRegisterInstantiator("vtkSpiderPlotActor", vtkInstantiatorvtkSpiderPlotActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTemporalDataSetCache", vtkInstantiatorvtkTemporalDataSetCacheNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTemporalInterpolator", vtkInstantiatorvtkTemporalInterpolatorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTemporalShiftScale", vtkInstantiatorvtkTemporalShiftScaleNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTemporalSnapToTimeStep", vtkInstantiatorvtkTemporalSnapToTimeStepNew);
  vtkInstantiator::UnRegisterInstantiator("vtkThinPlateSplineTransform", vtkInstantiatorvtkThinPlateSplineTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTransformToGrid", vtkInstantiatorvtkTransformToGridNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVRMLImporter", vtkInstantiatorvtkVRMLImporterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVectorText", vtkInstantiatorvtkVectorTextNew);
  vtkInstantiator::UnRegisterInstantiator("vtkVideoSource", vtkInstantiatorvtkVideoSourceNew);
  vtkInstantiator::UnRegisterInstantiator("vtkWeightedTransformFilter", vtkInstantiatorvtkWeightedTransformFilterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkXYPlotActor", vtkInstantiatorvtkXYPlotActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkX3DExporter", vtkInstantiatorvtkX3DExporterNew);
  vtkInstantiator::UnRegisterInstantiator("vtkExodusIICache", vtkInstantiatorvtkExodusIICacheNew);
  vtkInstantiator::UnRegisterInstantiator("vtkExodusIIReader", vtkInstantiatorvtkExodusIIReaderNew);
  vtkInstantiator::UnRegisterInstantiator("vtkDSPFilterDefinition", vtkInstantiatorvtkDSPFilterDefinitionNew);
  vtkInstantiator::UnRegisterInstantiator("vtkExodusModel", vtkInstantiatorvtkExodusModelNew);
  vtkInstantiator::UnRegisterInstantiator("vtkDSPFilterGroup", vtkInstantiatorvtkDSPFilterGroupNew);
  vtkInstantiator::UnRegisterInstantiator("vtkExodusReader", vtkInstantiatorvtkExodusReaderNew);

  
}

vtkHybridInstantiator::vtkHybridInstantiator()
{
  if(++vtkHybridInstantiator::Count == 1)
    { 
    vtkHybridInstantiator::ClassInitialize(); 
    }
}

vtkHybridInstantiator::~vtkHybridInstantiator()
{
  if(--vtkHybridInstantiator::Count == 0)
    { 
    vtkHybridInstantiator::ClassFinalize(); 
    }
}

// Number of translation units that include this class's header.
// Purposely not initialized.  Default is static initialization to 0.
unsigned int vtkHybridInstantiator::Count;

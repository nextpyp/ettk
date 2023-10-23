#include "vtkChartsInstantiator.h"
  
extern vtkObject* vtkInstantiatorvtkAxisNew();
extern vtkObject* vtkInstantiatorvtkAxisExtendedNew();
extern vtkObject* vtkInstantiatorvtkBlockItemNew();
extern vtkObject* vtkInstantiatorvtkBrushNew();
extern vtkObject* vtkInstantiatorvtkChartLegendNew();
extern vtkObject* vtkInstantiatorvtkChartHistogram2DNew();
extern vtkObject* vtkInstantiatorvtkChartMatrixNew();
extern vtkObject* vtkInstantiatorvtkChartParallelCoordinatesNew();
extern vtkObject* vtkInstantiatorvtkChartXYNew();
extern vtkObject* vtkInstantiatorvtkChartPieNew();
extern vtkObject* vtkInstantiatorvtkColorLegendNew();
extern vtkObject* vtkInstantiatorvtkPlotPieNew();
extern vtkObject* vtkInstantiatorvtkColorSeriesNew();
extern vtkObject* vtkInstantiatorvtkColorTransferFunctionItemNew();
extern vtkObject* vtkInstantiatorvtkColorTransferControlPointsItemNew();
extern vtkObject* vtkInstantiatorvtkCompositeControlPointsItemNew();
extern vtkObject* vtkInstantiatorvtkCompositeTransferFunctionItemNew();
extern vtkObject* vtkInstantiatorvtkContext2DNew();
extern vtkObject* vtkInstantiatorvtkContextActorNew();
extern vtkObject* vtkInstantiatorvtkContextClipNew();
extern vtkObject* vtkInstantiatorvtkContextInteractorStyleNew();
extern vtkObject* vtkInstantiatorvtkContextSceneNew();
extern vtkObject* vtkInstantiatorvtkContextTransformNew();
extern vtkObject* vtkInstantiatorvtkContextViewNew();
extern vtkObject* vtkInstantiatorvtkImageItemNew();
extern vtkObject* vtkInstantiatorvtkLookupTableItemNew();
extern vtkObject* vtkInstantiatorvtkPenNew();
extern vtkObject* vtkInstantiatorvtkPiecewiseControlPointsItemNew();
extern vtkObject* vtkInstantiatorvtkPiecewiseFunctionItemNew();
extern vtkObject* vtkInstantiatorvtkPiecewisePointHandleItemNew();
extern vtkObject* vtkInstantiatorvtkPlotBarNew();
extern vtkObject* vtkInstantiatorvtkPlotGridNew();
extern vtkObject* vtkInstantiatorvtkPlotHistogram2DNew();
extern vtkObject* vtkInstantiatorvtkPlotLineNew();
extern vtkObject* vtkInstantiatorvtkPlotStackedNew();
extern vtkObject* vtkInstantiatorvtkPlotParallelCoordinatesNew();
extern vtkObject* vtkInstantiatorvtkPlotPointsNew();
extern vtkObject* vtkInstantiatorvtkScatterPlotMatrixNew();
extern vtkObject* vtkInstantiatorvtkTooltipItemNew();


  
void vtkChartsInstantiator::ClassInitialize()
{
  
  vtkInstantiator::RegisterInstantiator("vtkAxis", vtkInstantiatorvtkAxisNew);
  vtkInstantiator::RegisterInstantiator("vtkAxisExtended", vtkInstantiatorvtkAxisExtendedNew);
  vtkInstantiator::RegisterInstantiator("vtkBlockItem", vtkInstantiatorvtkBlockItemNew);
  vtkInstantiator::RegisterInstantiator("vtkBrush", vtkInstantiatorvtkBrushNew);
  vtkInstantiator::RegisterInstantiator("vtkChartLegend", vtkInstantiatorvtkChartLegendNew);
  vtkInstantiator::RegisterInstantiator("vtkChartHistogram2D", vtkInstantiatorvtkChartHistogram2DNew);
  vtkInstantiator::RegisterInstantiator("vtkChartMatrix", vtkInstantiatorvtkChartMatrixNew);
  vtkInstantiator::RegisterInstantiator("vtkChartParallelCoordinates", vtkInstantiatorvtkChartParallelCoordinatesNew);
  vtkInstantiator::RegisterInstantiator("vtkChartXY", vtkInstantiatorvtkChartXYNew);
  vtkInstantiator::RegisterInstantiator("vtkChartPie", vtkInstantiatorvtkChartPieNew);
  vtkInstantiator::RegisterInstantiator("vtkColorLegend", vtkInstantiatorvtkColorLegendNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotPie", vtkInstantiatorvtkPlotPieNew);
  vtkInstantiator::RegisterInstantiator("vtkColorSeries", vtkInstantiatorvtkColorSeriesNew);
  vtkInstantiator::RegisterInstantiator("vtkColorTransferFunctionItem", vtkInstantiatorvtkColorTransferFunctionItemNew);
  vtkInstantiator::RegisterInstantiator("vtkColorTransferControlPointsItem", vtkInstantiatorvtkColorTransferControlPointsItemNew);
  vtkInstantiator::RegisterInstantiator("vtkCompositeControlPointsItem", vtkInstantiatorvtkCompositeControlPointsItemNew);
  vtkInstantiator::RegisterInstantiator("vtkCompositeTransferFunctionItem", vtkInstantiatorvtkCompositeTransferFunctionItemNew);
  vtkInstantiator::RegisterInstantiator("vtkContext2D", vtkInstantiatorvtkContext2DNew);
  vtkInstantiator::RegisterInstantiator("vtkContextActor", vtkInstantiatorvtkContextActorNew);
  vtkInstantiator::RegisterInstantiator("vtkContextClip", vtkInstantiatorvtkContextClipNew);
  vtkInstantiator::RegisterInstantiator("vtkContextInteractorStyle", vtkInstantiatorvtkContextInteractorStyleNew);
  vtkInstantiator::RegisterInstantiator("vtkContextScene", vtkInstantiatorvtkContextSceneNew);
  vtkInstantiator::RegisterInstantiator("vtkContextTransform", vtkInstantiatorvtkContextTransformNew);
  vtkInstantiator::RegisterInstantiator("vtkContextView", vtkInstantiatorvtkContextViewNew);
  vtkInstantiator::RegisterInstantiator("vtkImageItem", vtkInstantiatorvtkImageItemNew);
  vtkInstantiator::RegisterInstantiator("vtkLookupTableItem", vtkInstantiatorvtkLookupTableItemNew);
  vtkInstantiator::RegisterInstantiator("vtkPen", vtkInstantiatorvtkPenNew);
  vtkInstantiator::RegisterInstantiator("vtkPiecewiseControlPointsItem", vtkInstantiatorvtkPiecewiseControlPointsItemNew);
  vtkInstantiator::RegisterInstantiator("vtkPiecewiseFunctionItem", vtkInstantiatorvtkPiecewiseFunctionItemNew);
  vtkInstantiator::RegisterInstantiator("vtkPiecewisePointHandleItem", vtkInstantiatorvtkPiecewisePointHandleItemNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotBar", vtkInstantiatorvtkPlotBarNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotGrid", vtkInstantiatorvtkPlotGridNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotHistogram2D", vtkInstantiatorvtkPlotHistogram2DNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotLine", vtkInstantiatorvtkPlotLineNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotStacked", vtkInstantiatorvtkPlotStackedNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotParallelCoordinates", vtkInstantiatorvtkPlotParallelCoordinatesNew);
  vtkInstantiator::RegisterInstantiator("vtkPlotPoints", vtkInstantiatorvtkPlotPointsNew);
  vtkInstantiator::RegisterInstantiator("vtkScatterPlotMatrix", vtkInstantiatorvtkScatterPlotMatrixNew);
  vtkInstantiator::RegisterInstantiator("vtkTooltipItem", vtkInstantiatorvtkTooltipItemNew);

  
}
          
void vtkChartsInstantiator::ClassFinalize()
{ 

  vtkInstantiator::UnRegisterInstantiator("vtkAxis", vtkInstantiatorvtkAxisNew);
  vtkInstantiator::UnRegisterInstantiator("vtkAxisExtended", vtkInstantiatorvtkAxisExtendedNew);
  vtkInstantiator::UnRegisterInstantiator("vtkBlockItem", vtkInstantiatorvtkBlockItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkBrush", vtkInstantiatorvtkBrushNew);
  vtkInstantiator::UnRegisterInstantiator("vtkChartLegend", vtkInstantiatorvtkChartLegendNew);
  vtkInstantiator::UnRegisterInstantiator("vtkChartHistogram2D", vtkInstantiatorvtkChartHistogram2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkChartMatrix", vtkInstantiatorvtkChartMatrixNew);
  vtkInstantiator::UnRegisterInstantiator("vtkChartParallelCoordinates", vtkInstantiatorvtkChartParallelCoordinatesNew);
  vtkInstantiator::UnRegisterInstantiator("vtkChartXY", vtkInstantiatorvtkChartXYNew);
  vtkInstantiator::UnRegisterInstantiator("vtkChartPie", vtkInstantiatorvtkChartPieNew);
  vtkInstantiator::UnRegisterInstantiator("vtkColorLegend", vtkInstantiatorvtkColorLegendNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotPie", vtkInstantiatorvtkPlotPieNew);
  vtkInstantiator::UnRegisterInstantiator("vtkColorSeries", vtkInstantiatorvtkColorSeriesNew);
  vtkInstantiator::UnRegisterInstantiator("vtkColorTransferFunctionItem", vtkInstantiatorvtkColorTransferFunctionItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkColorTransferControlPointsItem", vtkInstantiatorvtkColorTransferControlPointsItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkCompositeControlPointsItem", vtkInstantiatorvtkCompositeControlPointsItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkCompositeTransferFunctionItem", vtkInstantiatorvtkCompositeTransferFunctionItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkContext2D", vtkInstantiatorvtkContext2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkContextActor", vtkInstantiatorvtkContextActorNew);
  vtkInstantiator::UnRegisterInstantiator("vtkContextClip", vtkInstantiatorvtkContextClipNew);
  vtkInstantiator::UnRegisterInstantiator("vtkContextInteractorStyle", vtkInstantiatorvtkContextInteractorStyleNew);
  vtkInstantiator::UnRegisterInstantiator("vtkContextScene", vtkInstantiatorvtkContextSceneNew);
  vtkInstantiator::UnRegisterInstantiator("vtkContextTransform", vtkInstantiatorvtkContextTransformNew);
  vtkInstantiator::UnRegisterInstantiator("vtkContextView", vtkInstantiatorvtkContextViewNew);
  vtkInstantiator::UnRegisterInstantiator("vtkImageItem", vtkInstantiatorvtkImageItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkLookupTableItem", vtkInstantiatorvtkLookupTableItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPen", vtkInstantiatorvtkPenNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPiecewiseControlPointsItem", vtkInstantiatorvtkPiecewiseControlPointsItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPiecewiseFunctionItem", vtkInstantiatorvtkPiecewiseFunctionItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPiecewisePointHandleItem", vtkInstantiatorvtkPiecewisePointHandleItemNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotBar", vtkInstantiatorvtkPlotBarNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotGrid", vtkInstantiatorvtkPlotGridNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotHistogram2D", vtkInstantiatorvtkPlotHistogram2DNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotLine", vtkInstantiatorvtkPlotLineNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotStacked", vtkInstantiatorvtkPlotStackedNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotParallelCoordinates", vtkInstantiatorvtkPlotParallelCoordinatesNew);
  vtkInstantiator::UnRegisterInstantiator("vtkPlotPoints", vtkInstantiatorvtkPlotPointsNew);
  vtkInstantiator::UnRegisterInstantiator("vtkScatterPlotMatrix", vtkInstantiatorvtkScatterPlotMatrixNew);
  vtkInstantiator::UnRegisterInstantiator("vtkTooltipItem", vtkInstantiatorvtkTooltipItemNew);

  
}

vtkChartsInstantiator::vtkChartsInstantiator()
{
  if(++vtkChartsInstantiator::Count == 1)
    { 
    vtkChartsInstantiator::ClassInitialize(); 
    }
}

vtkChartsInstantiator::~vtkChartsInstantiator()
{
  if(--vtkChartsInstantiator::Count == 0)
    { 
    vtkChartsInstantiator::ClassFinalize(); 
    }
}

// Number of translation units that include this class's header.
// Purposely not initialized.  Default is static initialization to 0.
unsigned int vtkChartsInstantiator::Count;

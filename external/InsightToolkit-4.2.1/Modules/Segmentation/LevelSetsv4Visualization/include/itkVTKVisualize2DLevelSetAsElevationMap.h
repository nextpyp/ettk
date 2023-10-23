/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef __itkVTKVisualize2DLevelSetAsElevationMap_h
#define __itkVTKVisualize2DLevelSetAsElevationMap_h

#include "itkVTKVisualizeImageLevelSet.h"

#include "itkImageToVTKImageFilter.h"

#include "vnl/vnl_math.h"

#include "vtkImageData.h"
#include "vtkMarchingSquares.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
#include "vtkProperty.h"
#include "vtkTextProperty.h"

namespace itk
{

template< class TInputImage, class TLevelSet >
class VTKVisualize2DLevelSetAsElevationMap :
    public VTKVisualizeImageLevelSet< TInputImage, ImageToVTKImageFilter< TInputImage > >
{
public:
  typedef VTKVisualize2DLevelSetAsElevationMap                                           Self;
  typedef VTKVisualizeImageLevelSet< TInputImage, ImageToVTKImageFilter< TInputImage > > Superclass;
  typedef SmartPointer< Self >                                                           Pointer;
  typedef SmartPointer< const Self >                                                     ConstPointer;

  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro(itkVTKVisualize2DLevelSetAsElevationMap, VTKVisualizeImageLevelSet);

  typedef typename Superclass::InputImageType     InputImageType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::SizeValueType  InputImageSizeValueType;

  typedef TLevelSet                         LevelSetType;
  typedef typename LevelSetType::Pointer    LevelSetPointer;

  void SetLevelSet( LevelSetType * levelSet );

protected:
  VTKVisualize2DLevelSetAsElevationMap();
  virtual ~VTKVisualize2DLevelSetAsElevationMap();

  virtual void PrepareVTKPipeline();

  void GenerateElevationMap();

private:
  VTKVisualize2DLevelSetAsElevationMap( const Self & ); // purposely not implemented
  void operator=( const VTKVisualize2DLevelSetAsElevationMap & ); // purposely not implemented

  LevelSetPointer           m_LevelSet;

  vtkSmartPointer< vtkPolyData >          m_Mesh;
  vtkSmartPointer< vtkScalarBarActor >    m_ScalarBarActor;
  vtkSmartPointer< vtkPolyDataMapper >    m_MeshMapper;
  vtkSmartPointer< vtkActor >             m_SurfaceActor;

  InputImageSizeType m_NumberOfSamples;

  double m_Constant;
  double m_MinValue;
  double m_MaxValue;

  bool   m_ColorValue;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVTKVisualize2DLevelSetAsElevationMap.hxx"
#endif

#endif // itkVTKVisualize2DLevelSetAsElevationMap_H

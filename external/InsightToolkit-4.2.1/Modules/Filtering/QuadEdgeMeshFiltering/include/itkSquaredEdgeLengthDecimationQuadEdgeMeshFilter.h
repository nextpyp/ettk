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
#ifndef __itkSquaredEdgeLengthDecimationQuadEdgeMeshFilter_h
#define __itkSquaredEdgeLengthDecimationQuadEdgeMeshFilter_h

#include "itkEdgeDecimationQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class SquaredEdgeLengthDecimationQuadEdgeMeshFilter
 * \brief
 * \ingroup ITKQuadEdgeMeshFiltering
 */
template< class TInput, class TOutput, class TCriterion >
class ITK_EXPORT SquaredEdgeLengthDecimationQuadEdgeMeshFilter:
  public EdgeDecimationQuadEdgeMeshFilter< TInput, TOutput, TCriterion >
{
public:
  typedef SquaredEdgeLengthDecimationQuadEdgeMeshFilter Self;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;
  typedef EdgeDecimationQuadEdgeMeshFilter<
    TInput, TOutput, TCriterion >                       Superclass;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(SquaredEdgeLengthDecimationQuadEdgeMeshFilter, EdgeDecimationQuadEdgeMeshFilter);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  typedef TInput                          InputMeshType;
  typedef typename InputMeshType::Pointer InputMeshPointer;

  typedef TOutput                                         OutputMeshType;
  typedef typename OutputMeshType::Pointer                OutputMeshPointer;
  typedef typename OutputMeshType::PointIdentifier        OutputPointIdentifier;
  typedef typename OutputMeshType::PointType              OutputPointType;
  typedef typename OutputMeshType::QEType                 OutputQEType;
  typedef typename OutputMeshType::EdgeCellType           OutputEdgeCellType;
  typedef typename OutputMeshType::CellsContainerIterator OutputCellsContainerIterator;

  typedef TCriterion                          CriterionType;
  typedef typename CriterionType::MeasureType MeasureType;

  typedef typename Superclass::PriorityType          PriorityType;
  typedef typename Superclass::PriorityQueueItemType PriorityQueueItemType;
  typedef typename Superclass::PriorityQueueType     PriorityQueueType;
  typedef typename Superclass::PriorityQueuePointer  PriorityQueuePointer;

  typedef typename Superclass::QueueMapType     QueueMapType;
  typedef typename Superclass::QueueMapIterator QueueMapIterator;

  typedef typename Superclass::OperatorType    OperatorType;
  typedef typename Superclass::OperatorPointer OperatorPointer;
protected:

  SquaredEdgeLengthDecimationQuadEdgeMeshFilter();
  virtual ~SquaredEdgeLengthDecimationQuadEdgeMeshFilter();

  /**
   * \brief Compute the measure value for iEdge
   * \param[in] iEdge
   * \return measure value, here the squared edge length
   */
  inline MeasureType MeasureEdge(OutputQEType *iEdge)
  {
    OutputPointIdentifier id_org = iEdge->GetOrigin();
    OutputPointIdentifier id_dest = iEdge->GetDestination();

    OutputPointType org = this->m_OutputMesh->GetPoint(id_org);
    OutputPointType dest = this->m_OutputMesh->GetPoint(id_dest);

    return static_cast< MeasureType >( org.SquaredEuclideanDistanceTo(dest) );
  }

  /**
   * \param[in] iEdge
   * \return the optimal point location
   */
  OutputPointType Relocate(OutputQEType *iEdge);

private:
  SquaredEdgeLengthDecimationQuadEdgeMeshFilter(const Self &);
  void operator=(const Self &);
};
}

#include "itkSquaredEdgeLengthDecimationQuadEdgeMeshFilter.hxx"
#endif

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
#ifndef __itkJointHistogramMutualInformationGetValueAndDerivativeThreader_hxx
#define __itkJointHistogramMutualInformationGetValueAndDerivativeThreader_hxx

#include "itkJointHistogramMutualInformationGetValueAndDerivativeThreader.h"

namespace itk
{

template< class TDomainPartitioner, class TImageToImageMetric, class TJointHistogramMetric >
void
JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >
::BeforeThreadedExecution()
{
  Superclass::BeforeThreadedExecution();

  TJointHistogramMetric * associate = dynamic_cast< TJointHistogramMetric * >( this->m_Associate );

  this->m_JointPDFInterpolatorPerThread.resize( this->GetNumberOfThreadsUsed() );
  this->m_FixedImageMarginalPDFInterpolatorPerThread.resize( this->GetNumberOfThreadsUsed() );
  this->m_MovingImageMarginalPDFInterpolatorPerThread.resize( this->GetNumberOfThreadsUsed() );

  for( ThreadIdType i = 0; i < this->GetNumberOfThreadsUsed(); ++i )
    {
    if( this->m_JointPDFInterpolatorPerThread[i].IsNull() )
      {
      this->m_JointPDFInterpolatorPerThread[i] = JointPDFInterpolatorType::New();
      }
    this->m_JointPDFInterpolatorPerThread[i]->SetInputImage( associate->m_JointPDF );
    if( this->m_FixedImageMarginalPDFInterpolatorPerThread[i].IsNull() )
      {
      this->m_FixedImageMarginalPDFInterpolatorPerThread[i] = MarginalPDFInterpolatorType::New();
      }
    this->m_FixedImageMarginalPDFInterpolatorPerThread[i]->SetInputImage( associate->m_FixedImageMarginalPDF );
    if( this->m_MovingImageMarginalPDFInterpolatorPerThread[i].IsNull() )
      {
      this->m_MovingImageMarginalPDFInterpolatorPerThread[i] = MarginalPDFInterpolatorType::New();
      }
    this->m_MovingImageMarginalPDFInterpolatorPerThread[i]->SetInputImage( associate->m_MovingImageMarginalPDF );
    }
}

template< class TDomainPartitioner, class TImageToImageMetric, class TJointHistogramMetric >
void
JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >
::AfterThreadedExecution()
{
  Superclass::AfterThreadedExecution();

  TJointHistogramMetric * associate = dynamic_cast< TJointHistogramMetric * >( this->m_Associate );
  // The Superclass does not generate a valid m_Value for this metric.  We have to calculate it
  // here, but only if there are 1 or more valid points. Otherwise the Superclass
  // will have already set a default value and issued a warning.
  if( associate->GetNumberOfValidPoints() > 0 )
    {
    associate->m_Value = associate->ComputeValue();
    }
}

template< class TDomainPartitioner, class TImageToImageMetric, class TJointHistogramMetric >
bool
JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >
::ProcessPoint( const VirtualIndexType &,
                const VirtualPointType &        virtualPoint,
                const FixedImagePointType &,
                const FixedImagePixelType &     fixedImageValue,
                const FixedImageGradientType &,
                const MovingImagePointType &,
                const MovingImagePixelType &    movingImageValue,
                const MovingImageGradientType & movingImageGradient,
                MeasureType &,
                DerivativeType &                localDerivativeReturn,
                const ThreadIdType              threadId ) const
{
  TJointHistogramMetric * associate = dynamic_cast< TJointHistogramMetric * >( this->m_Associate );

  // check that the moving image sample is within the range of the true min
  // and max, hence being within the moving image mask
  if ( movingImageValue < associate->m_MovingImageTrueMin )
    {
    return false;
    }
  else if ( movingImageValue > associate->m_MovingImageTrueMax )
    {
    return false;
    }
  /** the scalingfactor is the MI specific scaling of the image gradient and jacobian terms */
  InternalComputationValueType scalingfactor = NumericTraits< InternalComputationValueType >::Zero; // for scaling the jacobian terms

  JointPDFPointType jointPDFpoint;
  associate->ComputeJointPDFPoint( fixedImageValue, movingImageValue, jointPDFpoint );
  // Make sure the point is inside th joint pdf.
  if ( ! this->m_JointPDFInterpolatorPerThread[threadId]->IsInsideBuffer( jointPDFpoint ) )
    {
    return false;
    }
  InternalComputationValueType jointPDFValue = this->m_JointPDFInterpolatorPerThread[threadId]->Evaluate( jointPDFpoint );
  SizeValueType ind = 1;
  InternalComputationValueType dJPDF = this->ComputeJointPDFDerivative( jointPDFpoint, threadId , ind );
  typename MarginalPDFType::PointType mind;
  mind[0] = jointPDFpoint[ind];
  InternalComputationValueType movingImagePDFValue =
    this->m_MovingImageMarginalPDFInterpolatorPerThread[threadId]->Evaluate(mind);
  InternalComputationValueType dMmPDF =
    this->ComputeMovingImageMarginalPDFDerivative( mind , threadId );

  InternalComputationValueType term1 = NumericTraits< InternalComputationValueType >::Zero;
  InternalComputationValueType term2 = NumericTraits< InternalComputationValueType >::Zero;
  InternalComputationValueType eps = 1.e-16;
  if( jointPDFValue > eps &&  movingImagePDFValue > eps )
    {
    const InternalComputationValueType pRatio =
                            vcl_log(jointPDFValue)-vcl_log(movingImagePDFValue);
    term1 = dJPDF*pRatio;
    term2 = associate->m_Log2 * dMmPDF * jointPDFValue / movingImagePDFValue;
    scalingfactor =  ( term2 - term1 );
    }  // end if-block to check non-zero bin contribution
  else
    {
    scalingfactor = NumericTraits< InternalComputationValueType >::Zero;
    }

  /* Use a pre-allocated jacobian object for efficiency */
  FixedTransformJacobianType & jacobian =
    const_cast< FixedTransformJacobianType &   >(this->m_MovingTransformJacobianPerThread[threadId]);

  /** For dense transforms, this returns identity */
  associate->m_MovingTransform->ComputeJacobianWithRespectToParameters( virtualPoint, jacobian );

  for ( NumberOfParametersType par = 0; par < associate->GetNumberOfLocalParameters(); par++ )
    {
    InternalComputationValueType sum = NumericTraits< InternalComputationValueType >::Zero;
    for ( SizeValueType dim = 0; dim < TImageToImageMetric::MovingImageDimension; dim++ )
      {
      sum += scalingfactor * jacobian(dim, par) * movingImageGradient[dim];
      }
    localDerivativeReturn[par] = sum;
    }
  return true;
}

template< class TDomainPartitioner, class TImageToImageMetric, class TJointHistogramMetric >
typename JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >::InternalComputationValueType
JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >
::ComputeFixedImageMarginalPDFDerivative( const MarginalPDFPointType & margPDFpoint,
                                          const ThreadIdType threadID ) const
{
  TJointHistogramMetric * associate = dynamic_cast< TJointHistogramMetric * >( this->m_Associate );

  InternalComputationValueType offset = 0.5*this->m_JointPDFSpacing[0];
  InternalComputationValueType eps = this->m_JointPDFSpacing[0];
  MarginalPDFPointType         leftpoint = margPDFpoint;
  leftpoint[0] -= offset;
  MarginalPDFPointType  rightpoint = margPDFpoint;
  rightpoint[0] += offset;
  if (leftpoint[0] < eps )
    {
    leftpoint[0] = eps;
    }
  if (rightpoint[0] < eps )
    {
    rightpoint[0] = eps;
    }
  if (leftpoint[0] > 1.0 )
    {
    leftpoint[0] = 1.0;
    }
  if (rightpoint[0] > 1.0 )
    {
    rightpoint[0] = 1.0;
    }
  InternalComputationValueType delta = rightpoint[0]-leftpoint[0];
  if ( delta > NumericTraits< InternalComputationValueType >::Zero )
    {
    InternalComputationValueType deriv = this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->Evaluate(rightpoint) -
      this->m_ThreaderFixedImageMarginalPDFInterpolator[threadID]->Evaluate(leftpoint);
    return deriv/delta;
    }
  else
    {
    return NumericTraits< InternalComputationValueType >::Zero;
    }
}

template< class TDomainPartitioner, class TImageToImageMetric, class TJointHistogramMetric >
typename JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >::InternalComputationValueType
JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >
::ComputeMovingImageMarginalPDFDerivative( const MarginalPDFPointType & margPDFpoint,
                                           const ThreadIdType threadId ) const
{
  TJointHistogramMetric * associate = dynamic_cast< TJointHistogramMetric * >( this->m_Associate );

  InternalComputationValueType offset = 0.5*associate->m_JointPDFSpacing[0];
  InternalComputationValueType eps = associate->m_JointPDFSpacing[0];
  MarginalPDFPointType  leftpoint = margPDFpoint;
  leftpoint[0] -= offset;
  MarginalPDFPointType  rightpoint = margPDFpoint;
  rightpoint[0] += offset;
  if( leftpoint[0] < eps )
    {
    leftpoint[0] = eps;
    }
  if( rightpoint[0] < eps )
    {
    rightpoint[0] = eps;
    }
  if( leftpoint[0] > 1.0 )
    {
    leftpoint[0] = 1.0;
    }
  if( rightpoint[0] > 1.0  )
    {
    rightpoint[0] = 1.0;
    }
  InternalComputationValueType delta = rightpoint[0] - leftpoint[0];
  if ( delta > NumericTraits< InternalComputationValueType >::Zero )
    {
    InternalComputationValueType deriv =
      this->m_MovingImageMarginalPDFInterpolatorPerThread[threadId]->Evaluate(rightpoint) -
      this->m_MovingImageMarginalPDFInterpolatorPerThread[threadId]->Evaluate(leftpoint);
    return deriv/delta;
    }
  else
    {
    return NumericTraits< InternalComputationValueType >::Zero;
    }
}

template< class TDomainPartitioner, class TImageToImageMetric, class TJointHistogramMetric >
typename JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >::InternalComputationValueType
JointHistogramMutualInformationGetValueAndDerivativeThreader< TDomainPartitioner, TImageToImageMetric, TJointHistogramMetric >
::ComputeJointPDFDerivative( const JointPDFPointType & jointPDFpoint,
                             const ThreadIdType threadId,
                             const SizeValueType ind ) const
{
  TJointHistogramMetric * associate = dynamic_cast< TJointHistogramMetric * >( this->m_Associate );

  InternalComputationValueType offset = 0.5*associate->m_JointPDFSpacing[ind];
  InternalComputationValueType eps = associate->m_JointPDFSpacing[ind];
  JointPDFPointType  leftpoint = jointPDFpoint;
  leftpoint[ind] -= offset;
  JointPDFPointType  rightpoint = jointPDFpoint;
  rightpoint[ind] += offset;

  if (leftpoint[ind] < eps )
    {
    leftpoint[ind] = eps;
    }

  if (rightpoint[ind] < eps )
    {
    rightpoint[ind] = eps;
    }

  if (leftpoint[ind] > 1.0 )
    {
    leftpoint[ind] = 1.0;
    }

  if (rightpoint[ind] > 1.0 )
    {
    rightpoint[ind] = 1.0;
    }

  InternalComputationValueType delta = rightpoint[ind] - leftpoint[ind];
  InternalComputationValueType deriv = NumericTraits< InternalComputationValueType >::Zero;
  if ( delta > NumericTraits< InternalComputationValueType >::Zero )
    {
    deriv = this->m_JointPDFInterpolatorPerThread[threadId]->Evaluate(rightpoint)-
          this->m_JointPDFInterpolatorPerThread[threadId]->Evaluate(leftpoint);
    return deriv/delta;
    }
  return deriv;
}

} // end namespace itk

#endif

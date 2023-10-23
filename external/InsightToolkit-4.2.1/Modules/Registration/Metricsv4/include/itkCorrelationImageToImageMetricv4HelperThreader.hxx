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
#ifndef __itkCorrelationImageToImageMetricv4HelperThreader_hxx
#define __itkCorrelationImageToImageMetricv4HelperThreader_hxx

#include "itkCorrelationImageToImageMetricv4HelperThreader.h"

namespace itk
{
template<class TDomainPartitioner, class TImageToImageMetric, class TCorrelationMetric>
void CorrelationImageToImageMetricv4HelperThreader< TDomainPartitioner, TImageToImageMetric, TCorrelationMetric>
::BeforeThreadedExecution()
{
   Superclass::BeforeThreadedExecution();

   this->m_FixSumPerThread.resize(this->GetNumberOfThreadsUsed());
   this->m_MovSumPerThread.resize(this->GetNumberOfThreadsUsed());

    //---------------------------------------------------------------
    // Set initial values.
    for (ThreadIdType i = 0; i < this->GetNumberOfThreadsUsed(); i++)
      {
      this->m_FixSumPerThread[i] = NumericTraits<InternalComputationValueType>::Zero;
      this->m_MovSumPerThread[i] = NumericTraits<InternalComputationValueType>::Zero;
      }

}

template<class TDomainPartitioner, class TImageToImageMetric,
  class TCorrelationMetric>
void
CorrelationImageToImageMetricv4HelperThreader<TDomainPartitioner,
    TImageToImageMetric, TCorrelationMetric>::AfterThreadedExecution()
{

  // only need to dynamic cast asscociate to the specific metric class
  // to  access the inherited protected members from its parent class
  TCorrelationMetric * associate =
      dynamic_cast<TCorrelationMetric *>(this->m_Associate);

  /* Store the number of valid points the enclosing class \c
   * m_NumberOfValidPoints by collecting the valid points per thread. */
  associate->m_NumberOfValidPoints =
      NumericTraits<SizeValueType>::Zero;
  for (ThreadIdType i = 0; i < this->GetNumberOfThreadsUsed(); i++)
    {
    associate->m_NumberOfValidPoints += this->m_NumberOfValidPointsPerThread[i];
    }

  if (associate->m_NumberOfValidPoints <= 0 )
    {
    itkWarningMacro("collected only zero points");
    return;
    }

  InternalComputationValueType sumF = NumericTraits<InternalComputationValueType>::Zero;
  InternalComputationValueType sumM = NumericTraits<InternalComputationValueType>::Zero;

  for (size_t i = 0; i < this->m_MeasurePerThread.size(); i++)
    {
    sumF += this->m_FixSumPerThread[i];
    sumM += this->m_MovSumPerThread[i];
    }

  associate->m_AverageFix = sumF / associate->m_NumberOfValidPoints;
  associate->m_AverageMov = sumM / associate->m_NumberOfValidPoints;
}

template<class TDomainPartitioner, class TImageToImageMetric,
class TCorrelationMetric>
bool
CorrelationImageToImageMetricv4HelperThreader<TDomainPartitioner,
TImageToImageMetric, TCorrelationMetric>
::ProcessVirtualPoint( const VirtualIndexType & virtualIndex,
                       const VirtualPointType & virtualPoint,
                       const ThreadIdType threadID )
{
  FixedOutputPointType        mappedFixedPoint;
  FixedImagePixelType         mappedFixedPixelValue;
  FixedImageGradientType      mappedFixedImageGradient;
  MovingOutputPointType       mappedMovingPoint;
  MovingImagePixelType        mappedMovingPixelValue;
  MovingImageGradientType     mappedMovingImageGradient;
  bool                        pointIsValid = false;


  TCorrelationMetric * associate =
        dynamic_cast<TCorrelationMetric *>(this->m_Associate);

  /* Transform the point into fixed and moving spaces, and evaluate.
   * Different behavior with pre-warping enabled is handled transparently.
   * Do this in a try block to catch exceptions and print more useful info
   * then we otherwise get when exceptions are caught in MultiThreader. */
  try
    {
    pointIsValid = associate->TransformAndEvaluateFixedPoint( virtualIndex,
                                      virtualPoint,
                                      false,
                                      mappedFixedPoint,
                                      mappedFixedPixelValue,
                                      mappedFixedImageGradient );
    }
  catch( ExceptionObject & exc )
    {
    //NOTE: there must be a cleaner way to do this:
    std::string msg("Caught exception: \n");
    msg += exc.what();
    ExceptionObject err(__FILE__, __LINE__, msg);
    throw err;
    }
  if( !pointIsValid )
    {
    return pointIsValid;
    }

  try
    {
    pointIsValid = associate->TransformAndEvaluateMovingPoint( virtualIndex,
                                    virtualPoint,
                                    false,
                                    mappedMovingPoint,
                                    mappedMovingPixelValue,
                                    mappedMovingImageGradient );
    }
  catch( ExceptionObject & exc )
    {
    std::string msg("Caught exception: \n");
    msg += exc.what();
    ExceptionObject err(__FILE__, __LINE__, msg);
    throw err;
    }
  if( !pointIsValid )
    {
    return pointIsValid;
    }

  /* Do the specific calculations for values */
  try
    {
      this->m_FixSumPerThread[threadID] += mappedFixedPixelValue;
      this->m_MovSumPerThread[threadID] += mappedMovingPixelValue;
    }
  catch( ExceptionObject & exc )
    {
    //NOTE: there must be a cleaner way to do this:
    std::string msg("Exception in ProcessVirtualPoint:\n");
    msg += exc.what();
    ExceptionObject err(__FILE__, __LINE__, msg);
    throw err;
    }
  if( pointIsValid )
    {
    this->m_NumberOfValidPointsPerThread[threadID]++;
    }

  return pointIsValid;
}
} // end namespace itk

#endif

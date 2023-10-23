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
#ifndef __itkOtsuThresholdCalculator_hxx
#define __itkOtsuThresholdCalculator_hxx

#include "itkOtsuThresholdCalculator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_math.h"

namespace itk
{
template< class THistogram, class TOutput >
void
OtsuThresholdCalculator< THistogram, TOutput >
::GenerateData(void)
{
  const HistogramType * histogram = this->GetInput();
  // histogram->Print(std::cout);
  if ( histogram->GetTotalFrequency() == 0 )
    {
    itkExceptionMacro(<< "Histogram is empty");
    }
  ProgressReporter progress(this, 0, histogram->GetSize(0)*2 );
  if( histogram->GetSize(0) == 1 )
    {
    this->GetOutput()->Set( static_cast<OutputType>(histogram->GetMeasurement(0,0)) );
    }

  unsigned int size = histogram->GetSize(0);

  // normalize the frequencies
  std::vector< double > relativeFrequency;
  relativeFrequency.resize(size);
  double totalMean = 0.0;
  for ( unsigned int j = 0; j < size; j++ )
    {
    relativeFrequency[j] = histogram->GetFrequency(j,0);
    relativeFrequency[j] /= histogram->GetTotalFrequency();
    totalMean += ( j + 1 ) * relativeFrequency[j];
    progress.CompletedPixel();
    }

  // compute Otsu's threshold by maximizing the between-class
  // variance
  double freqLeft = relativeFrequency[0];
  double meanLeft = 1.0;
  double meanRight = ( totalMean - freqLeft ) / ( 1.0 - freqLeft );

  double maxVarBetween = freqLeft * ( 1.0 - freqLeft )
                         * vnl_math_sqr(meanLeft - meanRight);
  int maxBinNumber = 0;

  double freqLeftOld = freqLeft;
  double meanLeftOld = meanLeft;

  for ( unsigned int j = 1; j < size; j++ )
    {
    freqLeft += relativeFrequency[j];
    meanLeft = ( meanLeftOld * freqLeftOld
                 + ( j + 1 ) * relativeFrequency[j] ) / freqLeft;
    if ( freqLeft == 1.0 )
      {
      meanRight = 0.0;
      }
    else
      {
      meanRight = ( totalMean - meanLeft * freqLeft )
                  / ( 1.0 - freqLeft );
      }
    double varBetween = freqLeft * ( 1.0 - freqLeft )
                        * vnl_math_sqr(meanLeft - meanRight);

    // for portability - different compilers seem to produce results
    // that differ by one bin, presumable because of the precision of
    // calculating varBetween
    const double tolerance = 0.00001;
    if ( (varBetween - tolerance) > maxVarBetween )
      {
      maxVarBetween = varBetween;
      maxBinNumber = j;
      }

    // cache old values
    freqLeftOld = freqLeft;
    meanLeftOld = meanLeft;
    progress.CompletedPixel();
    }
  // should be this for backward compatibility
  // this->GetOutput()->Set( static_cast<OutputType>(
  // histogram->GetBinMin( 0, maxBinNumber + 1 ) ) );
  this->GetOutput()->Set( static_cast<OutputType>( histogram->GetMeasurement( maxBinNumber + 1, 0 ) ) );
}

} // end namespace itk

#endif

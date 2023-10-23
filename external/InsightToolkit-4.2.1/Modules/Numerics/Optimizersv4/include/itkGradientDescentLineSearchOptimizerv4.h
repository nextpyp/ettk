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
#ifndef __itkGradientDescentLineSearchOptimizerv4_h
#define __itkGradientDescentLineSearchOptimizerv4_h

#include "itkGradientDescentOptimizerv4.h"
#include "itkOptimizerParameterScalesEstimator.h"
#include "itkWindowConvergenceMonitoringFunction.h"

namespace itk
{
/** \class GradientDescentLineSearchOptimizerv4
 *  \brief Gradient descent optimizer with a golden section line search.
 *
 * GradientDescentLineSearchOptimizer implements a simple gradient descent optimizer
 * that is followed by a line search to find the best value for the learning rate.
 * At each iteration the current position is updated according to
 *
 * \f[
 *        p_{n+1} = p_n
 *                + \mbox{learningRateByGoldenSectionLineSearch}
                  \, \frac{\partial f(p_n) }{\partial p_n}
 * \f]
 *
 * Options are identical to the superclass's except for:
 *
 * options Epsilon, LowerLimit and UpperLimit that will guide
 * a golden section line search to find the optimal gradient update
 * within the range :
 *
 *   [ learningRate * LowerLimit , learningRate * UpperLimit ]
 *
 * where Epsilon sets the resolution of the search.  Smaller values
 * lead to additional computation time but better localization of
 * the minimum.
 *
 * By default, this optimizer will return the best value and associated
 * parameters that were calculated during the optimization.
 * See SetReturnBestParametersAndValue().
 *
 * \ingroup ITKOptimizersv4
 */
class ITK_EXPORT GradientDescentLineSearchOptimizerv4
  : public GradientDescentOptimizerv4
{
public:
  /** Standard class typedefs. */
  typedef GradientDescentLineSearchOptimizerv4  Self;
  typedef GradientDescentOptimizerv4            Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GradientDescentLineSearchOptimizerv4, GradientDescentOptimizerv4);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Derivative type */
  typedef Superclass::DerivativeType      DerivativeType;

  /** Metric type over which this class is templated */
  typedef Superclass::MeasureType                  MeasureType;
  typedef Superclass::InternalComputationValueType InternalComputationValueType;

  /** Type for the convergence checker */
  typedef itk::Function::WindowConvergenceMonitoringFunction<double> ConvergenceMonitoringType;

  /** The epsilon determines the accuracy of the line search
   *  i.e. the energy alteration that is considered convergent.
   */
  itkSetMacro( Epsilon , InternalComputationValueType );
  itkGetMacro( Epsilon , InternalComputationValueType );

  /** The upper and lower limit below determine the range
   *  of values over which the learning rate can be adjusted
   *  by the golden section line search.  The update can then
   *  occur in the range from the smallest change given by :
   *     NewParams = OldParams + LowerLimit * gradient
   *  to the largest change given by :
   *     NewParams = OldParams + UpperLimit * gradient
   *  Reasonable values might be 0 and 2.
   */
  itkSetMacro( LowerLimit , InternalComputationValueType );
  itkGetMacro( LowerLimit , InternalComputationValueType );
  itkSetMacro( UpperLimit , InternalComputationValueType );
  itkGetMacro( UpperLimit , InternalComputationValueType );
  itkSetMacro( MaximumLineSearchIterations , unsigned int );
  itkGetMacro( MaximumLineSearchIterations , unsigned int );

protected:

  /** Advance one Step following the gradient direction.
   * Includes transform update. */
  virtual void AdvanceOneStep(void);

  /** Default constructor */
  GradientDescentLineSearchOptimizerv4();

  /** Destructor */
  virtual ~GradientDescentLineSearchOptimizerv4();

  virtual void PrintSelf( std::ostream & os, Indent indent ) const;

  InternalComputationValueType GoldenSectionSearch( InternalComputationValueType a, InternalComputationValueType b, InternalComputationValueType c );

  InternalComputationValueType m_LowerLimit;
  InternalComputationValueType m_UpperLimit;
  InternalComputationValueType m_Phi;
  InternalComputationValueType m_Resphi;
  InternalComputationValueType m_Epsilon;

  /** Controls the maximum recursion depth for the golden section search */
  unsigned int      m_MaximumLineSearchIterations;
  /** Counts the recursion depth for the golden section search */
  unsigned int      m_LineSearchIterations;

private:

  GradientDescentLineSearchOptimizerv4( const Self & ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

};

} // end namespace itk

#endif

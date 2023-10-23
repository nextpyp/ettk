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
#ifndef __itkMultiGradientOptimizerv4_h
#define __itkMultiGradientOptimizerv4_h
#include "itkObjectToObjectOptimizerBase.h"
#include "itkGradientDescentOptimizerv4.h"

namespace itk
{
/** \class MultiGradientOptimizerv4
 *  \brief Multiple gradient-based optimizers are combined in order to perform a multi-objective optimization.
 *
 *  This optimizer will do a combined gradient descent optimization using whatever metric/optimizer gradient
 *  sub-optimizers are passed to it by the user.  The learning rate or scaleestimator for each sub-optimizer
 *  controls the relative weight of each metric in the optimization.  Denote the weights as \f$ w_1 \f$ and \f$ w_2 \f$ then
 *  the MultiGradientOptimizer will optimize \f$ \sum_i w_i Metric_i \f$ by using update rule:
 *
 *  \f[
 *    params_{new} = params_{old} + \frac{1}{N_{Metrics}} * ( \sum_i w_i Grad(Metric_i) )
 *  \f]
 *
 *  The test for this class illustrates the expected behavior.
 *
 * \ingroup ITKOptimizersv4
 */

class ITK_EXPORT MultiGradientOptimizerv4
  : public GradientDescentOptimizerv4
{
public:
  /** Standard class typedefs. */
  typedef MultiGradientOptimizerv4       Self;
  typedef GradientDescentOptimizerv4     Superclass;
  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiGradientOptimizerv4, GradientDescentOptimizerv4);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef itk::GradientDescentOptimizerv4            LocalOptimizerType;
  typedef itk::GradientDescentOptimizerv4::Pointer   LocalOptimizerPointer;
  typedef Superclass::ParametersType                 ParametersType;
  typedef ObjectToObjectOptimizerBase                OptimizerType;
  typedef OptimizerType::Pointer                     OptimizerPointer;
  typedef std::vector< LocalOptimizerPointer >       OptimizersListType;
  typedef OptimizersListType::size_type              OptimizersListSizeType;

  /** Stop condition return string type */
  typedef std::string                            StopConditionReturnStringType;

  /** Stop condition internal string type */
  typedef std::ostringstream                     StopConditionDescriptionType;

  /** Metric type over which this class is templated */
  typedef Superclass::MetricType                    MetricType;
  typedef MetricType::Pointer                       MetricTypePointer;

  /** Derivative type */
  typedef MetricType::DerivativeType                DerivativeType;

  /** Measure type */
  typedef Superclass::MeasureType                   MeasureType;
  typedef std::vector< MeasureType >                MetricValuesListType;

  /** Internal computation type, for maintaining a desired precision */
  typedef Superclass::InternalComputationValueType InternalComputationValueType;

  /** Get stop condition enum */
  itkGetConstReferenceMacro(StopCondition, StopConditionType);

  /** Set the number of iterations. */
  itkSetMacro(NumberOfIterations, SizeValueType);

  /** Get the number of iterations. */
  itkGetConstReferenceMacro(NumberOfIterations, SizeValueType);

  /** Get the current iteration number. */
  itkGetConstMacro(CurrentIteration, SizeValueType);

  /** Begin the optimization */
  virtual void StartOptimization(void);

  /** Stop optimization. The object is left in a state so the
   * optimization can be resumed by calling ResumeOptimization. */
  virtual void StopOptimization(void);

  /** Resume the optimization. Can be called after StopOptimization to
   * resume. The bulk of the optimization work loop is here. */
  virtual void ResumeOptimization();

  /** Get the reason for termination */
  virtual const StopConditionReturnStringType GetStopConditionDescription() const;

  /** Get the list of optimizers currently held.  */
  OptimizersListType & GetOptimizersList();

  /** Set the list of optimizers to combine */
  void SetOptimizersList(OptimizersListType & p);

  /** Get the list of metric values that we produced after the multi-objective search.  */
  const MetricValuesListType & GetMetricValuesList() const;

protected:

  /** Default constructor */
  MultiGradientOptimizerv4();
  virtual ~MultiGradientOptimizerv4();

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

  /* Common variables for optimization control and reporting */
  bool                          m_Stop;
  StopConditionType             m_StopCondition;
  StopConditionDescriptionType  m_StopConditionDescription;
  SizeValueType                 m_NumberOfIterations;
  SizeValueType                 m_CurrentIteration;
  OptimizersListType            m_OptimizersList;
  MetricValuesListType          m_MetricValuesList;
  MeasureType                   m_MinimumMetricValue;
  MeasureType                   m_MaximumMetricValue;

private:
  MultiGradientOptimizerv4( const Self & ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

};

} // end namespace itk

#endif

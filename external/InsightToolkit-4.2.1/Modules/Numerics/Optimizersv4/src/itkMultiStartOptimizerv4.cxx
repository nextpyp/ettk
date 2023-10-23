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

#include "itkMultiStartOptimizerv4.h"

namespace itk
{

//-------------------------------------------------------------------
MultiStartOptimizerv4
::MultiStartOptimizerv4()
{
  this->m_NumberOfIterations = static_cast<SizeValueType>(0);
  this->m_CurrentIteration   = static_cast<SizeValueType>(0);
  this->m_StopCondition      = MAXIMUM_NUMBER_OF_ITERATIONS;
  this->m_StopConditionDescription << this->GetNameOfClass() << ": ";

  this->m_BestParametersIndex= static_cast<ParameterListSizeType>(0);
  this->m_MaximumMetricValue=NumericTraits<MeasureType>::max();
  this->m_MinimumMetricValue = this->m_MaximumMetricValue;
  m_LocalOptimizer = NULL;
}

//-------------------------------------------------------------------
MultiStartOptimizerv4
::~MultiStartOptimizerv4()
{}

//-------------------------------------------------------------------
void
MultiStartOptimizerv4
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Number of iterations: " << this->m_NumberOfIterations  << std::endl;
  os << indent << "Current iteration: " << this->m_CurrentIteration << std::endl;
  os << indent << "Stop condition:"<< this->m_StopCondition << std::endl;
  os << indent << "Stop condition description: " << this->m_StopConditionDescription.str()  << std::endl;
}

//-------------------------------------------------------------------
MultiStartOptimizerv4::ParametersListType &
MultiStartOptimizerv4
::GetParametersList()
{
  return this->m_ParametersList;
}

/** Set the list of parameters over which to search */
void
MultiStartOptimizerv4
::SetParametersList(MultiStartOptimizerv4::ParametersListType & p)
{
  if( p != this->m_ParametersList )
    {
    this->m_ParametersList = p;
    this->Modified();
    }
}

/** Get the list of metric values that we produced after the multi-start search.  */
const MultiStartOptimizerv4::MetricValuesListType &
MultiStartOptimizerv4
::GetMetricValuesList() const
{
  return this->m_MetricValuesList;
}

//-------------------------------------------------------------------
MultiStartOptimizerv4::ParametersType
MultiStartOptimizerv4
::GetBestParameters( )
{
  return this->m_ParametersList[m_BestParametersIndex];
}


//-------------------------------------------------------------------
void
MultiStartOptimizerv4
::InstantiateLocalOptimizer()
{
  LocalOptimizerPointer optimizer = LocalOptimizerType::New();
  typedef LocalOptimizerType::InternalComputationValueType LearningType;
  optimizer->SetLearningRate( static_cast<LearningType>(1.e-1) );
  optimizer->SetNumberOfIterations( 25 );
  this->m_LocalOptimizer=optimizer;
}

//-------------------------------------------------------------------
const MultiStartOptimizerv4::StopConditionReturnStringType
MultiStartOptimizerv4
::GetStopConditionDescription() const
{
  return this->m_StopConditionDescription.str();
}

//-------------------------------------------------------------------
void
MultiStartOptimizerv4
::StopOptimization(void)
{
  itkDebugMacro( "StopOptimization called with a description - "
    << this->GetStopConditionDescription() );
  this->m_Stop = true;

  this->m_Metric->SetParameters( this->m_ParametersList[ this->m_BestParametersIndex ] );

  this->InvokeEvent( EndEvent() );
}

/**
 * Start and run the optimization
 */
void
MultiStartOptimizerv4
::StartOptimization()
{
  itkDebugMacro("StartOptimization");

  this->m_NumberOfIterations=this->m_ParametersList.size();
  this->m_MetricValuesList.clear();
  this->m_BestParametersIndex = static_cast<ParameterListSizeType>(0);
  this->m_MinimumMetricValue = this->m_MaximumMetricValue;

  /* Must call the superclass version for basic validation and setup */
  if ( this->m_NumberOfIterations > static_cast<SizeValueType>(0) )
    {
    Superclass::StartOptimization();
    }

  this->m_CurrentIteration = static_cast<SizeValueType>(0);

  if ( this->m_NumberOfIterations > static_cast<SizeValueType>(0) )
    {
    this->ResumeOptimization();
    }
}

/**
 * Resume optimization.
 */
void
MultiStartOptimizerv4
::ResumeOptimization()
{
  this->m_StopConditionDescription.str("");
  this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
  this->InvokeEvent( StartEvent() );

  this->m_Stop = false;
  while( ! this->m_Stop )
    {
    /* Compute metric value */
    try
      {
      this->m_Metric->SetParameters( this->m_ParametersList[ this->m_CurrentIteration ] );
      if (  this->m_LocalOptimizer )
        {
        this->m_LocalOptimizer->SetMetric( this->m_Metric );
        this->m_LocalOptimizer->StartOptimization();
        this->m_ParametersList[this->m_CurrentIteration] = this->m_Metric->GetParameters();
        }
      this->m_CurrentMetricValue = this->m_Metric->GetValue();
      this->m_MetricValuesList.push_back(this->m_CurrentMetricValue);
      }
    catch ( ExceptionObject & )
      {
      /** We simply ignore this exception because it may just be a bad starting point.
       *  We hope that other start points are better.
       */
      itkWarningMacro("An exception occurred in sub-optimization number " << this->m_CurrentIteration << ".  If too many of these occur, you may need to set a different set of initial parameters.");
      }

    if ( this->m_CurrentMetricValue <  this->m_MinimumMetricValue )
      {
      this->m_MinimumMetricValue=this->m_CurrentMetricValue;
      this->m_BestParametersIndex = this->m_CurrentIteration;
      }
    /* Check if optimization has been stopped externally.
     * (Presumably this could happen from a multi-threaded client app?) */
    if ( this->m_Stop )
      {
      this->m_StopConditionDescription << "StopOptimization() called";
      break;
      }

    this->InvokeEvent( IterationEvent() );

    /* Update and check iteration count */
    this->m_CurrentIteration++;
    if ( this->m_CurrentIteration >= this->m_NumberOfIterations )
      {
      this->m_StopConditionDescription << "Maximum number of iterations ("
                                 << this->m_NumberOfIterations
                                 << ") exceeded.";
      this->m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
      this->StopOptimization();
      break;
      }
    } //while (!m_Stop)
}

} //namespace itk

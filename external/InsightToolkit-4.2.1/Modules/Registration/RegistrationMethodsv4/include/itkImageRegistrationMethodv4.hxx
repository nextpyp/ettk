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
#ifndef __itkImageRegistrationMethodv4_hxx
#define __itkImageRegistrationMethodv4_hxx

#include "itkImageRegistrationMethodv4.h"

#include "itkDiscreteGaussianImageFilter.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkIterationReporter.h"
#include "itkJointHistogramMutualInformationImageToImageMetricv4.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkShrinkImageFilter.h"

namespace itk
{
/**
 * Constructor
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::ImageRegistrationMethodv4()
{
  this->SetNumberOfRequiredOutputs( 2 );

  this->m_CurrentLevel = 0;
  this->m_CurrentIteration = 0;
  this->m_CurrentMetricValue = 0.0;
  this->m_CurrentConvergenceValue = 0.0;
  this->m_IsConverged = false;

  this->m_CompositeTransform = CompositeTransformType::New();

  typedef IdentityTransform<RealType, ImageDimension> IdentityTransformType;

  typename IdentityTransformType::Pointer defaultFixedInitialTransform = IdentityTransformType::New();
  this->m_FixedInitialTransform = defaultFixedInitialTransform;

  typename IdentityTransformType::Pointer defaultMovingInitialTransform = IdentityTransformType::New();
  this->m_MovingInitialTransform = defaultMovingInitialTransform;

  typedef LinearInterpolateImageFunction<FixedImageType, RealType> DefaultFixedInterpolatorType;
  typename DefaultFixedInterpolatorType::Pointer fixedInterpolator = DefaultFixedInterpolatorType::New();
  this->m_FixedInterpolator = fixedInterpolator;

  typedef LinearInterpolateImageFunction<MovingImageType, RealType> DefaultMovingInterpolatorType;
  typename DefaultMovingInterpolatorType::Pointer movingInterpolator = DefaultMovingInterpolatorType::New();
  this->m_MovingInterpolator = movingInterpolator;

  typedef JointHistogramMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType> MetricForStageOneType;
  typename MetricForStageOneType::Pointer mutualInformationMetric = MetricForStageOneType::New();
  mutualInformationMetric->SetNumberOfHistogramBins( 20 );
  mutualInformationMetric->SetUseMovingImageGradientFilter( false );
  mutualInformationMetric->SetUseFixedImageGradientFilter( false );
  mutualInformationMetric->SetUseFixedSampledPointSet( false );
  this->m_Metric = mutualInformationMetric;

  typedef RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorForAffineTransformType;
  typename ScalesEstimatorForAffineTransformType::Pointer scalesEstimator = ScalesEstimatorForAffineTransformType::New();
  scalesEstimator->SetMetric( mutualInformationMetric );
  scalesEstimator->SetTransformForward( true );

  typedef GradientDescentOptimizerv4 DefaultOptimizerType;
  typename DefaultOptimizerType::Pointer optimizer = DefaultOptimizerType::New();
  optimizer->SetLearningRate( 1.0 );
  optimizer->SetNumberOfIterations( 1000 );
  optimizer->SetScalesEstimator( scalesEstimator );
  this->m_Optimizer = optimizer;

  // By default we set up an affine transform for a 3-level image registration.

  m_NumberOfLevels = 0;
  this->SetNumberOfLevels( 3 );

  this->m_OutputTransform = OutputTransformType::New();

  DecoratedOutputTransformPointer transformDecorator = DecoratedOutputTransformType::New().GetPointer();
  transformDecorator->Set( this->m_OutputTransform );
  this->ProcessObject::SetNthOutput( 0, transformDecorator );

  this->m_ShrinkFactorsPerLevel.SetSize( this->m_NumberOfLevels );
  this->m_ShrinkFactorsPerLevel[0] = 2;
  this->m_ShrinkFactorsPerLevel[1] = 1;
  this->m_ShrinkFactorsPerLevel[2] = 1;

  this->m_SmoothingSigmasPerLevel.SetSize( this->m_NumberOfLevels );
  this->m_SmoothingSigmasPerLevel[0] = 2;
  this->m_SmoothingSigmasPerLevel[1] = 1;
  this->m_SmoothingSigmasPerLevel[2] = 0;

  this->m_MetricSamplingStrategy = NONE;
  this->m_MetricSamplingPercentagePerLevel.SetSize( this->m_NumberOfLevels );
  this->m_MetricSamplingPercentagePerLevel.Fill( 1.0 );
}

template<typename TFixedImage, typename TMovingImage, typename TTransform>
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::~ImageRegistrationMethodv4()
{
}

/*
 * Initialize by setting the interconnects between components.
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
void
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::InitializeRegistrationAtEachLevel( const SizeValueType level )
{
  // Sanity checks

  if ( !this->m_Optimizer )
    {
    itkExceptionMacro( "The optimizer is not present." );
    }
  if ( !this->m_FixedInterpolator )
    {
    itkExceptionMacro( "The fixed image interpolator is not present." );
    }
  if ( !this->m_MovingInterpolator )
    {
    itkExceptionMacro( "The moving image interpolator is not present." );
    }
  if ( !this->m_Metric )
    {
    itkExceptionMacro( "The image metric is not present." );
    }

  this->m_CurrentIteration = 0;
  this->m_CurrentMetricValue = 0.0;
  this->m_CurrentConvergenceValue = 0.0;
  this->m_IsConverged = false;

  this->InvokeEvent( InitializeEvent() );

  // For each level, we adapt the current transform.  For many transforms, e.g.
  // affine, the base transform adaptor does not do anything.  However, in the
  // case of other transforms, e.g. the b-spline and displacement field transforms
  // the fixed parameters are changed to reflect an increase in transform resolution.
  // This could involve increasing the mesh size of the B-spline transform or
  // increase the resolution of the displacement field.

  if( this->m_TransformParametersAdaptorsPerLevel[level] )
    {
    this->m_TransformParametersAdaptorsPerLevel[level]->SetTransform( this->m_OutputTransform );
    this->m_TransformParametersAdaptorsPerLevel[level]->AdaptTransformParameters();
    }

  // At each resolution, we can
  //   1. subsample the reference domain (typically the fixed image) and/or
  //   2. smooth the fixed and moving images.

  typedef ShrinkImageFilter<FixedImageType, FixedImageType> ShrinkFilterType;
  typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
  shrinkFilter->SetShrinkFactors( this->m_ShrinkFactorsPerLevel[level] );
  shrinkFilter->SetInput( this->GetFixedImage() );
  shrinkFilter->Update();

  typedef DiscreteGaussianImageFilter<FixedImageType, FixedImageType> FixedImageSmoothingFilterType;
  typename FixedImageSmoothingFilterType::Pointer fixedImageSmoothingFilter = FixedImageSmoothingFilterType::New();
  fixedImageSmoothingFilter->SetUseImageSpacingOn();
  fixedImageSmoothingFilter->SetVariance( vnl_math_sqr( this->m_SmoothingSigmasPerLevel[level] ) );
  fixedImageSmoothingFilter->SetMaximumError( 0.01 );
  fixedImageSmoothingFilter->SetInput( this->GetFixedImage() );

  this->m_FixedSmoothImage = fixedImageSmoothingFilter->GetOutput();
  this->m_FixedSmoothImage->Update();
  this->m_FixedSmoothImage->DisconnectPipeline();

  typedef DiscreteGaussianImageFilter<MovingImageType, MovingImageType> MovingImageSmoothingFilterType;
  typename MovingImageSmoothingFilterType::Pointer movingImageSmoothingFilter = MovingImageSmoothingFilterType::New();
  movingImageSmoothingFilter->SetUseImageSpacingOn();
  movingImageSmoothingFilter->SetVariance( vnl_math_sqr( this->m_SmoothingSigmasPerLevel[level] ) );
  movingImageSmoothingFilter->SetMaximumError( 0.01 );
  movingImageSmoothingFilter->SetInput( this->GetMovingImage() );

  this->m_MovingSmoothImage = movingImageSmoothingFilter->GetOutput();
  this->m_MovingSmoothImage->Update();
  this->m_MovingSmoothImage->DisconnectPipeline();

  // Set-up the composite transform at initialization
  if( level == 0 )
    {
    this->m_CompositeTransform->ClearTransformQueue();

    // If the moving initial transform is a composite transform, unroll
    // it into m_CompositeTransform.  This is a temporary fix to accommodate
    // the lack of support for calculating the jacobian in the composite
    // transform.

    this->m_CompositeTransform->AddTransform( this->m_MovingInitialTransform );
    this->m_CompositeTransform->AddTransform( this->m_OutputTransform );
    this->m_CompositeTransform->FlattenTransformQueue();
    }
  this->m_CompositeTransform->SetOnlyMostRecentTransformToOptimizeOn();


  // Update the image metric

  this->m_Metric->SetFixedTransform( this->m_FixedInitialTransform );
  this->m_Metric->SetMovingTransform( this->m_CompositeTransform );
  this->m_Metric->SetFixedInterpolator( this->m_FixedInterpolator );
  this->m_Metric->SetMovingInterpolator( this->m_MovingInterpolator );
  this->m_Metric->SetFixedImage( this->m_FixedSmoothImage );
  this->m_Metric->SetMovingImage( this->m_MovingSmoothImage );
  this->m_Metric->SetVirtualDomainFromImage( shrinkFilter->GetOutput() );

  if( this->m_MetricSamplingStrategy != NONE )
    {
    this->SetMetricSamplePoints();
    }

  // Update the optimizer

  this->m_Optimizer->SetMetric( this->m_Metric );

  if( ( this->m_Optimizer->GetScales() ).Size() != this->m_OutputTransform->GetNumberOfLocalParameters() )
    {
    typedef typename OptimizerType::ScalesType ScalesType;
    ScalesType scales;
    scales.SetSize( this->m_OutputTransform->GetNumberOfLocalParameters() );
    scales.Fill( NumericTraits<typename ScalesType::ValueType>::OneValue() );
    this->m_Optimizer->SetScales( scales );
    }
}

/*
 * Start the registration
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
void
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::GenerateData()
{
  for( this->m_CurrentLevel = 0; this->m_CurrentLevel < this->m_NumberOfLevels; this->m_CurrentLevel++ )
    {
    this->InitializeRegistrationAtEachLevel( this->m_CurrentLevel );

    this->m_Metric->Initialize();
    this->m_Optimizer->StartOptimization();
    }
}

/**
 * Set the moving transform adaptors per stage
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
void
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::SetTransformParametersAdaptorsPerLevel( TransformParametersAdaptorsContainerType & adaptors )
{
  if( this->m_NumberOfLevels != adaptors.size() )
    {
    itkExceptionMacro( "The number of levels does not equal the number array size." );
    }
  else
    {
    itkDebugMacro( "Setting the transform parameter adaptors." );
    this->m_TransformParametersAdaptorsPerLevel = adaptors;
    this->Modified();
    }
}

/**
 * Get the moving transform adaptors per stage
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
const typename  ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>::TransformParametersAdaptorsContainerType &
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::GetTransformParametersAdaptorsPerLevel() const
{
  return this->m_TransformParametersAdaptorsPerLevel;
}

/**
 * Set the number of levels
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
void
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::SetNumberOfLevels( const SizeValueType numberOfLevels )
{
  if( this->m_NumberOfLevels != numberOfLevels )
    {
    this->m_NumberOfLevels = numberOfLevels;

    // Set default transform adaptors which don't do anything to the input transform
    // Similarly, fill in some default values for the shrink factors, smoothing sigmas,
    // and learning rates.

    this->m_TransformParametersAdaptorsPerLevel.clear();
    for( SizeValueType level = 0; level < this->m_NumberOfLevels; level++ )
      {
      typename TransformParametersAdaptorType::Pointer transformParametersAdaptor = TransformParametersAdaptorType::New();
      this->m_TransformParametersAdaptorsPerLevel.push_back( transformParametersAdaptor.GetPointer() );
      }

    this->m_ShrinkFactorsPerLevel.SetSize( this->m_NumberOfLevels );
    this->m_ShrinkFactorsPerLevel.Fill( 1 );

    this->m_SmoothingSigmasPerLevel.SetSize( this->m_NumberOfLevels );
    this->m_SmoothingSigmasPerLevel.Fill( 1.0 );

    this->m_MetricSamplingPercentagePerLevel.SetSize( this->m_NumberOfLevels );
    this->m_MetricSamplingPercentagePerLevel.Fill( 1.0 );

    this->Modified();
    }
}

/**
 * Get the metric samples
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
void
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::SetMetricSamplePoints()
{
  typename MetricSamplePointSetType::Pointer samplePointSet = MetricSamplePointSetType::New();
  samplePointSet->Initialize();

  typedef typename MetricSamplePointSetType::PointType SamplePointType;

  typedef typename MetricType::VirtualImageType         VirtualDomainImageType;
  typedef typename VirtualDomainImageType::RegionType   VirtualDomainRegionType;

  const VirtualDomainImageType * virtualImage = this->m_Metric->GetVirtualImage();
  const VirtualDomainRegionType & virtualDomainRegion = virtualImage->GetRequestedRegion();
  const typename VirtualDomainImageType::SpacingType virtualSpacing = virtualImage->GetSpacing();

  unsigned long sampleCount = virtualDomainRegion.GetNumberOfPixels();

  typedef typename Statistics::MersenneTwisterRandomVariateGenerator RandomizerType;
  typename RandomizerType::Pointer randomizer = RandomizerType::New();
  randomizer->SetSeed( 1234 );

  unsigned long index = 0;

  switch( this->m_MetricSamplingStrategy )
    {
    case REGULAR:
      {
      sampleCount = static_cast<unsigned long>( vcl_ceil( 1.0 / this->m_MetricSamplingPercentagePerLevel[this->m_CurrentLevel] ) );

      unsigned long count = 0;
      ImageRegionConstIteratorWithIndex<VirtualDomainImageType> It( virtualImage, virtualDomainRegion );
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        if( count % sampleCount == 0 )
          {
          SamplePointType point;
          virtualImage->TransformIndexToPhysicalPoint( It.GetIndex(), point );

          // randomly perturb the point within a voxel (approximately)
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            point[d] += randomizer->GetNormalVariate() / 3.0 * virtualSpacing[d];
            }
          samplePointSet->SetPoint( index, point );
          index++;
          }
        count++;
        }
      break;
      }
    case RANDOM:
      {
      sampleCount = static_cast<unsigned long>( static_cast<float>( sampleCount ) * this->m_MetricSamplingPercentagePerLevel[this->m_CurrentLevel] );

      ImageRandomConstIteratorWithIndex<VirtualDomainImageType> ItR( virtualImage, virtualDomainRegion );
      ItR.SetNumberOfSamples( sampleCount );
      for( ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR )
        {
        SamplePointType point;
        virtualImage->TransformIndexToPhysicalPoint( ItR.GetIndex(), point );

        // randomly perturb the point within a voxel (approximately)
        for ( unsigned int d = 0; d < ImageDimension; d++ )
          {
          point[d] += randomizer->GetNormalVariate() / 3.0 * virtualSpacing[d];
          }
        samplePointSet->SetPoint( index, point );
        index++;
        }
      break;
      }
    default:
      {
      itkExceptionMacro( "Invalid sampling strategy requested." );
      }
    }

  this->m_Metric->SetFixedSampledPointSet( samplePointSet );
  this->m_Metric->SetUseFixedSampledPointSet( true );
}

/*
 * PrintSelf
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
void
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Moving interpolator:" << std::endl;
  this->m_MovingInterpolator->Print( std::cout, indent );
  os << indent << "Fixed interpolator:" << std::endl;
  this->m_FixedInterpolator->Print( std::cout, indent );

  os << "Number of levels = " << this->m_NumberOfLevels << std::endl;

  os << indent << "Shrink factors: " << this->m_ShrinkFactorsPerLevel << std::endl;
  os << indent << "Smoothing sigmas: " << this->m_SmoothingSigmasPerLevel << std::endl;

  os << indent << "Metric sampling strategy: " << this->m_MetricSamplingStrategy << std::endl;

  os << indent << "Metric sampling percentage: ";
  for( unsigned int i = 0; i < this->m_NumberOfLevels; i++ )
    {
    os << this->m_MetricSamplingPercentagePerLevel[i] << " ";
    }
  os << std::endl;
}

/*
 *  Get output transform
 */
template<typename TFixedImage, typename TMovingImage, typename TTransform>
const typename ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>::DecoratedOutputTransformType *
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::GetOutput() const
{
  return static_cast<const DecoratedOutputTransformType *>( this->ProcessObject::GetOutput( 0 ) );
}

template<typename TFixedImage, typename TMovingImage, typename TTransform>
DataObject::Pointer
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::MakeOutput( DataObjectPointerArraySizeType output )
{
  switch ( output )
    {
    case 0:
      return static_cast<DataObject *>( DecoratedOutputTransformType::New().GetPointer() );
      break;
    default:
      itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
      return 0;
    }
}

template<typename TFixedImage, typename TMovingImage, typename TTransform>
void
ImageRegistrationMethodv4<TFixedImage, TMovingImage, TTransform>
::SetMetricSamplingPercentage( const RealType samplingPercentage )
{
  MetricSamplingPercentageArrayType samplingPercentagePerLevel;
  samplingPercentagePerLevel.SetSize( this->m_NumberOfLevels );
  samplingPercentagePerLevel.Fill( samplingPercentage );
  this->SetMetricSamplingPercentagePerLevel( samplingPercentagePerLevel );
}

} // end namespace itk

#endif

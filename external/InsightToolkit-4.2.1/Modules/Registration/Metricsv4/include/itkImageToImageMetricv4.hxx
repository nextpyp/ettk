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
#ifndef __itkImageToImageMetricv4_hxx
#define __itkImageToImageMetricv4_hxx

#include "itkImageToImageMetricv4.h"
#include "itkPixelTraits.h"
#include "itkDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkCentralDifferenceImageFunction.h"

namespace itk
{

template<class TFixedImage,class TMovingImage,class TVirtualImage>
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::ImageToImageMetricv4()
{
  /* Interpolators. Default to linear. */
  typedef LinearInterpolateImageFunction< FixedImageType,
                                          CoordinateRepresentationType >
                                                  FixedLinearInterpolatorType;
  typedef LinearInterpolateImageFunction< MovingImageType,
                                          CoordinateRepresentationType >
                                                  MovingLinearInterpolatorType;
  this->m_FixedInterpolator  = FixedLinearInterpolatorType::New();
  this->m_MovingInterpolator = MovingLinearInterpolatorType::New();

  /* Setup default gradient filter. It gets initialized with default
   * parameters during Initialize. */
  this->m_DefaultFixedImageGradientFilter  = DefaultFixedImageGradientFilter::New();
  this->m_DefaultMovingImageGradientFilter = DefaultMovingImageGradientFilter::New();
  this->m_FixedImageGradientFilter         = this->m_DefaultFixedImageGradientFilter;
  this->m_MovingImageGradientFilter        = this->m_DefaultMovingImageGradientFilter;

  /* Interpolators for image gradient filters */
  this->m_FixedImageGradientInterpolator  = FixedImageGradientInterpolatorType::New();
  this->m_MovingImageGradientInterpolator = MovingImageGradientInterpolatorType::New();

  /* Setup default gradient image function */
  typedef CentralDifferenceImageFunction<FixedImageType,
                                         CoordinateRepresentationType>
                                          FixedCentralDifferenceCalculatorType;
  typedef CentralDifferenceImageFunction<MovingImageType,
                                         CoordinateRepresentationType>
                                          MovingCentralDifferenceCalculatorType;
  typename FixedCentralDifferenceCalculatorType::Pointer
                  fixedCalculator       = FixedCentralDifferenceCalculatorType::New();
  fixedCalculator->UseImageDirectionOn();
  this->m_FixedImageGradientCalculator  = fixedCalculator;
  typename MovingCentralDifferenceCalculatorType::Pointer
                  movingCalculator      = MovingCentralDifferenceCalculatorType::New();
  movingCalculator->UseImageDirectionOn();
  this->m_MovingImageGradientCalculator = movingCalculator;

  /* Setup default options assuming dense-sampling */
  this->m_UseFixedImageGradientFilter  = true;
  this->m_UseMovingImageGradientFilter = true;
  this->m_UseFixedSampledPointSet      = false;

  this->m_FloatingPointCorrectionResolution = 1e6;
  this->m_UseFloatingPointCorrection = false;

  this->m_HaveMadeGetValueWarning = false;
  this->m_NumberOfSkippedFixedSampledPoints = 0;

  this->m_Value = NumericTraits<MeasureType>::max();
  this->m_DerivativeResult = NULL;
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::~ImageToImageMetricv4()
{
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::Initialize() throw ( itk::ExceptionObject )
{
  itkDebugMacro("Initialize entered");

  /* Verify things are connected */
  if ( this->m_FixedImage.IsNull() )
    {
    itkExceptionMacro(<< "FixedImage is not present");
    }
  if ( this->m_MovingImage.IsNull() )
    {
    itkExceptionMacro(<< "MovingImage is not present");
    }
  if ( this->m_FixedTransform.IsNull() )
    {
    itkExceptionMacro(<< "FixedTransform is not present");
    }
  if ( this->m_MovingTransform.IsNull() )
    {
    itkExceptionMacro(<< "MovingTransform is not present");
    }

  // If the image is provided by a source, update the source.
  if ( this->m_MovingImage->GetSource() )
    {
    this->m_MovingImage->GetSource()->Update();
    }

  // If the image is provided by a source, update the source.
  if ( this->m_FixedImage->GetSource() )
    {
    this->m_FixedImage->GetSource()->Update();
    }

  /* If a virtual image has not been set or created,
   * create one from fixed image settings */
  if( ! this->m_UserHasSetVirtualDomain )
    {
    /* Instantiate a virtual image, but do not call Allocate to allocate
     * the data, to save memory. We don't need data. We'll simply be iterating
     * over the image to get indecies and transform to points.
     * Note that it will be safer to have a dedicated VirtualImage class
     * that prevents accidental access of data. */
    /* Just copy information from fixed image */
    VirtualImagePointer image = VirtualImageType::New();
    image->CopyInformation( this->m_FixedImage );
    /* CopyInformation does not copy buffered region */
    image->SetBufferedRegion( this->m_FixedImage->GetBufferedRegion() );
    image->SetRequestedRegion( this->m_FixedImage->GetRequestedRegion() );
    this->SetVirtualDomainFromImage( image );
    }

  /*
   * Superclass Initialize.
   * Requires the above actions to already have been taken.
   */
  Superclass::Initialize();

  /* Map the fixed samples into the virtual domain and store in
   * a searpate point set. */
  if( this->m_UseFixedSampledPointSet )
    {
    this->MapFixedSampledPointSetToVirtual();
    }

  /* Inititialize interpolators. */
  itkDebugMacro("Initialize Interpolators");
  this->m_FixedInterpolator->SetInputImage( this->m_FixedImage );
  this->m_MovingInterpolator->SetInputImage( this->m_MovingImage );

  /* Setup for image gradient calculations. */
  if( ! this->m_UseFixedImageGradientFilter )
    {
    itkDebugMacro("Initialize FixedImageGradientCalculator");
    this->m_FixedImageGradientImage = NULL;
    this->m_FixedImageGradientCalculator->SetInputImage(this->m_FixedImage);
    }
  if( ! this->m_UseMovingImageGradientFilter )
    {
    itkDebugMacro("Initialize MovingImageGradientCalculator");
    this->m_MovingImageGradientImage = NULL;
    this->m_MovingImageGradientCalculator->SetInputImage(this->m_MovingImage);
    }

  /* Initialize default gradient image filters. */
  itkDebugMacro("InitializeDefaultFixedImageGradientFilter");
  this->InitializeDefaultFixedImageGradientFilter();
  itkDebugMacro("InitializeDefaultMovingImageGradientFilter");
  this->InitializeDefaultMovingImageGradientFilter();

  /* If user set to use a pre-calculated fixed gradient image,
   * then we need to calculate the gradient image.
   * We only need to compute once since the fixed transform isn't
   * optimized. */
  if ( this->m_UseFixedImageGradientFilter )
    {
    itkDebugMacro("Initialize: ComputeFixedImageGradientFilterImage");
    this->ComputeFixedImageGradientFilterImage();
    }

  /* Compute gradient image for moving image. Needed now for
   * derived classes that use it before InitializeForIteration is called.
   * It's also computed at begin of every iteration. */
  if( this->m_UseMovingImageGradientFilter )
    {
    itkDebugMacro("Initialize: ComputeMovingImageGradientFilterImage");
    this->ComputeMovingImageGradientFilterImage();
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
typename ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage>::MeasureType
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage>
::GetValue() const
{
  /* As long as this is done as a simple inefficient implementation, give
   * the user a warning. The intention is to provide a more efficient
   * implementation here eventually, that avoids derivative calcualtions.
   * Derived classes may override this method and provide their own
   * efficient implementations as appropriate. */
  if( ! this->m_HaveMadeGetValueWarning )
    {
    itkWarningMacro("Using ImageToImageMetricv4::GetValue which is a "
                    "temporary, inefficient implementation. " );
    this->m_HaveMadeGetValueWarning = true;
    }
  DerivativeType derivative;
  MeasureType value;
  this->GetValueAndDerivative( value, derivative );
  return value;
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage>
::GetDerivative( DerivativeType & derivative ) const
{
  MeasureType value;
  this->GetValueAndDerivative( value, derivative );
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage>
::GetValueAndDerivative( MeasureType & value, DerivativeType & derivative ) const
{
  this->m_DerivativeResult = &derivative;
  this->InitializeForIteration();

  // Do the threaded processing using the appropriate
  // GetValueAndDerivativeThreader. Results get written to
  // member vars.
  this->GetValueAndDerivativeExecute();

  value = this->m_Value;
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::GetValueAndDerivativeExecute() const
{
  if( this->m_UseFixedSampledPointSet ) // sparse sampling
    {
    SizeValueType numberOfPoints = this->GetNumberOfDomainPoints();
    if( numberOfPoints < 1 )
      {
      itkExceptionMacro("VirtualSampledPointSet must have 1 or more points.");
      }
    typename ImageToImageMetricv4GetValueAndDerivativeThreader< ThreadedIndexedContainerPartitioner, Self >::DomainType range;
    range[0] = 0;
    range[1] = numberOfPoints - 1;
    this->m_SparseGetValueAndDerivativeThreader->Execute( const_cast< Self* >(this), range );
    }
  else // dense sampling
    {
    this->m_DenseGetValueAndDerivativeThreader->Execute( const_cast< Self* >(this), this->GetVirtualRegion() );
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::InitializeForIteration() const
{
  if( this->m_DerivativeResult )
    {
    /* This size always comes from the active transform */
    const NumberOfParametersType globalDerivativeSize = this->GetNumberOfParameters();
    if( this->m_DerivativeResult->GetSize() != globalDerivativeSize )
      {
      this->m_DerivativeResult->SetSize( globalDerivativeSize );
      }
    /* Clear derivative final result. */
    this->m_DerivativeResult->Fill( NumericTraits< DerivativeValueType >::Zero );
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
bool
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateFixedPoint(
                         const VirtualIndexType & itkNotUsed(index),
                         const VirtualPointType & virtualPoint,
                         const bool computeImageGradient,
                         FixedImagePointType & mappedFixedPoint,
                         FixedImagePixelType & mappedFixedPixelValue,
                         FixedImageGradientType & mappedFixedImageGradient ) const
{
  bool pointIsValid = true;
  mappedFixedPixelValue = NumericTraits<FixedImagePixelType>::Zero;

  // map the point into fixed space
  mappedFixedPoint = this->m_FixedTransform->TransformPoint( virtualPoint );

  // check against the mask if one is assigned
  if ( this->m_FixedImageMask )
    {
    // Check if mapped point is within the support region of the fixed image
    // mask
    pointIsValid = this->m_FixedImageMask->IsInside( mappedFixedPoint );
    if( ! pointIsValid )
      {
      return pointIsValid;
      }
    }

  // Check if mapped point is inside image buffer
  pointIsValid = this->m_FixedInterpolator->IsInsideBuffer(mappedFixedPoint);
  if( ! pointIsValid )
    {
    return pointIsValid;
    }

  mappedFixedPixelValue = this->m_FixedInterpolator->Evaluate(mappedFixedPoint);
  if( computeImageGradient )
    {
    this->ComputeFixedImageGradientAtPoint( mappedFixedPoint,
                                     mappedFixedImageGradient );
    }

  return pointIsValid;
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
bool
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::TransformAndEvaluateMovingPoint(
                         const VirtualIndexType & itkNotUsed(index),
                         const VirtualPointType & virtualPoint,
                         const bool computeImageGradient,
                         MovingImagePointType & mappedMovingPoint,
                         MovingImagePixelType & mappedMovingPixelValue,
                         MovingImageGradientType & mappedMovingImageGradient ) const
{
  bool pointIsValid = true;
  mappedMovingPixelValue = NumericTraits<MovingImagePixelType>::Zero;

  // map the point into moving space
  mappedMovingPoint = this->m_MovingTransform->TransformPoint( virtualPoint );

  // check against the mask if one is assigned
  if ( this->m_MovingImageMask )
    {
    // Check if mapped point is within the support region of the fixed image
    // mask
    pointIsValid = this->m_MovingImageMask->IsInside( mappedMovingPoint );
    if( ! pointIsValid )
      {
      return pointIsValid;
      }
    }

  // Check if mapped point is inside image buffer
  pointIsValid = this->m_MovingInterpolator->IsInsideBuffer(mappedMovingPoint);
  if( ! pointIsValid )
    {
    return pointIsValid;
    }

  mappedMovingPixelValue = this->m_MovingInterpolator->Evaluate( mappedMovingPoint );
  if( computeImageGradient )
    {
    this->ComputeMovingImageGradientAtPoint( mappedMovingPoint, mappedMovingImageGradient );
    }

  return pointIsValid;
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageGradientAtPoint( const FixedImagePointType & mappedPoint,
                             FixedImageGradientType & gradient ) const
{
  if ( this->m_UseFixedImageGradientFilter )
    {
    gradient = m_FixedImageGradientInterpolator->Evaluate( mappedPoint );
    }
  else
    {
    // if not using the gradient image
    gradient = this->m_FixedImageGradientCalculator->Evaluate( mappedPoint );
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageGradientAtPoint(
                              const MovingImagePointType & mappedPoint,
                              MovingImageGradientType & gradient ) const
{
  if ( this->m_UseMovingImageGradientFilter )
    {
    gradient = m_MovingImageGradientInterpolator->Evaluate( mappedPoint );
    }
  else
    {
    // if not using the gradient image
    gradient = this->m_MovingImageGradientCalculator->Evaluate(mappedPoint);
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageGradientAtIndex(
                              const VirtualIndexType & index,
                              FixedImageGradientType & gradient ) const
{
  if ( this->m_UseFixedImageGradientFilter )
    {
    gradient = this->m_FixedImageGradientImage->GetPixel(index);
    }
  else
    {
    // if not using the gradient image
    gradient = this->m_FixedImageGradientCalculator->EvaluateAtIndex(index);
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageGradientAtIndex(
                              const VirtualIndexType & index,
                              MovingImageGradientType & gradient ) const
{
  if ( this->m_UseMovingImageGradientFilter )
    {
    gradient = this->m_MovingImageGradientImage->GetPixel(index);
    }
  else
    {
    // if not using the gradient image
    gradient = this->m_MovingImageGradientCalculator->EvaluateAtIndex(index);
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::ComputeFixedImageGradientFilterImage()
{
  this->m_FixedImageGradientFilter->SetInput( this->m_FixedImage );
  this->m_FixedImageGradientFilter->Update();
  this->m_FixedImageGradientImage = this->m_FixedImageGradientFilter->GetOutput();
  this->m_FixedImageGradientInterpolator->SetInputImage( this->m_FixedImageGradientImage );
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::ComputeMovingImageGradientFilterImage() const
{
  this->m_MovingImageGradientFilter->SetInput( this->m_MovingImage );
  this->m_MovingImageGradientFilter->Update();
  this->m_MovingImageGradientImage = this->m_MovingImageGradientFilter->GetOutput();
  this->m_MovingImageGradientInterpolator->SetInputImage( this->m_MovingImageGradientImage );
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::InitializeDefaultFixedImageGradientFilter()
{
  const typename FixedImageType::SpacingType & spacing = this->m_FixedImage->GetSpacing();
  double maximumSpacing = 0.0;
  for ( ImageDimensionType i = 0; i < FixedImageDimension; i++ )
    {
    if ( spacing[i] > maximumSpacing )
      {
      maximumSpacing = spacing[i];
      }
    }
  this->m_DefaultFixedImageGradientFilter->SetSigma( maximumSpacing );
  this->m_DefaultFixedImageGradientFilter->SetNormalizeAcrossScale( true );
  this->m_DefaultFixedImageGradientFilter->SetNumberOfThreads( this->GetMaximumNumberOfThreads() );
  this->m_DefaultFixedImageGradientFilter->SetUseImageDirection( true );
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::InitializeDefaultMovingImageGradientFilter()
{
  const typename MovingImageType::SpacingType & spacing = this->m_MovingImage->GetSpacing();
  double maximumSpacing = 0.0;
  for ( ImageDimensionType i = 0; i < MovingImageDimension; i++ )
    {
    if ( spacing[i] > maximumSpacing )
      {
      maximumSpacing = spacing[i];
      }
    }
  this->m_DefaultMovingImageGradientFilter->SetSigma(maximumSpacing);
  this->m_DefaultMovingImageGradientFilter->SetNormalizeAcrossScale(true);
  this->m_DefaultMovingImageGradientFilter->SetNumberOfThreads(this->GetMaximumNumberOfThreads());
  this->m_DefaultMovingImageGradientFilter->SetUseImageDirection(true);
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::SetMaximumNumberOfThreads( const ThreadIdType number )
{
  if( number != this->m_SparseGetValueAndDerivativeThreader->GetMaximumNumberOfThreads() )
    {
    this->m_SparseGetValueAndDerivativeThreader->SetMaximumNumberOfThreads( number );
    this->Modified();
    }
  if( number != this->m_DenseGetValueAndDerivativeThreader->GetMaximumNumberOfThreads() )
    {
    this->m_DenseGetValueAndDerivativeThreader->SetMaximumNumberOfThreads( number );
    this->Modified();
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
ThreadIdType
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::GetMaximumNumberOfThreads() const
{
  if( this->m_UseFixedSampledPointSet )
    {
    return this->m_SparseGetValueAndDerivativeThreader->GetMaximumNumberOfThreads();
    }
  return  this->m_DenseGetValueAndDerivativeThreader->GetMaximumNumberOfThreads();
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
ThreadIdType
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::GetNumberOfThreadsUsed() const
{
  if( this->m_UseFixedSampledPointSet )
    {
    return this->m_SparseGetValueAndDerivativeThreader->GetNumberOfThreadsUsed();
    }
  else
    {
    return this->m_DenseGetValueAndDerivativeThreader->GetNumberOfThreadsUsed();
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage>
::MapFixedSampledPointSetToVirtual()
{
  this->m_VirtualSampledPointSet = VirtualPointSetType::New();
  this->m_VirtualSampledPointSet->Initialize();

  typedef typename FixedSampledPointSetType::PointsContainer PointsContainer;
  typename PointsContainer::ConstPointer
    points = this->m_FixedSampledPointSet->GetPoints();
  typename PointsContainer::ConstIterator fixedIt = points->Begin();

  typename FixedTransformType::InverseTransformBasePointer
    inverseTransform = this->m_FixedTransform->GetInverseTransform();
  if( inverseTransform.IsNull() )
    {
    itkExceptionMacro("Unable to get inverse transform for mapping sampled "
                      " point set.");
    }

  this->m_NumberOfSkippedFixedSampledPoints = 0;
  SizeValueType virtualIndex = 0;
  while( fixedIt != points->End() )
    {
    typename FixedSampledPointSetType::PointType point = inverseTransform->TransformPoint( fixedIt.Value() );
    typename VirtualImageType::IndexType tempIndex;
    /* Verify that the point is valid. We may be working with a resized virtual domain,
     * and a fixed sampled point list that was created before the resizing. */
    if( this->TransformPhysicalPointToVirtualIndex( point, tempIndex ) )
      {
      this->m_VirtualSampledPointSet->SetPoint( virtualIndex, point );
      virtualIndex++;
      }
    else
      {
      this->m_NumberOfSkippedFixedSampledPoints++;
      }
    ++fixedIt;
    }
  if( this->m_VirtualSampledPointSet->GetNumberOfPoints() == 0 )
    {
    itkExceptionMacro("The virtual sampled point set has zero points because "
                      "no fixed sampled points were within the virtual "
                      "domain after mapping. There are no points to evaulate.");
    }
}


template<class TFixedImage,class TMovingImage,class TVirtualImage>
SizeValueType
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage>
::GetNumberOfDomainPoints() const
{
  if( this->m_UseFixedSampledPointSet )
    {
    //The virtual sampled point set holds the actual points
    // over which we're evaluating over.
    return this->m_VirtualSampledPointSet->GetNumberOfPoints();
    }
  else
    {
    typename VirtualImageType::RegionType region = this->GetVirtualRegion();
    return region.GetNumberOfPixels();
    }
}

template<class TFixedImage,class TMovingImage,class TVirtualImage>
void
ImageToImageMetricv4<TFixedImage, TMovingImage, TVirtualImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "ImageToImageMetricv4: " << std::endl
     << indent << "GetUseFixedImageGradientFilter: " << this->GetUseFixedImageGradientFilter() << std::endl
     << indent << "GetUseMovingImageGradientFilter: " << this->GetUseMovingImageGradientFilter() << std::endl
     << indent << "UseFloatingPointCorrection: " << this->GetUseFloatingPointCorrection() << std::endl
     << indent << "FloatingPointCorrectionResolution: " << this->GetFloatingPointCorrectionResolution() << std::endl;

  if( this->GetFixedImage() != NULL )
    {
    os << indent << "FixedImage: " << this->GetFixedImage() << std::endl;
    }
  else
    {
    os << indent << "FixedImage is NULL." << std::endl;
    }
  if( this->GetMovingImage() != NULL )
    {
    os << indent << "MovingImage: " << this->GetMovingImage() << std::endl;
    }
  else
    {
    os << indent << "MovingImage is NULL." << std::endl;
    }
  if( this->GetFixedTransform() != NULL )
    {
    os << indent << "FixedTransform: " << this->GetFixedTransform() << std::endl;
    }
  else
    {
    os << indent << "FixedTransform is NULL." << std::endl;
    }
  if( this->GetMovingTransform() != NULL )
    {
    os << indent << "MovingTransform: " << this->GetMovingTransform() << std::endl;
    }
  else
    {
    os << indent << "MovingTransform is NULL." << std::endl;
    }
  if( this->GetFixedImageMask() != NULL )
    {
    os << indent << "FixedImageMask: " << this->GetFixedImageMask() << std::endl;
    }
  else
    {
    os << indent << "FixedImageMask is NULL." << std::endl;
    }
  if( this->GetMovingImageMask() != NULL )
    {
    os << indent << "MovingImageMask: " << this->GetMovingImageMask() << std::endl;
    }
  else
    {
    os << indent << "MovingImageMask is NULL." << std::endl;
    }
}

}//namespace itk

#endif

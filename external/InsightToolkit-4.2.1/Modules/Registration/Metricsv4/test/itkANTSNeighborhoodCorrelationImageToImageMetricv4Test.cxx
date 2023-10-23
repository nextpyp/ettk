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

#include "itkTranslationTransform.h"

#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkTestingMacros.h"

/**
 * Test program for ANTSNeighborhoodCorrelationImageToImageMetricv4,
 * using a synthectic image and initial displacement.
 *
 */

template<class ImagePointerType, class DerivativeType>
void ANTSNeighborhoodCorrelationImageToImageMetricv4Test_PrintDerivativeAsVectorImage(ImagePointerType image, DerivativeType &derivative, itk::SizeValueType vecdim){

    typedef typename ImagePointerType::ObjectType ImageType;
    typename ImageType::RegionType imageRegion = image->GetBufferedRegion();

    // only display the first slice
    itk::SizeValueType dim0 = imageRegion.GetSize()[0];
    itk::SizeValueType dim1 = imageRegion.GetSize()[1];

    typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType it(image, imageRegion);
    it.GoToBegin();
    itk::SizeValueType cnt = 0;
    for (itk::SizeValueType ycnt = 0; ycnt < dim1; ycnt++)
      {
      for (itk::SizeValueType xcnt = 0; xcnt < dim0; xcnt++)
        {
        std::cout << '[';
        for(itk::SizeValueType d = 0; d < vecdim-1; d++)
          {
          std::cout << derivative[cnt*vecdim + d] <<",";
          }
          std::cout << derivative[cnt * vecdim + vecdim-1] << ']' << "\t";
          ++it;
          ++cnt;
        }
        std::cout << std::endl;
      }
    return;
}


template<class ImageType>
void ANTSNeighborhoodCorrelationImageToImageMetricv4Test_PrintImage(ImageType *imageP) {

    typedef typename ImageType::ConstPointer ImageConstPointerType;
    ImageConstPointerType image = imageP;

    typename ImageType::RegionType imageRegion = image->GetBufferedRegion();

    // only display the first slice
    itk::SizeValueType dim0 = imageRegion.GetSize()[0];
    itk::SizeValueType dim1 = imageRegion.GetSize()[1];

    typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType it(image, imageRegion);
    it.GoToBegin();

    for (itk::SizeValueType ycnt = 0; ycnt < dim1; ycnt++)
      {
      for (itk::SizeValueType xcnt = 0; xcnt < dim0; xcnt++)
        {
        std::cout << it.Get() << "\t";
        ++it;
        }
        std::cout << std::endl;
      }
    return;
}

template<class ImagePointerType>
void ANTSNeighborhoodCorrelationImageToImageMetricv4Test_PrintImage(const ImagePointerType &image) {

    typedef typename ImagePointerType::ObjectType ImageType;
    typename ImageType::RegionType imageRegion = image->GetBufferedRegion();

    // only display the first slice
    itk::SizeValueType dim0 = imageRegion.GetSize()[0];
    itk::SizeValueType dim1 = imageRegion.GetSize()[1];

    typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType it(image, imageRegion);
    it.GoToBegin();

    for (itk::SizeValueType ycnt = 0; ycnt < dim1; ycnt++)
      {
      for (itk::SizeValueType xcnt = 0; xcnt < dim0; xcnt++)
        {
        std::cout << it.Get() << "\t";
        ++it;
        }
        std::cout << std::endl;
      }
    return;
}

int itkANTSNeighborhoodCorrelationImageToImageMetricv4Test( int, char ** const )
{

//    MultiThreader::SetGlobalMaximumNumberOfThreads(1);

    const itk::SizeValueType ImageDimension = 2;

    typedef itk::Image<double, ImageDimension>  ImageType;
    typedef ImageType::Pointer                  ImagePointerType;
    typedef ImageType::RegionType               RegionType;

    typedef itk::Vector<double, ImageDimension>     VectorType;
    typedef itk::Image<VectorType, ImageDimension>  VectorImageType;

    typedef itk::Transform<double, ImageDimension>  TransformType;
    typedef itk::IdentityTransform<double, ImageDimension>
                                                    IdentityTransformType;
    typedef itk::CompositeTransform<double, ImageDimension>
                                                    CompositeTransformType;
    typedef itk::TranslationTransform<double, ImageDimension>
                                                    TranslationTransformType;
    typedef itk::DisplacementFieldTransform<double, ImageDimension>
                                                    DisplacementTransformType;
    typedef DisplacementTransformType::DisplacementFieldType FieldType;

    IdentityTransformType::Pointer transformFId = IdentityTransformType::New();

    IdentityTransformType::Pointer transformMId = IdentityTransformType::New();
    if(transformMId.IsNull())
      {
      std::cerr << "transformMId == NULL" << std::endl;
      return EXIT_FAILURE;
      }
    DisplacementTransformType::Pointer transformMdisplacement = DisplacementTransformType::New();
    TranslationTransformType::Pointer transformMtranslation =
            TranslationTransformType::New();
    TranslationTransformType::Pointer transformMtranslation2 =
            TranslationTransformType::New();
    CompositeTransformType::Pointer transformMComp =
            CompositeTransformType::New();
    CompositeTransformType::Pointer transformFComp =
            CompositeTransformType::New();


    const itk::SizeValueType imageSize = 6;

    ImageType::SizeType size;
    size.Fill(imageSize);
    ImageType::IndexType index;
    index.Fill(0);
    ImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(index);
    ImageType::SpacingType spacing;
    spacing.Fill(1.0);
    ImageType::PointType origin;
    origin.Fill(0);
    ImageType::DirectionType direction;
    direction.SetIdentity();

    /* Create simple test images. */
    ImageType::Pointer fixedImage = ImageType::New();
    fixedImage->SetRegions(region);
    fixedImage->SetSpacing(spacing);
    fixedImage->SetOrigin(origin);
    fixedImage->SetDirection(direction);
    fixedImage->Allocate();

    ImageType::Pointer movingImage = ImageType::New();
    movingImage->SetRegions(region);
    movingImage->SetSpacing(spacing);
    movingImage->SetOrigin(origin);
    movingImage->SetDirection(direction);
    movingImage->Allocate();

    /* Fill images */
    itk::ImageRegionIterator<ImageType> itFixed(fixedImage, region);
    itFixed.GoToBegin();
    itk::SizeValueType count = 1;
    while (!itFixed.IsAtEnd())
      {
      itFixed.Set(count * count);
      count++;
      ++itFixed;
      }
    itk::ImageRegionIteratorWithIndex<ImageType> itMoving(movingImage, region);
    itMoving.GoToBegin();
    count = 1;
    while (!itMoving.IsAtEnd())
      {
      itMoving.Set(count * count);
      count++;
      ++itMoving;
      }

    VectorType zero;
    float def_value = 2.5;

    def_value = -0.5;
    zero.Fill(def_value);
    FieldType::Pointer field = FieldType::New();
    field->SetRegions(fixedImage->GetLargestPossibleRegion());
    field->SetSpacing(fixedImage->GetSpacing());
    field->SetOrigin(fixedImage->GetOrigin());
    field->SetDirection(fixedImage->GetDirection());
    field->Allocate();
    field->FillBuffer(zero);

    FieldType::Pointer fieldInv = FieldType::New();

    zero.Fill(def_value * (-1.0));
    fieldInv->SetRegions(fixedImage->GetLargestPossibleRegion());
    fieldInv->SetSpacing(fixedImage->GetSpacing());
    fieldInv->SetOrigin(fixedImage->GetOrigin());
    fieldInv->SetDirection(fixedImage->GetDirection());
    fieldInv->Allocate();
    fieldInv->FillBuffer(zero);

    zero.Fill(def_value * (1.0));
    transformMtranslation->Translate(zero);
    zero.Fill(def_value * (1.0));
    transformMtranslation2->Translate(zero);

    transformMdisplacement->SetDisplacementField(field);
    transformMdisplacement->SetInverseDisplacementField(fieldInv);

    transformMComp->AddTransform(transformMtranslation);
    transformFComp->AddTransform(transformFId);

    typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<ImageType, ImageType> MetricType;

    typedef MetricType::Pointer MetricTypePointer;
    MetricTypePointer metric = MetricType::New();

    itk::Size<ImageDimension> neighborhood_radius;
    neighborhood_radius.Fill(1);

    metric->SetRadius(neighborhood_radius);

    metric->SetFixedImage(fixedImage);
    metric->SetMovingImage(movingImage);

    metric->SetFixedTransform(transformFId);
    metric->SetMovingTransform(transformMdisplacement);

    std::cout << "fixedImage:" << std::endl;
    ANTSNeighborhoodCorrelationImageToImageMetricv4Test_PrintImage(fixedImage);

    std::cout << "movingImage:" << std::endl;
    ANTSNeighborhoodCorrelationImageToImageMetricv4Test_PrintImage(movingImage);

    /* Initialize. */
    try
      {
      std::cout << "Calling Initialize..." << std::endl;
      metric->Initialize();
      }
    catch (itk::ExceptionObject & exc)
      {
      std::cerr << "Caught unexpected exception during Initialize: " << exc;
      std::cerr << "Test FAILED." << std::endl;
      return EXIT_FAILURE;
      }

    // Evaluate
    MetricType::MeasureType valueReturn1;
    MetricType::DerivativeType derivativeReturn;
    try
      {
      std::cout << "Calling GetValueAndDerivative..." << std::endl;
      metric->GetValueAndDerivative(valueReturn1, derivativeReturn);
      }
    catch (itk::ExceptionObject & exc)
      {
      std::cerr << "Caught unexpected exception during GetValueAndDerivative: "
                << exc;
      std::cerr << "Test FAILED." << std::endl;
      return EXIT_FAILURE;
      }

    /* Re-initialize. */
    try
      {
      std::cout << "Calling Initialize..." << std::endl;
      metric->Initialize();
      }
    catch (itk::ExceptionObject & exc)
      {
      std::cerr << "Caught unexpected exception during re-initialize: " << exc;
      std::cerr << "Test FAILED." << std::endl;
      return EXIT_FAILURE;
      }

    // Evaluate with GetValue
    MetricType::MeasureType valueReturn2;
    try
      {
      std::cout << "Calling GetValue..." << std::endl;
      valueReturn2 = metric->GetValue();
      }
    catch (itk::ExceptionObject & exc)
      {
      std::cerr << "Caught unexpected exception during GetValue: "
                << exc;
      std::cerr << "Test FAILED." << std::endl;
      return EXIT_FAILURE;
      }

  // Test same value returned by different methods
  std::cout << "Check Value return values..." << std::endl;
  if( valueReturn1 != valueReturn2 )
    {
    std::cerr << "Results for Value don't match: " << valueReturn1
              << ", " << valueReturn2 << std::endl;
    }

  std::cout << "Test passed." << std::endl;
  std::cout << "transformMdisplacement parameters" << std::endl;
  std::cout << transformMdisplacement->GetParameters() << std::endl;
  ANTSNeighborhoodCorrelationImageToImageMetricv4Test_PrintImage(transformMdisplacement->GetDisplacementField());

  std::cout << "derivative of moving transform:" << std::endl;
  std::cout << derivativeReturn << std::endl;
  std::cout << std::endl << "derivative of moving transform as a field:" << std::endl;
  ANTSNeighborhoodCorrelationImageToImageMetricv4Test_PrintDerivativeAsVectorImage(fixedImage, derivativeReturn, ImageDimension);

  // Test that non-overlapping images will generate a warning
  // and return max value for metric value.
  DisplacementTransformType::ParametersType parameters( transformMdisplacement->GetNumberOfParameters() );
  parameters.Fill( static_cast<DisplacementTransformType::ParametersValueType>(1000.0) );
  transformMdisplacement->SetParameters( parameters );
  MetricType::MeasureType expectedMetricMax, valueReturn;
  expectedMetricMax = itk::NumericTraits<MetricType::MeasureType>::max();
  std::cout << "Testing non-overlapping images. Expect a warning:" << std::endl;
  metric->GetValueAndDerivative( valueReturn, derivativeReturn );
  if( metric->GetNumberOfValidPoints() != 0 || valueReturn != expectedMetricMax )
    {
    std::cerr << "Failed testing for non-overlapping images. " << std::endl
              << "  Number of valid points: " << metric->GetNumberOfValidPoints() << std::endl
              << "  Metric value: " << valueReturn << std::endl
              << "  Expected metric max value: " << expectedMetricMax << std::endl;
    }

  // Test that an exception is thrown if sampling is set
  std::cout << "Testing use with sampled point set, expect a warning: " << std::endl;
  metric->UseFixedSampledPointSetOn();
  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  std::cout << "Test PASSED." << std::endl;
  return EXIT_SUCCESS;

}
